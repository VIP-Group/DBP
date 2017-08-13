% =========================================================================
% Decentralized DOWNLINK simulator for the paper
% "Decentralized Baseband Processing for Massive MU-MIMO Systems"
% -------------------------------------------------------------------------
% Revision history:
%
%   - aug-13-2017  v0.1   cs: simplified and commented code for GitHub
%
% -------------------------------------------------------------------------
% (c) 2017 Christoph Studer; e-mail: studer@cornell.edu
% -------------------------------------------------------------------------
% If you are using the simulator (or parts of it) for a publication, then
% you MUST cite our paper:
%
% K. Li, R. Sharan, Y. Chen, T. Goldstein, J. R. Cavallaro, and C. Studer,
% "Decentralized Baseband Processing for Massive MU-MIMO Systems",
% IEEE J. Emerging and Sel. Topics in Circuits and Systems (JETCAS)
% to appear in 2017
%
% and clearly mention this in your paper.
%=========================================================================

function res = DBP_precoder_sim(varargin)

% -- set up default/custom parameters

if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.runId = 0; % simulation ID (used to reproduce results)
    par.U = 8; % user antennas
    par.B = 128; % BS antennas
    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 1000; % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list = 0:4:20; % list of SNR [dB] values to be simulated
    % select precoder to be used 
    %   centralized   : `ZF`, `MRC'
    %   decentralized : 'DP' as described in Algorithm 3
    par.precoder = 'DP'; % select precoder scheme 'MRC', ''
    par.CHEST = 'on'; % channel estimation 'on' or 'off'
    par.plot = 'on'; % plot results 'on' or 'off'
    par.save = 'on'; % save results: 'on' or 'off'
    
    % parameters for DBP (see paper)
    par.C = 8; % number of clusters   
    par.vers = 'SxS1'; % inverse: 'UxU1' or 'SxS1' or 'UxU2' or 'SxS2'
    par.maxiter = 5; % maximum algorithm iterations
    par.rho = 1; % tuning parameter: regularizer (only for ADMM)
    par.gamma = 1; % tuning parameters: step size (only for ADMM)
    par.epsilon = 0.0; % precoder accuracy (0 = ZF precoding)
    
else
    
    disp('use custom simulation settings and parameters...')
    par = varargin{1}; % only argument is par structure
    
end

% -- initialization

% use runId random seed (enables reproducibility)
rng(par.runId);

% generate reasonable filename
par.simName = ['ERR_DL_' num2str(par.B) 'x' num2str(par.U) '_' par.mod '_' par.precoder '_rho' num2str(par.rho) '_gamma' num2str(par.gamma) ]; % simulation name (used for saving results)

% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
        
end

% extract average symbol energy
par.Es = mean(abs(par.symbols).^2);

% normalize so that transmit vector s has unit norm
par.symbols = par.symbols/sqrt(mean(abs(par.symbols).^2))/sqrt(par.U);

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% extract cluster size
par.S = par.B/par.C; % cluster size

% initialize result arrays (detector x SNR)
if ~strcmp(par.precoder,'DP')
    par.maxiter = 1;
    par.distr = 0;
else
    par.distr = 1;
end

res.par = par; % save parameter structure
res.VER = zeros(par.maxiter,length(par.SNRdB_list)); % vector error rate
res.SER = zeros(par.maxiter,length(par.SNRdB_list)); % symbol error rate
res.BER = zeros(par.maxiter,length(par.SNRdB_list)); % bit error rate

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.U,par.Q,par.trials);

% trials loop
tic
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).'; % create unit-norm transmit vector
    
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.U,1)+1i*randn(par.U,1));
    H = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B));
    NH = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B)); % used for CHEST
    
    % TX side processing
    switch (par.precoder)
        case 'ZF' % ZF beamforming
            x = pinv(H)*s;
        case 'MRC' % MRC beamforming
            x = H'*(s./sum(abs(H).^2,2));
        case 'DP' % decentralized precoding
            x = DP(par,H,s);
        otherwise,
            error('par.precoder type not defined.')
    end
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
        
        % iteration loop
        for l = 1:par.maxiter
            
            par.Ex = norm(x(:,l),2)^2;
            
            % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
            N0 = par.Ex*10^(-par.SNRdB_list(k)/10);
            
            % channel estimation (CHEST)
            switch (par.CHEST)
                case 'on'
                    Hest = H + sqrt(N0/par.Ex)*NH; % error happens in uplink
                otherwise % assume perfect CHEST
                    Hest = H;
            end
            
            % transmit over noiseless channel
            Hx = Hest*x(:,l);
            
            % transmit data over noisy channel
            y = Hx+sqrt(N0)*n;
            
            [~,idxhat] = min(abs(y*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
            bithat = par.bits(idxhat,:);
            
            % -- compute error and complexity metrics
            err = (idx~=idxhat);
            res.VER(l,k) = res.VER(l,k) + any(err);
            res.SER(l,k) = res.SER(l,k) + sum(err)/par.U;
            res.BER(l,k) = res.BER(l,k) + sum(sum(bits(:,:,t)~=bithat))/(par.U*par.Q);
            
        end % iteration loop
        
    end % SNR loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
end % trials loop

% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.time_elapsed = time_elapsed + toc;

% -- save final results (res structure)
if strcmp(par.save,'on')
    if ~exist('results','dir')
        mkdir results
    end
    save([ 'results' filesep par.simName ],'res');
end

% -- show results (generates fairly nice Matlab plot)
if strcmp(par.plot,'on')
    marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
    h=figure(1);
    plot_list = unique(round(logspace(0,log10(par.maxiter),7)));
    for d=1:length(plot_list)
        if d==1
            semilogy(par.SNRdB_list,res.SER(plot_list(d),:),marker_style{d},'LineWidth',2)
            hold on
        else
            semilogy(par.SNRdB_list,res.SER(plot_list(d),:),marker_style{d},'LineWidth',2)
        end
    end
    hold off
    grid on
    xlabel('average SNR per receive antenna [dB]','FontSize',12)
    ylabel('uncoded symbol error rate (SER)','FontSize',12)
    axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1])
    if par.distr
        l = cell(1,length(plot_list));
        for i = 1:length(plot_list)
            l{i} = ['Iteration ' num2str(plot_list(i))];
        end
        legend(l,'Fontsize',18)
    end
    set(gca,'FontSize',12)
end

end

%% decentralized precoder
function [x] = DP(par,H,s)

%  -- initialize
x = zeros(par.B,par.maxiter); % output, each column corresponds to one iteration
lambda_c = zeros(par.U,par.C);
x_c = zeros(par.S,par.C);
w_c = zeros(par.U,par.C);
Hx_c = zeros(par.U,par.C);
H_c = zeros(par.U,par.S,par.C);

% important for fast convergence (reasonable initial guess)
z_c = max(par.U/par.B,1/par.C)*s*ones(1,par.C);

% -- preprocessing
for c=1:par.C
    H_c(:,:,c) = H(:,par.S*(c-1)+1:par.S*c); % get the appropriate part of H
    switch par.vers
        case 'SxS1'
            Ainv(:,:,c) = (H_c(:,:,c)'*H_c(:,:,c) + (1/par.rho)*eye(par.S))\(H_c(:,:,c)'); % SxS inverse
        case 'SxS2'
            Ainv(:,:,c) = inv(H_c(:,:,c)'*H_c(:,:,c) + (1/par.rho)*eye(par.S)); % SxS inverse
        case 'UxU1'
            Binv(:,:,c) = H_c(:,:,c)'/(H_c(:,:,c)*H_c(:,:,c)' + (1/par.rho)*eye(par.U)); % UxU inverse
        case 'UxU2'
            Binv(:,:,c) = inv(H_c(:,:,c)*H_c(:,:,c)' + (1/par.rho)*eye(par.U)); % UxU inverse
        otherwise
            error('mode not defined')
    end
end

% -- start iteration
for l = 1:par.maxiter
    
    % cluster-wise equalization
    for c=1:par.C
        switch par.vers
            case 'SxS1'
                x_c(:,c) = Ainv(:,:,c)*(z_c(:,c) + lambda_c(:,c)); % SxS inverse
            case 'SxS2'
                x_c(:,c) = Ainv(:,:,c)*(H_c(:,:,c)'*(z_c(:,c) + lambda_c(:,c))); % SxS inverse
            case 'UxU1'
                x_c(:,c) = Binv(:,:,c)*(z_c(:,c) + lambda_c(:,c)); % UxU inverse
            case 'UxU2'
                x_c(:,c) = H_c(:,:,c)'*(Binv(:,:,c)*(z_c(:,c) + lambda_c(:,c))); % UxU inverse
            otherwise
                error('mode not defined')
        end
        Hx_c(:,c) = H_c(:,:,c)*x_c(:,c);
        w_c(:,c) = Hx_c(:,c)-lambda_c(:,c);
    end
    
    % consensus step
    w_avg = s-sum(w_c,2);
    w_norm = norm(w_avg,2);
    w_avg = (max(0,1-par.epsilon/w_norm)*1/par.C)*w_avg; % projection
    
    % cluster-wise update
    for c=1:par.C
        z_c(:,c) = w_c(:,c)+w_avg;
        lambda_c(:,c) = lambda_c(:,c) - par.gamma*(Hx_c(:,c)-z_c(:,c));
    end
    
    x(:,l) = x_c(:); % vectorize output
    
end

end
