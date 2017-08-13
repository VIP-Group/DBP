% =========================================================================
% Decentralized UPLINK simulator for the paper 
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

function res = DBP_detector_sim(varargin)

% -- set up default/custom parameters

if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.runID = 0; % simulation ID (used to reproduce results)
    par.B = 128; % receive antennas
    par.U = 8; % transmit antennas (set not larger than MR!)
    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 1000; % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list = -4:4:16; % list of SNR [dB] values to be simulated
    % select data detector to be used 
    %   centralized   : `ZF`, `uMMSE', 'SIMO'
    %   decentralized : CG-based `DUCG_ZF`, `DUCG_MMSE`
    %                   ADMM-based 'DZF', `DMMSE', 'DBOX'
    par.detector = 'DMMSE'; 
    par.vers = 'SxS2'; % inverse: 'UxU1', 'SxS1', 'SxS2' (only for ADMM)
    par.CHEST = 'on'; % channel estimation errors 'on' or 'off'
    par.plot = 'on'; % plot results? 'on' or 'off'
    par.save = 'on'; % save results? 'on' or 'off'    
    
    % parameters for DBP (see paper)
    par.C = 8; % number of clusters
    par.maxiter = 5; % maximum algorithm iterations (for CG and ADMM)
    par.rho = 7; % tuning parameter: regularizer (only for ADMM)
    par.gamma = 2; % tuning parameters: step size (only for ADMM)
    %
else
    
    disp('use custom simulation settings and parameters...')
    par = varargin{1}; % only argument is par structure
    
end

% -- initialization

% use runId random seed (enables reproducibility)
rng(par.runID);

% generate reasonable filename
par.simName = ['ERR_UL_' num2str(par.B) 'x' num2str(par.U) '_' par.mod '_' par.detector '_rho' num2str(par.rho) '_gamma' num2str(par.gamma) ]; % simulation name (used for saving results)

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

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% initialize result arrays (detector x SNR)
switch (par.detector)
    case {'DZF','DMMSE','DBOX','DUCG_MMSE','DUCG_ZF'}
        par.distr = 1;
    otherwise
        par.distr = 0;
        par.maxiter = 1;
end

res.par = par; % store param array
res.VER = zeros(par.maxiter,length(par.SNRdB_list)); % vector error rate
res.SER = zeros(par.maxiter,length(par.SNRdB_list)); % symbol error rate
res.BER = zeros(par.maxiter,length(par.SNRdB_list)); % bit error rate

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.U,par.Q,par.trials);

%initialize parameters for ADMM

% trials loop
disp('run simulation...')
tic
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.B,1)+1i*randn(par.B,1));
    H = sqrt(0.5)*(randn(par.B,par.U)+1i*randn(par.B,par.U));
    NH = sqrt(0.5)*(randn(par.B,par.U)+1i*randn(par.B,par.U)); % used for CHEST
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
        
        % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
        N0 = par.U*par.Es*10^(-par.SNRdB_list(k)/10);
        
        % transmit data over noisy channel
        y = x+sqrt(N0)*n;
        
        % model channel estimation errors (CHEST)
        switch (par.CHEST)
            case 'on'
                Hest = H + sqrt(N0/par.U/par.Es)*NH; % errors proportional to SNR
            otherwise % assume perfect CSI
                Hest = H;
        end
        
        % select algorithms
        switch (par.detector) 
            case 'ZF', % zero-forcing detection
                [idxhat,bithat] = ZF(par,Hest,y);
            case 'MMSE', % unbiased MMSE detector
                [idxhat,bithat] = MMSE(par,Hest,y,N0);
            case 'SIMO', % SIMO lower bound
                [idxhat,bithat] = SIMO(par,Hest,y,s); % also pass true signal
            case {'DUCG_MMSE','DUCG_ZF'} % conjuage gradients
                [idxhat,bithat] = DUCG(par,Hest,y,N0);
            case {'DZF','DMMSE','DBOX'} % ADMM based method
                [idxhat,bithat] = DU(par,Hest,y,N0);
            otherwise,
                error('par.detector type not defined.')
        end
        
        % -- compute error metrics
        for l = 1:size(idxhat,2)
            [VER,SER,BER] = getError(par,idxhat(:,l),bithat(:,:,l),idx,bits(:,:,t));
            res.VER(l,k) = res.VER(l,k) + VER;
            res.SER(l,k) = res.SER(l,k) + SER;
            res.BER(l,k) = res.BER(l,k) + BER;
        end
                
    end % SNR loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
    
end % trials loop

% -- normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.time_elapsed = time_elapsed + toc;

% -- save final results (par and res structure)
if strcmp(par.save,'on')
    if ~exist('results','dir')
        mkdir results
    end
    save([ 'results' filesep par.simName ],'res');
end

% -- show results (generates fairly nice Matlab plot)

if strcmp(par.plot,'on')
    marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
    plot_list = unique(round(logspace(0,log10(par.maxiter),7)));
    figure(1)
    for d=1:length(plot_list)
        if d==1
            semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
            hold on
        else
            semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
        end
    end
    hold off
    grid on
    xlabel('average SNR per receive antenna [dB]','FontSize',12)
    ylabel('uncoded bit error rate (BER)','FontSize',12)
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

% -- data detector functions

%% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y)
xhat = H\y;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% unbiased MMSE detector (MMSE)
function [idxhat,bithat] = MMSE(par,H,y,N0)
W = (H'*H+(N0/par.Es)*eye(par.U))\(H');
xhat = W*y;
G = real(diag(W*H));
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% SIMO lower bound
function [idxhat,bithat] = SIMO(par,Hest,y,sTrue)
y_tilde = y-Hest*sTrue;
% -- MMSE detection main loop
for n=1:par.U    
    % interference cancellation (with known data)
    y_SIMO = y_tilde + Hest(:,n)*sTrue(n);    
    % do optimal SIMO detection for interference-free system
    xhat(n,1) = Hest(:,n)'*y_SIMO/norm(Hest(:,n),2)^2;    
end
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% find nearest neighbors
function [idxhat,bithat] = getEstimate(par,xhat)
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% get error
function [VER,SER,BER] = getError(par,idxhat,bithat,idx,bits)
err = (idx~=idxhat);
VER = any(err);
SER = sum(err)/par.U;
BER = sum(sum(bits~=bithat))/(par.U*par.Q);
end

%% ADMM-based decentralized uplink detectors
function [idxhat,bithat] = DU(par,H,y,N0)

idxhat = zeros(par.U,par.maxiter);
bithat = zeros(par.U,par.Q,par.maxiter);
lambda = zeros(par.U,par.C);
z_c = zeros(par.U,par.C);
s = zeros(par.U,1);
S = par.B/par.C;

switch (par.detector)
    case 'DMMSE'
        reg = par.C+N0/(par.Es*par.rho);
    otherwise % Do not regularize for ZF and BOX
        reg = par.C;
end

% -- preprocessing
for c=1:par.C
    H_c(:,:,c) = H(S*(c-1)+1:S*c,:); %get the appropriate part of H
    y_c = y(S*(c-1)+1:S*c);
    switch (par.vers)
        case 'SxS1'
            Ainv(:,:,c) = inv(H_c(:,:,c)*H_c(:,:,c)'+eye(S)*par.rho);
            y_reg(:,c) = H_c(:,:,c)'*(Ainv(:,:,c)*y_c);
        case 'SxS2'
            Ainv(:,:,c) = inv(H_c(:,:,c)*H_c(:,:,c)'+eye(S)*par.rho);
            Gtmp = H_c(:,:,c)'*Ainv(:,:,c);
            y_reg(:,c) = Gtmp*y_c;
            G(:,:,c) = Gtmp*H_c(:,:,c);
        case 'UxU1'
            Binv(:,:,c) = inv(H_c(:,:,c)'*H_c(:,:,c)+eye(par.U)*par.rho);
            y_reg(:,c) = Binv(:,:,c)*(H_c(:,:,c)'*y_c);
        otherwise
            error('par.vers not defined')
    end
end

% -- detection loop
for l = 1:par.maxiter
    
    for c = 1:par.C %local minimization step
        switch (par.vers)
            case 'SxS1'
                z_c(:,c) = y_reg(:,c) + ( (s - lambda(:,c)) - H_c(:,:,c)'*(Ainv(:,:,c)*(H_c(:,:,c)*(s - lambda(:,c)))) );
            case 'SxS2'
                z_c(:,c) = y_reg(:,c) + ( (s - lambda(:,c)) - G(:,:,c)*(s - lambda(:,c))) ;
            case 'UxU1'
                z_c(:,c) = y_reg(:,c) + par.rho*(Binv(:,:,c)*(s - lambda(:,c)));
        end
    end
    
    switch (par.detector)
        case 'DBOX'
            s = sum((z_c + lambda),2)/(reg); %global averaging step
            s = projinf(par,s,max(real(par.symbols))); % experimental box regularizer
        otherwise
            s = sum((z_c + lambda),2)/(reg); %global averaging step
    end
    
    % update lagrange multiplier
    for c = 1:par.C
        lambda(:,c) = lambda(:,c) + par.gamma*(z_c(:,c) - s);
    end
    [idxhat(:,l),bithat(:,:,l)] = getEstimate(par,s);
end
end

% project onto alpha infinity-tilde-norm ball
function sproj = projinf(par,s,alpha)
sr = real(s);
idxr = abs(sr)>alpha;
sr(idxr) = sign(sr(idxr))*alpha;
si = imag(s);
idxi = abs(si)>alpha;
si(idxi) = sign(si(idxi))*alpha;
if strcmp(par.mod,'BPSK')
    sproj = sr;
else
    sproj = sr + 1i*si;
end

end % ADMM methods

%% decentralized uplink via decentralized CG
function [idxhat,bithat] = DUCG(par,H,y,N0)

% initialization
S = par.B/par.C;
H_c = zeros(S,par.U,par.C);
yMRC_c = zeros(par.U,par.C);

% distributed preprocessing (iteration 1)
for c=1:par.C
    H_c(:,:,c) = H(S*(c-1)+1:S*c,:); % get the appropriate part of H
    yMRC_c(:,c) = H_c(:,:,c)'*y(S*(c-1)+1:S*c); % compute local MRC
end

% centralized processing
r = sum(yMRC_c,2); % sum row-wise
p = r ;
rsold = r'*r;
x = zeros(par.U,1);
Ap_c = zeros(par.U,par.C);

% inner loop = equalization stage (iterations 2,3,...)
for l=1:par.maxiter
    
    % decentralized matrix processing
    for c=1:par.C
        Ap_c(:,c) = H_c(:,:,c)'*(H_c(:,:,c)*p);
    end
    
    % centralized processing (MMSE)
    switch  par.detector
        case 'DUCG_MMSE'
            Ap = sum(Ap_c,2) + (N0/par.Es)*p;
        case 'DUCG_ZF'
            Ap = sum(Ap_c,2);
    end
    
    % conventional CG updates
    alpha = rsold/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rsnew = r'*r;
    p = r + (rsnew/rsold)*p;
    rsold = rsnew;
    
    % get estimates
    [idxhat(:,l),bithat(:,:,l)] = getEstimate(par,x);
end

end % DUCG