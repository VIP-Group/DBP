=============================================================
Simulator for decentralized DOWNLINK as detailed in the paper 
"Decentralized Baseband Processing for Massive MU-MIMO Systems" 
-------------------------------------------------------------
(c) 2017 Christoph Studer; e-mail: studer@cornell.edu 
=============================================================

# Important information:

If you use this simulator or parts of it, then you MUST cite our paper: 

K. Li, R. Sharan, Y. Chen, T. Goldstein, J. R. Cavallaro, and C. Studer, "Decentralized Baseband Processing for Massive MU-MIMO Systems", IEEE J. Emerging and Sel. Topics in Circuits and Systems (JETCAS), to appear in 2017

and clearly mention this in your paper. 

More information about our research can be found at: http://vip.ece.cornell.edu

# How to start a simulation:

For decentralized DOWNLINK in massive MU-MIMO wireless systems:

Simply run  

>> DBP_precoder_sim

which starts a simulation in a 128 BS antenna, 8 user massive MU-MIMO system with 16-QAM modulation and C=8 clusters. The default algorithm is decentralized ZF precoding with 5 iterations as detailed in Algorithm 3 of the paper; we set the tuning parameters \rho=1 and \gamma=1 to achieve good performance in the 5th iteration. The simulator runs with a set of predefined simulation parameters. You can provide your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example).

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

The following precoders are supported by the simulator.
  - ZF : centralized zero-forcing (ZF) precoder
  - MRC : (de-)centralized MRC precoder
  - DP : decentralized ZF precoder using ADMM
  
Note that for the decentralized precoders, there are four modes defined by par.vers, which determine the mode of the matrix inversion as described in the paper. 


