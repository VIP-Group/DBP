# Simulator for Decentralized Baseband Processing (DBP)
(c) 2017 Christoph Studer 
e-mail: studer@cornell.edu 

More information about our research can be found at [http://vip.ece.cornell.edu].

### Important information 

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

K. Li, R. Sharan, Y. Chen, T. Goldstein, J. R. Cavallaro, and C. Studer, "Decentralized Baseband Processing for Massive MU-MIMO Systems", IEEE J. Emerging and Sel. Topics in Circuits and Systems (JETCAS), to appear in 2017

and clearly mention this in your paper.  

### How to start a simulation:

For a downlink simulation, go to the downlink/ folder and simply run

```sh
DBP_precoder_sim
```

which starts a simulation in a massive MIMO system with 8 users and 128 base station antennas using 16-QAM. The simulator runs with predefined parameters. You can specify your own system, algorithm, and simulation parameters by passing your own "par" structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with different parameters, then please refer to the MATLAB code for other parameter settings.

For an uplink simulation, go to the uplink/ folder and simply run

```sh
DBP_detector_sim
```
 	
We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
* Version 0.2 (Jul 24 2018) - studer@cornell.edu  - improved README.md file
