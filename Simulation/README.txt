Consists of two scripts;

O3model.cpp is the main class, outputting the simulated data.
Integrator.h is the simulator, storing a lattice and containing various methods for the cluster algorithm and observations.

Can be compiled by "g++ O3model.cpp -o O3sim", and run by "./O3sim 64 64 10000 100 1.5 5468945689". Resulting in a simulation for L = T = 64, with 10000 updates and a measurement for every 100 updates. beta = 1.5 and a random seed in the end.   
