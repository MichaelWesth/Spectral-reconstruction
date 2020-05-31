#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include "integrator.h"
//#include "observables.h"

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc != 8) {
        cerr << "usage: " << argv[0] << " <T> <L> <nsim> <nmeas> <beta> <test> <seed> " << endl;
        return 1;
    }

    // --------------------------------Store the input parameters------------------------------
    int T;       stringstream(argv[1]) >> T;        //Size
    int L;       stringstream(argv[2]) >> L;        //Size
    int nsim;    stringstream(argv[3]) >> nsim;     //Number of updates, after thermalization
    int nmeas;   stringstream(argv[4]) >> nmeas;    //Number of observations averaged for a measurement
    double beta; stringstream(argv[5]) >> beta;     //Coupling constant
    string test; stringstream(argv[6]) >> test;     //Used for different tests, Only one implemented here
    int seed;    stringstream(argv[7]) >> seed;     //USed to replicate data

    int cluster_size;

    // --------------------------------Initialize the algorithm------------------------------
    clusterAlg alg(T, L, seed, beta, 0);
    cerr << T << " " << L << " " << 3 << endl;
    int T_half = T/2.0;
    
    // Time the algorithm for fun
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    
    cerr << test << endl;
    

    //This is for 2point correlation
    vector<double> tw_pt(T_half, 0);
    vector<double> dat(2,0);

    double avg_c_size = 0;
    string path_tw_pt = "2_pt_obs.dat";
    ofstream out_tw_pt(path_tw_pt);

    out_tw_pt << nsim/nmeas << " " << 0 << " " << T_half << endl;
    // -------------------Run simulation 10% extra to reach a truly random configuration--------
    for (int i = 0; i <= nsim*1.1; i++){
        // update the lattice with 2 clusters
        alg.update(2); 

        if (i > nsim*0.1) {
            avg_c_size += alg.c_size(2);
            
            // Get the observation and average it
            for (int j = 0; j < T_half; j++) 
                tw_pt[j] += alg.two_point_obs(j, 0);

            // Store the average for eac nmeas'th observation
            if (i % nmeas == 0) {
                
                avg_c_size = 0;
                for (int j = 0; j < T_half; j++) {
                    tw_pt[j] /= (double)nmeas; 
                    out_tw_pt << tw_pt[j] << endl;
                    tw_pt[j] = 0;
                }
            }
        }
    }
    out_tw_pt.close();

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    cerr << "time for " << nsim << " simulations: " << time << endl;
    
    return 0;
}
