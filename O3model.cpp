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

    int T;       stringstream(argv[1]) >> T;
    int L;       stringstream(argv[2]) >> L;
    int nsim;    stringstream(argv[3]) >> nsim;
    int nmeas;   stringstream(argv[4]) >> nmeas;
    double beta; stringstream(argv[5]) >> beta;
    string test; stringstream(argv[6]) >> test;
    int seed;    stringstream(argv[7]) >> seed;

    int cluster_size;

    clusterAlg alg(T, L, seed, beta, 0);
    cerr << T << " " << L << " " << 3 << endl;
    int T_half = T/2.0;
    int T_frth = T/2.0;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    
    cerr << test << endl;
    
    if (test == "data") {
        //This is for 3point correlation
        vector<vector<double>> thr_pt(T_frth, vector<double>(T_frth, 0));
        //This is for 2point correlation
        vector<double> tw_pt(T_half, 0);
        vector<double> dat(2,0);

        double avg_c_size = 0;
        string path_tw_pt = "2_pt_obs.dat";
        string path_thr_pt = "3_pt_obs.dat";
        ofstream out_tw_pt(path_tw_pt);
        ofstream out_thr_pt(path_thr_pt);
    
        out_tw_pt << nsim/nmeas << " " << 0 << " " << T_half << endl;
        out_thr_pt << nsim/nmeas << " " << 0 << " " << T_frth << endl;
        for (int i = 0; i <= nsim*1.1; i++){
            alg.update(2);
            if (i > nsim*0.1) {
                avg_c_size += alg.c_size(2);
                
                for (int j = 0; j < T_half; j++) 
                    tw_pt[j] += alg.two_point_obs(j, 0);
                
                for (int j = 0; j < T_frth; j++) {
                    for (int k = 0; k < T_frth; k++) {    
                        dat = alg.three_point_obs(j+k, k, 0);
                        thr_pt[j][k] += dat[0] - dat[1];
                    }
                }
                if (i % nmeas == 0) {
                    cerr << "Average cluster size " << avg_c_size/(double)(nmeas) << endl;
                    avg_c_size = 0;
                    for (int j = 0; j < T_half; j++) {
                        tw_pt[j] /= (double)nmeas; 
                        out_tw_pt << tw_pt[j] << endl;
                        tw_pt[j] = 0;
                    }
                    for (int j = 0; j < T_frth; j++) {
                        for (int k = 0; k < T_frth; k++) {
                            thr_pt[j][k] /= (double)nmeas; 
                            out_thr_pt << beta*thr_pt[j][k] << endl;
                            thr_pt[j][k] = 0;
                        }
                    }
                }
            }
        }
        out_tw_pt.close();
        out_thr_pt.close();

    } else if (test == "g") {

        string path = "g_func.dat";
        ofstream out(path);
    
        out << nsim/nmeas << " " << 0 << " " << T_frth << endl;
        
        vector<double> pt_2(T_frth, 0);
        vector<double> pt_3(T_frth, 0);
        vector<double> dat(2,0);
        double avg_c_size = 0;

        for (int i = 0; i <= nsim*1.1; i++){
            alg.update(2);
            if (i > nsim*0.1) {
                avg_c_size += alg.c_size(1);
                
                for (int j = 0; j < T_frth; j++) {
                    dat = alg.three_point_obs(j+1, 0, -j);
                    
                    pt_2[j] += alg.two_point_obs(2*j+1, 0);
                    pt_3[j] += (dat[0] - dat[1]);
                }
                if (i % nmeas == 0) {
                    //alg.printCluster();
                    //cerr << "Average cluster size " << avg_c_size/(double)(nmeas) << endl;
                    avg_c_size = 0;
                    for (int j = 0; j < T_frth; j++) {
                        double tmp = -beta*pt_3[j]/(2.0*pt_2[j]);
                        out << tmp << endl;
                        pt_2[j] = 0;
                        pt_3[j] = 0;
                    }
                }
            }
        }
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    cerr << "time for " << nsim << " simulations: " << time << endl;
    
    return 0;
}
