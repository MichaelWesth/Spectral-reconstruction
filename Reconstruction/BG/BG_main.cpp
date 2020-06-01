#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>

#include "bg_recon_arb_mwh.h"
#include "../read_dats.h"

using namespace std;

double lambda, E0, energy; //0.039393815930658;
int ndigits = 1024, tmin, tmax, rows, T;

vector<vector<double>> mat;
vector<vector<double>> cov_mat;
vector<vector<double>> observ;
vector<vector<double>> cut_data;

vector<double> ave;
vector<double> error;
vector<double> qvec;

string path = "level_00.dat";

void getRho(Chimera::BGRecon &bg, analysis &read, const double &min, const double &max, const int &num) {
    for (int i = 0; i <= num; i++){
        double om = min + i / ((double)(num)) * (max - min);
        bg.setE(om);
        
        vector<double> rho_dat(observ.size(), 0);
        for (int j = 0; j < observ.size(); j++){
            rho_dat[j] = bg.getRho(observ[j]);
        }

        vector<double> rho_boot(1000, 0);
        read.bootstrap(rho_dat, rho_boot, 1000);
    
        double err = (rho_boot[840] - rho_boot[160])/2;
        double rho;
        read.average(rho, rho_dat);
        cout << om << " " << rho << " " << err << endl;
    }
}

int main(int argc, char** argv)
{
    string test, kern;

    double min, max, num;

    if (argc != 13)
    {
        cerr << "usage: " << argv[0] << " <lambda> <E0> <E> <tmin> <tmax> <T> <test> <observ> <min> <max> <num>" << endl;
        return 1;
    }

    // Defines the parameters
    stringstream(argv[1])  >> lambda;
    stringstream(argv[2])  >> E0;
    stringstream(argv[3])  >> energy;
    stringstream(argv[4])  >> tmin;
    stringstream(argv[5])  >> tmax;
    stringstream(argv[6])  >> T;
    stringstream(argv[7])  >> test;
    stringstream(argv[8])  >> path;
    stringstream(argv[9])  >> kern;
    stringstream(argv[10])  >> min;
    stringstream(argv[11]) >> max;
    stringstream(argv[12]) >> num;    

    int seed = 3454980158;

    rows = tmax - tmin;
    mat = vector<vector<double>>(rows, vector<double>(rows, 0));
    cov_mat = vector<vector<double>>(rows, vector<double>(rows, 0));
    observ = vector<vector<double>>(rows, vector<double>(rows, 0));
    ave = vector<double>(rows, 0);
    error = vector<double>(rows, 0);
    qvec = vector<double>(rows, 0);
    double scale = 1;

    //if lambda is non-zero, calculate average and covariance matrix
    analysis read;  
    if (lambda != 0) {      
        read.readFile(path, observ);
        cut_data = vector<vector<double>>(rows, vector<double>(read.smpls, 0));
        for (int i = 0; i < rows; i++) 
            for (int j = 0; j < read.smpls; j++) 
                cut_data[i][j] = observ[i][j+tmin];
        
        observ = vector<vector<double>>(read.smpls, vector<double>(rows, 0));
        for (int i = 0; i < rows; i++) 
            for (int j = 0; j < read.smpls; j++) 
                observ[j][i] = cut_data[i][j];

        for(int i = 0; i < rows; i++) {
            read.average(ave[i], cut_data[i]);
            //cerr << ave[i] << endl;
        }

        read.covariance(ave, cov_mat, cut_data, rows);
        
        read.init_rands(seed);
        
        scale = cov_mat[0][0];
    }
    //print out the size for correctness
    cerr << ave.size() << " " << cov_mat.size() << endl;

    Chimera::BGRecon bg(lambda, E0, ndigits, tmin, tmax, T, cov_mat, kern, "no", scale);

    //other functionals can be added, here only rho is used.
    if (test == "rho") {
        getRho(bg, read, min, max, num);
    }
     
    return 0;
}
