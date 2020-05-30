#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>

#include "hlt_therm_mwh.h"
#include "../read_dats.h"

using namespace std;

double lambda, epsilon, E0, energy; //0.039393815930658;
int ndigits = 1024, tmin, tmax, T, rows;
int seed = 3454980158;

vector<vector<double>> mat;
vector<vector<double>> cov_mat;
vector<vector<double>> observ;
vector<vector<double>> cut_data;
vector<double> ave;
vector<double> error;
vector<double> qvec;

string path = "level_00.dat";

void getW(Chimera::HLTRecon &hlt, double min, double max, int num)
{
    
    for (int i = 0; i < num; i++){
        double om = min + i / ((double)(num - 1)) * (max - min);
        double sigma = 1/(sqrt(2*M_PI)*hlt.getW(om, om)) ;
        double W = hlt.getW(energy, om) ;
        cout << om << " " << W << " " << om << " " << sigma << endl;
    }
}

void getRho(Chimera::HLTRecon &hlt, analysis &read, const double &min, const double &max, const int &num) {
    int m = 1000;
    vector<vector<double>> En(2, vector<double>(m, 0));
    for (int i = 0; i <= num; i++){
        double om = min + i / ((double)(num)) * (max - min);
        hlt.setE(om);
        
        vector<double> rho_dat(observ.size(), 0);
        for (int j = 0; j < observ.size(); j++){
            rho_dat[j] = hlt.getRho(observ[j]);
        }

        vector<double> rho_boot(m, 0);
        read.bootstrap(rho_dat, rho_boot, m);
        /*
        for (int j = 0; j < m; j++){
            if (En[0][j] < rho_boot[j]){
                En[0][j] = rho_boot[j];
                En[1][j] = om;
            }
            cout << om << " " << rho_boot[j] << " " << 0 << endl;
        }
        */
        double err = (rho_boot[int(m*0.84)] - rho_boot[int(m*0.16)])/2;
        double rho;
        read.average(rho, rho_dat);
        cout << om << " " << rho << " " << err << endl;
    
    }
}

void getDelta(Chimera::HLTRecon &hlt, const double &min, const double &max, const int &num) {
    
    vector<double> delta(num+1,0);
    hlt.setE(energy);
    hlt.getDelta(min, max, num, delta);
    for (int i = 0; i <= num; i++){
        double om = min + i / ((double)(num)) * (max - min);
        double del_eps = epsilon/((energy - om)*(energy-om)+epsilon*epsilon); 
        cout << om << " " << delta[i] << " " << del_eps << endl;
    }
}

void getOptEps(Chimera::HLTRecon &hlt, double min, double max, int num)
{
    for (int i = 0; i < num; i++){
        double om = min + i / ((double)(num - 1)) * (max - min);
        hlt.setEps(om);
        double eps = hlt.getOptEps();
        cout << om << " " << sqrt(eps) << endl;
    }
}

int main(int argc, char** argv)
{
    string test, kernel;
    double min, max, num;

    if (argc != 14)
    {
        cerr << "usage: " << argv[0] << " <lambda> <eps> <E0> <E> <tmin> <tmax> <T> <test> <observ> <kernel> <min> <max> <num>" << endl;
        return 1;
    }

    stringstream(argv[1])  >> lambda;
    stringstream(argv[2])  >> epsilon;
    stringstream(argv[3])  >> E0;
    stringstream(argv[4])  >> energy;
    stringstream(argv[5])  >> tmin;
    stringstream(argv[6])  >> tmax;
    stringstream(argv[7])  >> T;
    stringstream(argv[8])  >> test;
    stringstream(argv[9])  >> path;
    stringstream(argv[10])  >> kernel;
    stringstream(argv[11])  >> min;
    stringstream(argv[12]) >> max;
    stringstream(argv[13]) >> num;    

    rows = tmax - tmin;
    mat = vector<vector<double>>(rows, vector<double>(rows, 0));
    cov_mat = vector<vector<double>>(rows, vector<double>(rows, 0));
    observ = vector<vector<double>>(rows, vector<double>(rows, 0));
    
    ave = vector<double>(rows, 0);
    error = vector<double>(rows, 0);
    qvec = vector<double>(rows, 0);
    double scale = 1;
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
        }

        read.covariance(ave, cov_mat, cut_data, rows);
        
        read.init_rands(seed);
        
        scale = cov_mat[0][0];
    }

    cerr << ave.size() << " " << cov_mat.size() << endl;

    Chimera::HLTRecon hlt(lambda, epsilon, E0, ndigits, tmin, tmax, T, cov_mat, kernel, "no", true, scale);
    
    if (test == "W")
        getW(hlt, min, max, num);
    else if (test == "rho")
        getRho(hlt, read, min, max, num);
    else if (test == "delta")
        getDelta(hlt, min, max, num);
    else if (test == "eps")
        getOptEps(hlt, min, max, num);
    return 0;
}
