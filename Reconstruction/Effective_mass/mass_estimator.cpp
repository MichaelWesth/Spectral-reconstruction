#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
#include "../read_dats.h"

using namespace std;

// ----------------------------------Newton method------------------
//The function used
double func(const double& x, const double& c, const double t) {
    double f = c*cosh(x*(t+1)) - cosh(x*t);
    return f;
}
// The derivative of the function
double func_prime(const double& x, const double& c, const double t) {
    double f_prime = c*sinh(x*(t+1))*(t+1) - sinh(x*t)*t;
    return f_prime;
}
//The Newton method, used to find roots.
double newton(double x_0, const double& c, const double& t, const double tol, const int max_iter){
    double x;
    int i = 0;
    while (abs(func(x_0, c, t)) > tol && i < max_iter){
        x = x_0 - func(x_0, c, t)/func_prime(x_0, c, t);
        x_0 = x;
        i++;
    }
    
    return x;
}

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " <path> <kern>" << endl;
        return 1;
    }
    // -------------------------------Defines parameters

    string path;    stringstream(argv[1]) >> path; //Path to data
    string kern;    stringstream(argv[2]) >> kern; //kernel type, (zero, therm or avg)
    vector<vector<double>> data;
    analysis read;


    int seed = 3454980158;
    read.readFile(path, data);
    read.init_rands(seed);
    
    

    int m = 1000;
    vector<double> avg(data.size(), 0);
    vector<vector<double>> m_eff1(data.size(), vector<double>(m, 0));
    vector<vector<double>> boot(data.size(), vector<double>(m, 0));
    
    // Outputs the average of the data, with statistical errors. 
    if (kern == "avg"){
        for (int i = 0; i < data.size(); i++) {
            read.bootstrap(data[i], boot[i], m);
            read.average(avg[i], boot[i]);
            double err = (boot[i][int(m*0.84)] - boot[i][int(m*0.16)])/2.0;
            cout << i << " " << avg[i] << " " << err << endl;
        }
    // Outputs the effective mass with statistical errors
    } else {
        for (int i = 0; i < data.size(); i++) {
            read.bootstrap(data[i], boot[i], m); //calculate m bootstrap samples.
        }    
        for (int i = 0; i < (data.size()-1); i++) {
            cerr << "now calculating " << i << endl;
            for (int j = 0; j < m; j++){
                double c = boot[i][j]/boot[i+1][j];
                if (kern == "therm") {
                    double t = i - 160 ;
                    m_eff1[i][j] = newton(1, c, t, 10e-12, 1000);
                } else if (kern == "zero")
                    m_eff1[i][j] = log(c);
            }
            sort(m_eff1[i].begin(), m_eff1[i].end());
            double err = (m_eff1[i][int(m*0.84)] - m_eff1[i][int(m*0.16)])/2.0;
            cout << i << " " << m_eff1[i][int(m/2)] << " " << err << endl;
        }
    }
    return 0;
}