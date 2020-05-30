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

    
    /*
    void LambdaTest(Chimera::BGRecon &hlt, double &n_in) {

        std::string path_out = "lambda.dat";
        std::ofstream out(path_out);

        double func = 0;
        hlt.setE(energy);
        lambda = 0;
        while (lambda < .01) {
            hlt.setLambda(lambda, n_in);
            func = hlt.getFunc();
            //std::cout << lambda << ", " << rho_hat << std::endl;
            out << lambda << " " << func << "\n";
            lambda += 10e-6;
        }

        double opt_lambda = hlt.getOptLambda();  
        std::cout << "opt_lambda: " << opt_lambda << std::endl;

        out.close();

    }
    
    void EnergyTest(Chimera::BGRecon &hlt, int Enum, double Emax, double Emin) {

        std::string path_out = "Energy.dat";
        std::ofstream out(path_out);

        double rho_hat = 0;
        hlt.setE(energy);

        for (int i = 0; i < Enum; i++) {
            double om = (double) (i)/((double) (Enum - 1))*(Emax - Emin)+Emin;
            
            //double tmp = 2*epsilon*(energy-om)/pow(pow(energy-om,2)+epsilon*epsilon,2);
            double tmp = epsilon/(pow(energy-om,2)+epsilon*epsilon);
            
            out << om << " " << hlt.getDelta(om)  << " " << tmp << "\n";

        }

        out.close();
    }

    void EpsTest(Chimera::BGRecon &hlt, double &n_in, int EpsNum, double epsMax, double epsMin) {

        std::string path_out = "Eps.dat";
        std::ofstream out(path_out);

        hlt.setE(energy);

        for (int i = 0; i < EpsNum; i++) {
            double om = (double) (i)/((double) (EpsNum - 1))*(epsMax-epsMin) + epsMin;
            double eps = hlt.getEpsEff(om, n_in);
            out << om << " " << sqrt(eps) << " " << hlt.getFunc() <<  "\n";

        }

        out.close();

    }

    void QTest(Chimera::BGRecon &hlt, std::vector<double> &vec, int Num, double Max, double Min) {

        std::string path_out = "Qs.dat";
        std::ofstream out(path_out);

        double om = 0;
        for (int i = 0; i < Num; i++) {
            
            if (Num > 1)
                om = (double) (i)/((double) (Num - 1))*(Max-Min) + Min;
            else 
                om = Min;
            //hlt.setEpsilon(om);
            //hlt.setE(om);
            hlt.getQ(vec); 
            for (int j = 0; j < tmax-tmin+1; j++) {
                double val = vec[j];
                //double val = vec[j]*exp(-energy*(j+tmin));
                
                out << om << " " << j + tmin << " " << val << "\n";
            } 
        }

        out.close();
    }

    void VTest(Chimera::BGRecon &hlt, int Num, double Max, double Min, int t) {

        std::string path_out2 = "X_Y.dat";
        std::ofstream out2(path_out2);

        std::string path_out = "Vs.dat";
        std::ofstream out1(path_out);

        hlt.setE(energy);
        double val;
        double eps = 0, om = 0;
        for (int i = 0; i < Num; i++) {
            om = (double) (i)/((double) (Num - 1))*(Max-Min) + Min;
            hlt.setE(om);
            for (int j = 0; j < Num; j++) {
                eps = (double) (j)/((double) (Num - 1))*(Max-Min) + Min;
                hlt.setEpsilon(eps);
                hlt.getV(val, t);
                if (j == Num - 1) 
                    out1 << val;
                else
                    out1 << val << " ";
                    
            } 
            out1 << "\n";
            out2 << om << " " << om << "\n";
        }

        out1.close();
        out2.close();
    }
    
    void getInvM(Chimera::BGRecon &hlt, std::vector<std::vector<double>> &mat) {

        std::string path_out = "invM.dat";
        std::string path_2 = "X_Y.dat";
        std::ofstream out(path_out);
        std::ofstream out2(path_2);
        hlt.getInvW(mat);
        int nt = tmax - tmin + 1; 
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nt; j++) {
                //double val = pow(-1,i+j)*log(abs(mat[i][j]));
                double val = mat[i][j];
                if (j == nt-1) 
                    out << val;
                else
                    out << val << " ";
                    
            } 
            out2 << i << " " << i << "\n";
            out << "\n";
        }

        out.close();
    }
    */
    
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
    cerr << ave.size() << " " << cov_mat.size() << endl;

    Chimera::BGRecon bg(lambda, E0, ndigits, tmin, tmax, T, cov_mat, kern, "no", scale);

    if (test == "rho") {
        getRho(bg, read, min, max, num);
    }
     
    return 0;
}
