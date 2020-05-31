#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>

class analysis
{
    std::mt19937 gen;
    std::uniform_int_distribution<int> dist;

    //Used to split data
    size_t split(const std::string &txt, std::vector<std::string> &strs, char ch)
    {
        size_t pos = txt.find(ch);
        size_t initialPos = 0;
        strs.clear();

        // Decompose statement
        while (pos != std::string::npos)
        {
            strs.push_back(txt.substr(initialPos, pos - initialPos));
            initialPos = pos + 1;

            pos = txt.find(ch, initialPos);
        }

        // Add the last one
        strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

        return strs.size();
    }

public:
    //parameters from data
    int tmin, tmax, smpls, rows;

    //Initialize random generators, used in bootstrap
    void init_rands(const long int& seed){
        gen = std::mt19937(seed);
        dist = std::uniform_int_distribution<int>(0, smpls-1);
    }

    // read data files
    void readFile(std::string &path, std::vector<std::vector<double>> &data) {
        std::string line;
        std::ifstream file(path+".dat");
        int i = -1; int j = 0;
        if (file.is_open())
        {
            while (getline(file, line)) {
                std::vector<std::string> vec;

                split(line, vec, ' ');
                // read first line
                if (i == -1) { 
                    smpls = stoi(vec[0]);
                    tmin = stoi(vec[1]);
                    tmax = stoi(vec[2]);
                    rows = tmax - tmin;
                    data = std::vector<std::vector<double>>(rows, std::vector<double>(smpls, 0));
                    i++;
                
                // read the rest
                } else {
                    data[i][j] = stod(vec[1]); //store the second column

                    if ((i+1) % (rows) == 0){    
                        j++;
                        i = 0;
                    } else {
                        i++;
                    }
                    if (j == smpls) {
                        break;
                    }
                    
                }
            }
        }
    }

    // Calculate the covariance matrix
    void covariance(std::vector<double> &mean, std::vector<std::vector<double>> &mat, std::vector<std::vector<double>> &data, const int cols)
    {
        mat = std::vector<std::vector<double>>(cols, std::vector<double>(cols, 0));
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < cols; j++) {
                for (int k = 0; k < smpls; k++) {
                    mat[i][j] += (data[i][k]  - mean[i]) * (data[j][k] - mean[j]);
                }
                mat[i][j] /= (double)(smpls);
            }
        }
    }

    //Calculate the average
    void average(double &mean, const std::vector<double> &data)
    {
        mean = 0;
        for (int i = 0; i < data.size(); i++)
            mean += data[i];
        
        mean /= data.size();
    }

    // Make m bootstrap samples.
    void bootstrap(const std::vector<double>& data, std::vector<double>& out, const int& m){
        std::vector<double> tmp_data(smpls, 0); //store selected values
        out = std::vector<double> (m, 0); //store final values
        int selected = 0;
        
        
        for (int n = 0; n < m; n++) {
            for (int i = 0; i < data.size(); i++){
                selected = dist(gen);
                tmp_data[i] = data[selected];
            }

            average(out[n], tmp_data);
        }
        //sort the sample
        std::sort(out.begin(), out.end());
    }
};
