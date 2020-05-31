/*********************************************
* integrator.h
* - simulates O3 model 
* - using cluster
* - observables will be used in spectral recon
* 
************************************************/

#ifndef __integrator_h__  // Used to ensure only one instance is running
#define __integrator_h__

#include <iostream>
#include <stack>
#include <array>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <complex>

class clusterAlg{

    const int L; const int T; const int d = 3; // lattice size and dimension
    const double beta, p0; //lattice parameters

    //-----------Cluster parameters ----

    int cluster_size_r = 0; 
    int cluster_size_u = 0; 

    std::vector<std::vector<std::vector<double>>> sigmas; // lattice
    
    std::vector<std::vector<bool>> cluster_r;  // clusters
    std::vector<std::vector<bool>> cluster_u;
    
    std::stack<std::vector<int>> neighbors;   

    std::vector<std::vector<std::complex<double>>> avg_r;
    std::vector<std::vector<std::complex<double>>> avg_u;
    
    std::vector<std::vector<double>> avg_all;

    std::vector<double> r;
    std::vector<double> u;

    // Random number generators
    std::mt19937 gen;  
    std::uniform_real_distribution<double> dist;
    std::uniform_int_distribution<int> Ldist;
    std::uniform_int_distribution<int> Tdist;


    // method for dot product
    double dot(const std::vector<double> &v1, const std::vector<double> &v2){
        double tmp = 0;
        for (int i = 0; i < d; i++)
            tmp += v1[i]*v2[i];
        
        return tmp;
    }

    // method for complex dot product
    std::complex<double> dot_complex(const std::vector<std::complex<double>> &v1, const std::vector<double> &v2){
        std::complex<double> tmp = 0;
        for (int i = 0; i < d; i++)
            tmp += v1[i]*v2[i];
        
        return tmp;
    }

    // method for adding vectors
    std::vector<double> add(std::vector<double>& v1, std::vector<double>& v2) {
        std::vector<double> tmp(d,0);
        for (int i = 0; i < d; i++) 
            tmp[i] = v1[i] + v2[i]; 
        
        return tmp;
    }

    // calculates the spatial averages in a cluster 
    void calcAve(std::vector<std::vector<bool>>& cluster, std::vector<std::vector<std::complex<double>>>& avg){
        avg = std::vector<std::vector<std::complex<double>>>(T,std::vector<std::complex<double>>(d,0));
        
        for (int i = 0; i < T; i++) 
            for (int j = 0; j < L; j++)
                if (cluster[j][i]){
                    // Adding a momentum projection or fourie mode of the spins, by multiplying by exp(ia)
                    double a = p0*(double)j;
                    avg[i][0] += sigmas[j][i][0]*std::complex<double>(cos(a),sin(a));
                    avg[i][1] += sigmas[j][i][1]*std::complex<double>(cos(a),sin(a));
                    avg[i][2] += sigmas[j][i][2]*std::complex<double>(cos(a),sin(a));
                }
    }

    // calculates the spatial average in the whole lattice
    void calcAve_all(std::vector<std::vector<double>>& avg){
        avg = std::vector<std::vector<double>>(T,std::vector<double>(d,0));
        
        for (int i = 0; i < T; i++) 
            for (int j = 0; j < L; j++)
                avg[i] = add(avg[i], sigmas[j][i]);
            
    }        

    // produces a uniformly distributed unit vecter in 3 dimensions
    void unit_vec(std::vector<double>& v){
        std::uniform_real_distribution<double> phi_dist(0,2.0*M_PI);
        std::uniform_real_distribution<double> cos_theta_dist(-1,1);

        double phi=phi_dist(gen);
        double cos_theta=cos_theta_dist(gen);
        double sin_theta=sqrt(1 - cos_theta*cos_theta);    
        
        v[0] = sin_theta*cos(phi);
        v[1] = sin_theta*sin(phi);
        v[2] = cos_theta;
    }

    // produces an orthogonal vector
    void orthogonal_vec(const std::vector<double>& v1, std::vector<double>& v2){
        std::vector<double> tmp(d, 0);
        tmp[0] = v1[1]*v2[2]-v1[2]*v2[1];
        tmp[1] = v1[2]*v2[0]-v1[0]*v2[2];        
        tmp[2] = v1[0]*v2[1]-v1[1]*v2[0];
        
        double norm = std::sqrt(dot(tmp, tmp));
        tmp[0] /= norm;
        tmp[1] /= norm;
        tmp[2] /= norm;
        
        v2 = tmp;
    }

    public:
        clusterAlg(const int &T_in, const int &L_in, const int &seed, const double &beta_in, const double &p0_in)
            : T(T_in), L(L_in), beta(beta_in), p0(p0_in)
        {
            std::cerr << beta << " " << T << " " << L << " " << d << std::endl;
            sigmas = std::vector<std::vector<std::vector<double>>>(L, std::vector<std::vector<double>>(T, std::vector<double>(d, 0)));
            cluster_r = std::vector<std::vector<bool>>(L, std::vector<bool>(T,false));
            cluster_u = std::vector<std::vector<bool>>(L, std::vector<bool>(T,false));

            r = std::vector<double>(d,0);
            u = std::vector<double>(d,0);
            
            // Simple method for filling the lattice with unit vectors 
            for (int i = 0; i < L; i++)
                for (int j = 0; j < T; j++)
                    sigmas[i][j][0] = 1;

            gen = std::mt19937(seed);
            dist = std::uniform_real_distribution<double>(-1.0, 1.0);
            Ldist = std::uniform_int_distribution<int>(0, L-1);
            Tdist = std::uniform_int_distribution<int>(0, T-1);
        }
        
        // calculates the number of spins in the cluster
        int c_size(const int& cluster_count){
            int tmp = 0;
            if (cluster_count == 1)
                tmp += cluster_size_r;
            else if (cluster_count == 2) 
                tmp += cluster_size_r + cluster_size_u;
            return tmp;
        }
        
        // updates the lattice according to the cluster algorithm
        void update(const int& cluster_count){
            cluster_r = std::vector<std::vector<bool>>(L, std::vector<bool>(T,false));
            unit_vec(r);
            
            // ------------------ first find the cluster C_r -----------

            getCluster(cluster_r, r, cluster_size_r);
            
            // and calculate the spatial averages in th cluster
            calcAve(cluster_r, avg_r);
            // if a second cluster is used, do the same for that
            if (cluster_count == 2){
                cluster_u = std::vector<std::vector<bool>>(L, std::vector<bool>(T,false));
                unit_vec(u);
                orthogonal_vec(r, u);
                getCluster(cluster_u, u, cluster_size_u);
                calcAve(cluster_u, avg_u);
                calcAve_all(avg_all);
                //update the clusters
                updateCluster(cluster_u, u);
            }
            updateCluster(cluster_r, r);
        }

        void getCluster(std::vector<std::vector<bool>>& cluster, std::vector<double>& reflec, int& size) {
            //First cluster with random uniform cluster
            std::stack<std::array<int,2>> cluster_buffer;

            //choose random site
            int i_start = Ldist(gen);
            int j_start = Tdist(gen);
            
            cluster[i_start][j_start] = true;
            size = 1;        

            cluster_buffer.push({i_start, j_start}); //list of unchecked spins.

            while (cluster_buffer.size() > 0){
                const auto& cur_site = cluster_buffer.top(); //save the top element

                int i = cur_site[0];
                int j = cur_site[1];
                
                cluster_buffer.pop(); //remove top element

                double sxr = dot(sigmas[i][j], reflec);

                //add neighbours
                std::vector<std::array<int, 2> > neighbours;
                
                int im = (i-1+L)%L; //up 
                int ip = (i+1+L)%L; //down 
                int jm = (j-1+T)%T; //left
                int jp = (j+1+T)%T; //right
                
                neighbours.push_back({im, j });
                neighbours.push_back({ip, j });
                neighbours.push_back({i , jm});
                neighbours.push_back({i , jp});

                //check for each neighbour
                for (const auto& n: neighbours) {
                    int ni = n[0];
                    int nj = n[1];

                    if (cluster[ni][nj])
                        continue;

                    double syr = dot(sigmas[ni][nj], reflec);
                    //check if it should be added to the cluster
                    if(sxr*syr > 0) {
                            
                        if (std::abs(dist(gen)) < (1 - std::exp(-2*beta*sxr*syr))){
                            
                            size++;
                            cluster_buffer.push({ni, nj}); //Add to unchecked list if accepted
                            cluster[ni][nj] = true; //add checked spin to cluster
                        }
                    }
                }
            } 
        }

        void updateCluster(std::vector<std::vector<bool>>& cluster, std::vector<double>& reflec){ //Updating the spins in the referenced cluster, using the reflection
            double ortho = 0;
            for (int i = 0; i < L; i++){
                for (int j = 0; j < T; j++) {
                    if(cluster[i][j]){
                        ortho = 2 * dot(sigmas[i][j], reflec);
                        for (int h = 0; h < d; h++)
                            sigmas[i][j][h] -= ortho * reflec[h];
                    }
                }
            }
        }

        void printCluster(){ //Used to visualize the cluster in the lattice. Done by 1 if spin is in cluster and 0 if not.
            for(int i = 0; i < L; i++){
                for (int j = 0; j < T; j++){
                    if (cluster_u[i][j])
                        std::cout << 1 << " ";
                    else
                        std::cout << 0 << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        // calculate the 2point correlator
        double two_point_obs(const int& t2, const int& t1) {
            std::complex<double> tmp = 0;
            for (int j = 0; j < T; j++) {
                tmp += dot_complex(avg_r[(j+t2+T)%T],r) * dot_complex(avg_r[(j+t1+T)%T],r)/(double)cluster_size_r;
                tmp += dot_complex(avg_u[(j+t2+T)%T],u) * dot_complex(avg_u[(j+t1+T)%T],u)/(double)cluster_size_u;
            }
            tmp *= 1.5;
            return tmp.real(); // Only the real value is used.
        }

};  

#endif 