/*********************************************
* hlt_recon_arb.h
* - uses the arb implementation 
* - Assumes e^{-Et} kernel (maybe cosh later)
* - only 1-d spectral recon for now
* 
************************************************/

#ifndef __bg_recon_arb_h__
#define __bg_recon_arb_h__

#include <iostream>
#include <vector>
#include <string>
#include <acb.h>
#include <acb_mat.h>
#include <acb_hypgeom.h>
#include <arb_hypgeom.h>

namespace Chimera {
    class BGRecon {

    // Define parameters
    private:
        arb_t lambda, E0, Ecur;
        int ndigits, tmin, tmax, T;
        arb_mat_t cov, M, R, q;
        std::string ker_type, constraint;
        double scale = 1;
// method to compute LU-decomposition.
void ludcmp(arb_mat_t& sol, arb_mat_t& mat, arb_mat_t& vec) {
    int ec = arb_mat_solve(sol, mat, vec, ndigits);
    if (ec == 0)
        throw std::string("matrix could not be inverted.");
}

// Defines the matrix M
void calcM() {

    arb_t ct0;
    arb_init(ct0);
    arb_set_d(ct0, scale);

    arb_t one_m_lambda, ea, ta, tmp1, tmp2, sum, ediff;
    
    arb_init(ta);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(sum);
    arb_init(ediff);
    arb_init(one_m_lambda);
    arb_init(ea);
    
    arb_neg(ea, E0);
    arb_exp(ea, ea, ndigits);

    arb_one(one_m_lambda);
    arb_sub(one_m_lambda, one_m_lambda, lambda, ndigits);
    int nt = tmax - tmin;
    arb_mat_init(M, nt, nt);
    
    arb_sub(ediff, E0, Ecur, ndigits);
    for (int i = 0; i < nt; i++) {
        for (int j = 0; j < nt; j++) {
            //For general part
            arb_set_d(ta, i + j + 2*tmin);

            arb_pow(tmp1, ea, ta, ndigits);
            arb_set_d(tmp2, 2);
            
            arb_div(sum, tmp2, ta, ndigits);

            arb_addmul(sum, ediff, tmp2, ndigits);
            arb_div(sum, sum, ta, ndigits);

            arb_addmul(sum, ediff, ediff, ndigits);
            arb_div(sum, sum, ta, ndigits);

            arb_mul(sum, sum, tmp1, ndigits);
            arb_mul(sum, sum, one_m_lambda, ndigits);

            // For the extra terms in thermal kernel
            if (ker_type == "therm") {
                // ------------------ second term
                arb_set_d(ta, 2*T - i - j - 2*tmin);

                arb_set_d(tmp2, 2);
                
                arb_div(tmp1, tmp2, ta, ndigits);
                arb_addmul(tmp1, ediff, tmp2, ndigits);
                arb_div(tmp1, tmp1, ta, ndigits);

                arb_addmul(tmp1, ediff, ediff, ndigits);
                arb_div(tmp1, tmp1, ta, ndigits);

                arb_pow(tmp2, ea, ta, ndigits);

                arb_mul(tmp1, tmp1, tmp2, ndigits);
                arb_mul(tmp1, tmp1, one_m_lambda, ndigits);

                arb_add(sum, sum, tmp1, ndigits);
                // ---------------- third term
                arb_set_d(ta, T + j - i);

                arb_set_d(tmp2, 2);
                
                arb_div(tmp1, tmp2, ta, ndigits);
                arb_addmul(tmp1, ediff, tmp2, ndigits);
                arb_div(tmp1, tmp1, ta, ndigits);

                arb_addmul(tmp1, ediff, ediff, ndigits);
                arb_div(tmp1, tmp1, ta, ndigits);

                arb_pow(tmp2, ea, ta, ndigits);

                arb_mul(tmp1, tmp1, tmp2, ndigits);
                arb_mul(tmp1, tmp1, one_m_lambda, ndigits);

                arb_add(sum, sum, tmp1, ndigits);

                // ---------------- fourth term
                arb_set_d(ta, T + i - j);

                arb_set_d(tmp2, 2);
                
                arb_div(tmp1, tmp2, ta, ndigits);
                arb_addmul(tmp1, ediff, tmp2, ndigits);
                arb_div(tmp1, tmp1, ta, ndigits);

                arb_addmul(tmp1, ediff, ediff, ndigits);
                arb_div(tmp1, tmp1, ta, ndigits);

                arb_pow(tmp2, ea, ta, ndigits);

                arb_mul(tmp1, tmp1, tmp2, ndigits);
                arb_mul(tmp1, tmp1, one_m_lambda, ndigits);

                arb_add(sum, sum, tmp1, ndigits);
            }
            // add the covariance matrix
            arb_div(tmp1,arb_mat_entry(cov,i,j), ct0,ndigits);
            
            arb_addmul(sum, lambda, tmp1, ndigits);
            arb_set(arb_mat_entry(M, i, j), sum);
        }
    }

    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(sum);
    arb_clear(ediff);
    arb_clear(ta);
    arb_clear(ea);
    arb_clear(one_m_lambda);
}

// calculate the R-vector
void calcR(){
    arb_t re1,re2,ta, ea;

    arb_init(re1);
    arb_init(re2);
    arb_init(ta);
    arb_init(ea);

    arb_neg(ea, E0);
    arb_exp(ea, ea, ndigits);

    int n = tmax-tmin;
    for (int i=0; i<n; i++) {
        arb_set_d(ta, i+tmin);
        arb_pow(re1, ea, ta, ndigits);

        arb_div(re1, re1, ta, ndigits);

        if (ker_type == "therm"){
            arb_set_d(ta, T - i - tmin);
            arb_pow(re2, ea, ta, ndigits);

            arb_div(re2, re2, ta, ndigits);
            arb_add(re1, re1, re2, ndigits);
        }

        arb_set(arb_mat_entry(R,i,0), re1);
    }

    arb_clear(re1);
    arb_clear(re2);
    arb_clear(ta);
    arb_clear(ea);
}

   public:
        BGRecon(const double& lambda_in, const double& E_0_in, int &ndigits_in, 
	int &tmin_in, int &tmax_in, int &T_in, const std::vector<std::vector<double>> &cmat_in, 
	const std::string& k_in, const std::string& c_in, const double& n_in) : 
	ndigits(ndigits_in), tmin(tmin_in), tmax(tmax_in), T(T_in), ker_type(k_in), constraint(c_in), scale(n_in) {
            
            std::cerr << "nbits = "<< ndigits <<std::endl;
            
            // convert values from double to ARB
            arb_init(lambda);
            arb_init(E0);
            
            arb_set_d(lambda,lambda_in);
            arb_set_d(E0,E_0_in);
            
            int n = tmax-tmin;
            arb_mat_init(cov,n,n);
            for (int i=0;i<n;i++)
                for (int j=0;j<=i;j++) {
                    arb_set_d(arb_mat_entry(cov,i,j),cmat_in[i][j]);
                    arb_set_d(arb_mat_entry(cov,j,i),cmat_in[j][i]);
                }

            arb_init(Ecur);
            arb_mat_init(R,n,1);
            calcM();
            calcR();
            arb_mat_init(q,n,1);
        }

        ~BGRecon() {
            arb_clear(lambda);
            arb_clear(E0);
            arb_clear(Ecur);

            arb_mat_clear(cov);
            arb_mat_clear(M);
            arb_mat_clear(R);
            arb_mat_clear(q);
        }
        
        // Calculates the q-vector for a new energy
        void setE(const double& energy){
            arb_set_d(Ecur, energy);
            calcM();
            calcR();
            
            ludcmp(q, M, R);

            int nt = tmax - tmin;

            arb_t tmp;
            arb_init(tmp);

            arb_zero(tmp);
            for(int i = 0; i < nt; i++){
                arb_addmul(tmp, arb_mat_entry(q, i, 0), arb_mat_entry(R, i, 0), ndigits);
            }
            for(int i = 0; i < nt; i++){
                arb_div(arb_mat_entry(q, i, 0), arb_mat_entry(q, i, 0), tmp, ndigits);
            }

        }
        
        //Calculate rho
        double getRho(const std::vector<double> &c_t) {
            int nt = tmax - tmin;
          
            arb_t rhat,ct;
            arb_init(rhat);
            arb_init(ct);
            arb_zero(rhat);
            for (int i = 0; i < nt; i++) {
                arb_set_d(ct,c_t[i]);
                arb_addmul(rhat,arb_mat_entry(q,i,0),ct,ndigits);
            }
            double rho_hat = arf_get_d(arb_midref(rhat),ARF_RND_NEAR);
            arb_clear(rhat);
            arb_clear(ct);

            return rho_hat;
        }
        
        //calculates the smearing function
        double getDelta(double omega) {
            arb_t ta, tmp, temp1, sum, earg;
            arb_init(ta);
            arb_init(tmp);
            arb_init(temp1);
            arb_init(sum);
            arb_init(earg);
            
            int nt = tmax-tmin+1;
            //double Emin = arf_get_d(arb_midref(E0),ARF_RND_NEAR);
            //double epsE = (Emax - Emin) / (double)(numE-1);
            arb_set_d(tmp, omega);
            arb_zero(sum);
            for (int t = 0; t < nt; t++) {
                arb_set_ui(ta, t + tmin);
                arb_neg(temp1, ta);
                arb_mul(earg, temp1, tmp, ndigits);
                arb_exp(temp1, earg, ndigits);
                arb_addmul(sum, temp1, arb_mat_entry(q, t, 0), ndigits);
            }
            arb_clear(temp1);
            arb_clear(tmp);
            arb_clear(ta);
            arb_clear(earg);
            arb_clear(sum); 
            return arf_get_d(arb_midref(sum), ARF_RND_NEAR);

        }
    };        
}
#endif
