/*********************************************
* hlt_recon_arb.h
* - uses the arb implementation 
* - Assumes e^{-Et} + e^{-E(T-t)} kernel
* - only 1-d spectral recon for now
* 
************************************************/

#ifndef __hlt_therm_mwh_h__
#define __hlt_therm_mwh_h__

#include <iostream>
#include <vector>
#include <string>
#include <acb.h>
#include <acb_mat.h>
#include <acb_hypgeom.h>
#include <arb_hypgeom.h>


namespace Chimera {
    class HLTRecon {
    private:
        arb_t lambda, epsilon, E0, Ecur;
        int ndigits, tmin, tmax, T;
        arb_mat_t cov, M, v, q, invM;
        std::string ker_type, constraint;
        bool real;
    
    // LU decomposition
    void ludcmp(arb_mat_t& sol, arb_mat_t& mat, arb_mat_t& vec) {
        int ec = arb_mat_solve(sol, mat, vec, ndigits);
        if (ec == 0)
            throw std::string("matrix could not be inverted.");
    }


// Calculates the matrix M
void calcM(const double& c_zero)
{
    int n=tmax-tmin;
    arb_mat_init(M,n,n);

    arb_t ct0;
    arb_init(ct0);
    arb_set_d(ct0,c_zero);

    arb_t one_m_lambda, ta, ea, ma, tmp1, tmp2;
    arb_init(one_m_lambda);
    arb_init(ta);
    arb_init(ea);
    arb_init(ma);
    arb_init(tmp1);    
    arb_init(tmp2);

    arb_neg(ea, E0);
    arb_exp(ea, ea, ndigits); // exp(-E_0)

    arb_one(ma);
    arb_sub(one_m_lambda, ma, lambda, ndigits);

    for (int i = 0; i<n; i++)
        for (int j = 0; j<=i; j++) {
            arb_set_ui(ta, 2*tmin + i + j);
            arb_pow(tmp1, ea, ta, ndigits); //exp(-t+t') E_0 )

            arb_div(ma, tmp1, ta, ndigits); //first term of four
            if (ker_type == "therm") {
                arb_set_ui(ta, T + i - j);
                arb_pow(tmp2, ea, ta, ndigits);
                arb_div(tmp2, tmp2, ta, ndigits);

                arb_add(ma, ma, tmp2, ndigits); //add the second term

                arb_set_ui(ta, T - i + j);
                arb_pow(tmp2, ea, ta, ndigits);
                arb_div(tmp2, tmp2, ta, ndigits);

                arb_add(ma, ma, tmp2, ndigits); //add the third term

                arb_set_ui(ta, 2*T - 2*tmin - i - j);
                arb_pow(tmp2, ea, ta, ndigits);
                arb_div(tmp2, tmp2, ta, ndigits);

                arb_add(ma, ma, tmp2, ndigits); //add the last term
            }
            arb_mul(ma,one_m_lambda,ma,ndigits);
            arb_div(tmp1,arb_mat_entry(cov,i,j),ct0,ndigits);
            arb_mul(tmp2, tmp1, lambda,ndigits);
            arb_add(arb_mat_entry(M,i,j),tmp2,ma,ndigits);  //add lambda of the covariance

            if (i!=j)
                arb_set(arb_mat_entry(M,j,i),arb_mat_entry(M,i,j));
        }

    arb_clear(one_m_lambda);
    arb_clear(ta);
    arb_clear(ea);
    arb_clear(ma);
    arb_clear(ct0);
    arb_clear(tmp1);
    arb_clear(tmp2);
}

// Calculates the vector V
void calcV(const int& t, arb_mat_t& vec) {
    arb_t re1,re2,im,ta,lambda_m_one;
    arb_init(re1);
    arb_init(re2);
    arb_init(im);
    arb_init(ta);
    arb_init(lambda_m_one);

    arb_one(re1);
    arb_sub(lambda_m_one, lambda, re1, ndigits);

    acb_t arg1, arg2, expc, e_int, prod, a_one;
    acb_init(arg1);
    acb_init(arg2);
    acb_init(expc);
    acb_init(e_int);
    acb_init(prod);
    acb_init(a_one);
    acb_one(a_one);

    int n = tmax-tmin;
    for (int i=0; i<n; i++) {
        if (t == 0)
            arb_set_ui(ta, tmin+i);
        else 
            arb_set_ui(ta, T - tmin-i);
            
        arb_mul(re2,Ecur,ta,ndigits);
        arb_neg(re1,re2);//re1=-E t

        arb_mul(im,epsilon,ta,ndigits);//im=eps t
        arb_set(re2,re1);
        arb_addmul(re2,E0,ta,ndigits);//re2=-E t+E0 t

        acb_set_arb_arb(arg1, re1, im);
        acb_set_arb_arb(arg2, re2, im);
        acb_exp(expc, arg1, ndigits); //expc=e^{-E t+i eps t}
        acb_hypgeom_expint(e_int, a_one, arg2, ndigits);
        acb_mul(prod, expc, e_int, ndigits);//prod = expc E1(arg2)
        
        if (real)
            acb_get_imag(re1, prod);
        else 
            acb_get_real(re1, prod);

        arb_mul(arb_mat_entry(vec,i,0),lambda_m_one,re1,ndigits);
    }
    acb_clear(arg1);
    acb_clear(arg2);
    acb_clear(expc);
    acb_clear(e_int);
    acb_clear(prod);
    acb_clear(a_one);

    arb_clear(re1);
    arb_clear(re2);
    arb_clear(im);
    arb_clear(ta);
    arb_clear(lambda_m_one);
}

//used for the thermal kernel
void calcThermV(){
    if(ker_type == "therm") {
        int n = tmax-tmin;
        arb_mat_t vec1, vec2;
        arb_mat_init(vec1, n, 1);
        arb_mat_init(vec2, n, 1);
        
        calcV(0, vec1);
        calcV(1, vec2);

        arb_mat_add(v, vec1, vec2, ndigits);
        
        arb_mat_clear(vec1);
        arb_mat_clear(vec2);
    } else {
        calcV(0, v);
    }

}

   public:
        HLTRecon(const double& lambda_in, const double& epsilon_in, 
        const double& E_0_in, int &ndigits_in, int &tmin_in, int &tmax_in, int &T_in, 
        const std::vector<std::vector<double>> &cmat_in, const std::string& k_in, 
        const std::string& c_in, bool re_in, const double& n_in) : ndigits(ndigits_in), tmin(tmin_in), 
        tmax(tmax_in), T(T_in), ker_type(k_in), real(re_in), constraint(c_in) {
             
            std::cerr << "nbits = "<<ndigits<< ", tmin = " << tmin << ", tmax = " << tmax <<std::endl;
            // Initialize and convert parameters from double to ARB
            arb_init(lambda);
            arb_init(epsilon);
            arb_init(E0);
            arb_init(Ecur);
            
            arb_set_d(lambda,lambda_in);
            arb_set_d(epsilon,epsilon_in);
            arb_set_d(E0,E_0_in);
            
            int n = tmax-tmin;
            double tmp = 0;
            arb_mat_init(cov,n,n);
            

            if (lambda_in != 0) {
                for (int i=0;i<n;i++)
                    for (int j=0;j<=i;j++) {
                        tmp = cmat_in[i][j];
                        arb_set_d(arb_mat_entry(cov,i,j), tmp);
                        if(i != j)
                            arb_set_d(arb_mat_entry(cov,j,i), tmp);
                    }
            } else 
                arb_mat_zero(cov);

            arb_mat_init(M, n, n);
            arb_mat_init(invM, n, n);
            arb_mat_init(v,n,1);
            arb_mat_init(q,n,1);

            calcM(n_in);            
            calcThermV();


            arb_mat_inv(invM, M, ndigits);
            arb_mat_solve(q, M, v, ndigits);
        }

        ~HLTRecon() {
            
            arb_clear(lambda);
            arb_clear(epsilon);
            arb_clear(E0);
            arb_clear(Ecur);

            arb_mat_clear(cov);
            arb_mat_clear(M);
            arb_mat_clear(invM);
            arb_mat_clear(v);
            arb_mat_clear(q);
        }
        
        //find q-vector for new energy
        void setE(const double& energy) {
            arb_set_d(Ecur,energy);
            if(constraint == "no") {
                calcThermV();
            } else if(constraint == "de") {
                //calcDerivativeV();
            }
            arb_mat_solve(q, M, v, ndigits);
        }

        //find q-vector for new epsilon
        void setEps(const double& eps) {
            arb_set_d(epsilon, eps);
            if(constraint == "no") {
                calcThermV();
            } else if(constraint == "de") {
                //calcDerivativeV();
            }
            arb_mat_solve(q, M, v, ndigits);
        }

        //find most effective epsilon
        double getOptEps(){
            int nt = tmax - tmin;

            arb_t norm, sum1, sum2, frac, tmp1, tmp2, eE0, Ediff;

            arb_init(norm);
            arb_init(sum1);
            arb_init(sum2);
            arb_init(frac);
            arb_init(tmp1);
            arb_init(tmp2);
            arb_init(eE0);
            arb_init(Ediff);

            arb_zero(sum1);
            arb_zero(sum2);

            arb_neg(eE0, E0);
            arb_exp(eE0, eE0, ndigits);

            arb_sub(Ediff, Ecur, E0, ndigits);
            for (int t = 0; t < nt; t++) {
                for (int s = 0; s < nt; s++){
                    arb_set_d(tmp1, t + s + 2*tmin);
                    arb_pow(norm, eE0, tmp1, ndigits);
                    arb_mul(frac, arb_mat_entry(q, t, 0), arb_mat_entry(q, s, 0), ndigits);
                    arb_mul(frac, frac, norm, ndigits); 
                    arb_div(frac, frac, tmp1, ndigits);
                    
                    arb_mul(norm, Ediff, tmp1, ndigits);
                    arb_set_d(tmp2, 1);
                    arb_add(norm, norm, tmp2, ndigits);
                    arb_set_d(tmp2, 2);
                    arb_mul(norm, norm, tmp2, ndigits);

                    arb_div(norm, norm, tmp1, ndigits);
                    arb_div(norm, norm, tmp1, ndigits);

                    arb_addmul(norm, Ediff, Ediff, ndigits);

                    arb_mul(norm, frac, norm, ndigits);

                    arb_add(sum1, sum1, frac, ndigits);
                    arb_add(sum2, sum2, norm, ndigits);

                }
            }

            arb_div(norm, sum2, sum1, ndigits);

            double val = arf_get_d(arb_midref(norm),ARF_RND_NEAR);

            arb_clear(norm);
            arb_clear(sum1);
            arb_clear(sum2);
            arb_clear(frac);
            arb_clear(tmp1);
            arb_clear(tmp2);
            arb_clear(eE0);
            arb_clear(Ediff);

            return val;
        }

        // find the smearing function delta
        void getDelta(const double &min_in, const double &max_in, const int &num, std::vector<double>& res){
            int nt = tmax - tmin; 
            res = std::vector<double>(num+1,0);
            
            arb_t expE, tmp, min, max, delta;
            
            arb_init(expE);
            arb_init(tmp);
            arb_init(min);
            arb_init(max);            
            arb_init(delta);
            

            for (int i = 0; i <= num; i++){
                arb_set_d(min, min_in);
                arb_set_d(max, max_in);

                arb_set_d(tmp, num);
                
                arb_sub(expE, max, min, ndigits);
                arb_div(tmp, expE, tmp, ndigits);
                
                arb_set_d(max, i);
                arb_mul(expE, max, tmp, ndigits);
                arb_add(tmp, min, expE, ndigits);        

                arb_exp(expE, tmp, ndigits);
                
                arb_zero(delta);
                for (int j = 0; j < nt; j++) {
                    arb_pow_ui(tmp, expE, j+tmin, ndigits);
                    arb_div(tmp, arb_mat_entry(q, j, 0), tmp, ndigits);
                    arb_add(delta, delta, tmp, ndigits);
                }
                res[i] = arf_get_d(arb_midref(delta),ARF_RND_NEAR); 
            }    
            arb_clear(expE);
            arb_clear(tmp);
            arb_clear(min);
            arb_clear(max);            
            arb_clear(delta);
        }

        //find the reconstruction
        double getRho(const std::vector<double> &c_t) {
            int nt = tmax - tmin;
            
            arb_t rhat,ct, tmp;
            arb_init(rhat);
            arb_init(ct);
            arb_init(tmp);

            arb_zero(rhat);
                    
            for (int i = 0; i < nt; i++) {
                arb_set_d(ct,c_t[i]);
                arb_mul(tmp, ct, arb_mat_entry(q, i, 0), ndigits);
                arb_add(rhat, rhat, tmp, ndigits);
                //std::cerr << arf_get_d(arb_midref(arb_mat_entry(q, i, 0)),ARF_RND_NEAR) << std::endl;
            }
            
            double rho_hat = arf_get_d(arb_midref(rhat),ARF_RND_NEAR);
            arb_clear(rhat);
            arb_clear(ct);
            arb_clear(tmp);

            return rho_hat;
        }

        //find the W-function
        double getW(const double& om1, const double& om2) {
            int nt = tmax - tmin;
            
            arb_t ea1, ea2, tmp1, tmp2, tmp3, sum;

            arb_init(ea1);
            arb_init(ea2);
            arb_init(tmp1);
            arb_init(tmp2);
            arb_init(tmp3);
            arb_init(sum);
            
            arb_set_d(ea1, om1);
            arb_set_d(ea2, om2);

            arb_neg(tmp1, ea1);
            arb_neg(tmp2, ea2);

            arb_exp(ea1, tmp1, ndigits);
            arb_exp(ea2, tmp2, ndigits);

            arb_mat_inv(invM, M, ndigits);

            for (int i = 0; i < nt; i++){
                for (int j = 0; j < nt; j++) {
                    arb_pow_ui(tmp3, ea1, i + tmin, ndigits);
                    arb_pow_ui(tmp2, ea2, j + tmin, ndigits);
                    if (ker_type == "therm") {
                        arb_pow_ui(tmp1, ea1, T - i - tmin, ndigits);
                        arb_add(tmp3, tmp3, tmp1, ndigits);

                        arb_pow_ui(tmp1, ea2, T - j - tmin, ndigits);
                        arb_add(tmp2, tmp1, tmp2, ndigits);
                    }
 
                    arb_mul(tmp3, tmp2, tmp3, ndigits);
                    arb_addmul(sum, tmp3, arb_mat_entry(invM, i, j), ndigits);
                }
            }

            arb_one(tmp1);
            arb_sub(tmp2, tmp1, lambda, ndigits);

            arb_mul(tmp3, sum, tmp2, ndigits);

            double W = arf_get_d(arb_midref(tmp3), ARF_RND_NEAR);
            return W;
        }

    };        
}
#endif
