#ifndef INPUT_FIELD_H
#define INPUT_FIELD_H
#include <Eigen/Dense>
#include <fftw3.h>
#include <math.h> 
//#include "my_func/include.hpp"
CONST C;
struct E_field{
    Eigen::ArrayXcd et_thz;

    
    Eigen::ArrayXcd uw_thz;


    Eigen::ArrayXcd k_thz;
    double N_scale; //normalized the E field to reduce computational effort
    int N_f;
    int N_f_thz;
    int N_half;
    double df;
    double dt;
    
    // electric field in time domain E=E_0exp(-t^2/tau^2)
    // electric field in frequency domian E=E_0*tau*sqrt(pi)*exp(-tau^2*delta_w^2/4)
    // electric field in angular frequency domain=E=E_0*tau*sqrt(1/2)*exp(-tau^2*delta_w^2/4)
    template<typename T1,typename T2>
     E_field(T1 &tf,T2 & chi,double & tau1,double & tau2,double &E_0,int &s){

        
         N_f=tf.N_f;
         N_f_thz=tf.N_f_thz;
         N_half=N_f/2;
         dt=tf.dt;
         df=tf.df;

         k_thz=chi.n_thz*tf.omega_thz/C.c;


         
         //for thz harmonics
 
         double w_0=2*C.pi*tf.f_c_thz;
         uw_thz=2*C.pi*E_0*pow(s*1e12/w_0,s)*(s/w_0)*pow(tf.omega_thz/1e12,s)*exp(-s*tf.omega_thz/w_0)/tgamma(s+1);
         uw_thz.head(N_half)= uw_thz.tail(N_half).reverse();
         

        
        
         
    };
    template<typename B1>
    void time_domain(B1& fftw,double &z);
};


template<typename B1>
void E_field::time_domain(B1& fftw, double &z){
    

    

// fftw.shift(&fftw.E_thz[0],&fftw.E_thz[N_half-1],&fftw.E_thz[N_half]);
// for fftw c to r, the complex array is N+1 one diginit more, and the 
// first element is kept to be zero
#pragma omp parallel for
     for(int i=0;i<N_half;i++)
          {
      
           fftw.E_thzcr[1+i][0]=(uw_thz[i+N_half]*exp(C.II*k_thz[i+N_half]*z)).real();
           fftw.E_thzcr[1+i][1]=(uw_thz[i+N_half]*exp(C.II*k_thz[i+N_half]*z)).imag();
              
           
          };

            fftw.E_thzcr[0][0]=0;
            fftw.E_thzcr[0][1]=0;

   
   fftw_execute_dft_c2r(fftw.c2r,&fftw.E_thzcr[0],&fftw.E_t_thz[0]);
   


}

#endif
