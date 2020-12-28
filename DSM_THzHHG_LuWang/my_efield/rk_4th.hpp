#pragma once
#ifndef RK_4TH_H
#define RK_4TH_H

#include "../my_func/include.hpp"
#include<complex.h>
#include<fftw3.h>
#include <cmath>

struct rk_4th_order {
  
     double dz; 
     Eigen::ArrayXcd uw_thz_or;

     Eigen::ArrayXcd rk_thz1,rk_thz2,rk_thz3,rk_thz4;
   
     double df;
     double dt;
     int N_f;
     int N_half;
     int N_f_thz;
     std::complex<double> tmp;
     double tmp_r;

     low_storage_RK  rk_low;
     
  
     template<typename B1,typename B2>
     rk_4th_order(B2& tf,B1 & z_struc):dz(z_struc.dz){
         
          rk_low.set_low_storage_RK(tf,dz);
    
          dt=tf.dt;
          df=tf.df;
          N_f=tf.N_f;
          N_f_thz=tf.N_f_thz;
          N_half=N_f/2;
               
       
          rk_thz1=Eigen::ArrayXcd::Zero(tf.N_f_thz);
          rk_thz2=Eigen::ArrayXcd::Zero(tf.N_f_thz);
          rk_thz3=Eigen::ArrayXcd::Zero(tf.N_f_thz);
          rk_thz4=Eigen::ArrayXcd::Zero(tf.N_f_thz);

        
    }
    ~rk_4th_order(){}
    
    template<typename B1,typename B2,typename B3,typename B4,typename B5>
    void   method(B1& E,B2& fftw,B3& chi, double z,B4 &tf,B5& mypython);


};

template<typename B1,typename B2,typename B3,typename B4,typename B5>
void rk_4th_order::method(B1& E,B2& fftw,B3& chi, double z,B4& tf,B5& mypython)
 {     
     


          //-----------------------------------
           //calculate E field in time domain; 
           //rk1
          //-----------------------------------
            //normalized the amplitude to the N_scale times smaller
      

            uw_thz_or=E.uw_thz;
            E.time_domain(fftw,z);
            rk_low.polarization_thz3hg(fftw,E,chi,tf,z,mypython);
            rk_thz1=dz*rk_low.P_thz;
            
        
           // rk2

            E.uw_thz=uw_thz_or+0.5*rk_thz1;
            E.time_domain(fftw,z);
            rk_low.polarization_thz3hg(fftw,E,chi,tf,z+0.5*dz,mypython);

            rk_thz2=dz*rk_low.P_thz;
           // rk3

            E.uw_thz=uw_thz_or+0.5*rk_thz2;
            E.time_domain(fftw,z);
            rk_low.polarization_thz3hg(fftw,E,chi,tf,z+0.5*dz,mypython);

            rk_thz3=dz*rk_low.P_thz;
           // rk4

            E.uw_thz=uw_thz_or+rk_thz3;
            E.time_domain(fftw,z);
            rk_low.polarization_thz3hg(fftw,E,chi,tf,dz+z,mypython);

            rk_thz4=dz*rk_low.P_thz;
//             -------------------------------------
//             update 
//             -------------------------------------------------
// 

              
            
            E.uw_thz=uw_thz_or+(rk_thz1+2*rk_thz2+2*rk_thz3+rk_thz4)/6.0;



};




#endif
