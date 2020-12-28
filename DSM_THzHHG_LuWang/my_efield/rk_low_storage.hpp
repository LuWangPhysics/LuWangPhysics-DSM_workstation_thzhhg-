#ifndef RK_LOW_SOTRAGE_H
#define RK_LOW_SOTRAGE_H

#include "../my_func/include.hpp"
#include<complex.h>
#include<fftw3.h>
#include <cmath>


struct low_storage_RK {
    
     Eigen::Array3d ARK;
     Eigen::Array3d BRK;
  
     double dz;
     Eigen::ArrayXcd P_ir1;
     Eigen::ArrayXcd P_ir2;
     Eigen::ArrayXcd P_thz;

     Eigen::ArrayXcd rk_ir1;
     Eigen::ArrayXcd rk_ir2;
     Eigen::ArrayXcd rk_thz;
   
     double df;
     double dt;
     int N_f;
     int N_half;
     int N_f_thz;
     std::complex<double> tmp;
     double tmp_r;
//constructor without input
     low_storage_RK(){}
     template<typename B1,typename B2>
     low_storage_RK(B2& tf,B1 & z_struc):dz(z_struc.dz){
          dt=tf.dt;
          df=tf.df;
          N_f=tf.N_f;
          N_f_thz=tf.N_f_thz;
          N_half=N_f/2;

         ARK<<0,(double)-5/9,(double)-153/128;
               
         BRK<<(double)1/3,(double)15/16,(double)8/15;
      
          rk_ir1=Eigen::ArrayXcd::Zero(tf.N_f);
          rk_ir2=Eigen::ArrayXcd::Zero(tf.N_f);
          rk_thz=Eigen::ArrayXcd::Zero(tf.N_f_thz);
          P_ir1=Eigen::ArrayXcd::Zero(tf.N_f);
          P_ir2=Eigen::ArrayXcd::Zero(tf.N_f);
          P_thz=Eigen::ArrayXcd::Zero(tf.N_f_thz);
        
    }
    ~low_storage_RK(){}
    template<typename B1,typename B2>
    void  set_low_storage_RK(B2& tf,B1 & z_in){
          dz=z_in;
          dt=tf.dt;
          df=tf.df;
          N_f=tf.N_f;
          N_f_thz=tf.N_f_thz;
          N_half=N_f/2;

         ARK<<0,(double)-5/9,(double)-153/128;
               
         BRK<<(double)1/3,(double)15/16,(double)8/15;
      
          rk_ir1=Eigen::ArrayXcd::Zero(tf.N_f);
          rk_ir2=Eigen::ArrayXcd::Zero(tf.N_f);
          rk_thz=Eigen::ArrayXcd::Zero(tf.N_f_thz);
          P_ir1=Eigen::ArrayXcd::Zero(tf.N_f);
          P_ir2=Eigen::ArrayXcd::Zero(tf.N_f);
          P_thz=Eigen::ArrayXcd::Zero(tf.N_f_thz);
        
    };
    template<typename B1,typename B2,typename B3,typename B4,typename B5>
    void   method(B1& E,B2& fftw,B3& chi, double z,B4 &tf,B5& mypython);
    template<typename B2,typename B1,typename B3,typename B4,typename B5>
    void polarization_eei(B2& fftw, B1& E,B3& chi,B4 & tf,double z,B5& mypython);
    template<typename B2,typename B1,typename B3,typename B4,typename B5>
    void polarization_thz3hg(B2& fftw, B1& E,B3& chi,B4 & tf,double z,B5& mypython);
};

template<typename B1,typename B2,typename B3,typename B4,typename B5>
void low_storage_RK::method(B1& E,B2& fftw,B3& chi, double z,B4& tf,B5& mypython)
 {      
     
        for(int  m=0;m<3;m++)
       {     
            //-----------------------------------
            //calculate E field in time domain;
            //-----------------------------------
           //normalized the amplitude to the N_scale times smaller

            E.time_domain(fftw,z);

            //-----------------------------------
            //calculate polarisation
            //-----------------------------------
            polarization_eei(fftw,E,chi,tf,z,mypython);
           
           // polarization_thz3hg(fftw,E,chi,tf,z,mypython);
     
            //------------------------------------
            //update using low-storage rk method
            //------------------------------------
  
            rk_ir1*=ARK[m];
            rk_ir1+=dz*P_ir1;
            rk_ir2*=ARK[m];
            rk_ir2+=dz*P_ir2;
            rk_thz*=ARK[m];
            rk_thz+=dz*P_thz;
    
            E.uw_ir1+=BRK[m]*rk_ir1;
            E.uw_ir2+=BRK[m]*rk_ir2;
            E.uw_thz+=BRK[m]*rk_thz;

            
        };

       
};





#endif
