#pragma once
#ifndef J_TEST_H
#define J_TEST_H
#include "../my_func/include.hpp"
template<typename B1,typename B2,typename B3,typename B4,typename B5,typename B6,typename B7>
void J_func(B1& tf, B2&C, B3&E, B4&mypython,B5&rk_low,B6&chi,B7& fftw,double & z,double& dz){
    
                
           
            rk_low.set_low_storage_RK(tf,dz);
          
            rk_low.polarization_eei(fftw,E,chi,tf,z,mypython);
          
               
            std::string name="J_thz"+std::to_string((int)(z*1e9))+"nm";
            Eigen::ArrayXcd phase=exp(C.II*E.k_thz*z);
            Eigen::ArrayXcd J_temp=(-C.II*C.eps*tf.omega_thz)*(pow(chi.n_thz,2)-1)*E.uw_thz*phase-rk_low.P_thz*2*chi.n_thz*C.c*C.eps*phase;
            save_array(J_temp,mypython.save_path+name+".txt");
            name="J_thz1"+std::to_string((int)(z*1e9))+"nm";
             J_temp=(-C.II*C.eps*tf.omega_thz)*(pow(chi.n_thz,2)-1)*E.uw_thz*phase;
            save_array(J_temp,mypython.save_path+name+".txt");
            
            name="J_ir1"+std::to_string((int)(z*1e9))+"nm";
            phase=exp(C.II*E.k_ir1*z);
            J_temp=(-C.II*C.eps*tf.omega_p1)*(pow(chi.n_p1,2)-1)*E.uw_ir1*phase-rk_low.P_ir1*2*chi.n_p1*C.c*C.eps*phase;
            save_array(J_temp,mypython.save_path+name+".txt");
            
            name="J_ir11"+std::to_string((int)(z*1e9))+"nm";
            J_temp=(-C.II*C.eps*tf.omega_p1)*(pow(chi.n_p1,2)-1)*E.uw_ir1*phase;
            save_array(J_temp,mypython.save_path+name+".txt");
            
            name="J_ir2"+std::to_string((int)(z*1e9))+"nm";
            phase=exp(C.II*E.k_ir2*z);
            J_temp=(-C.II*C.eps*tf.omega_p2)*(pow(chi.n_p2,2)-1)*E.uw_ir2*phase-rk_low.P_ir2*2*chi.n_p2*C.c*C.eps*phase;
            save_array(J_temp,mypython.save_path+name+".txt");
            name="J_ir21"+std::to_string((int)(z*1e9))+"nm";
            J_temp=(-C.II*C.eps*tf.omega_p2)*(pow(chi.n_p2,2)-1)*E.uw_ir2*phase;
            save_array(J_temp,mypython.save_path+name+".txt");
    
}

#endif
