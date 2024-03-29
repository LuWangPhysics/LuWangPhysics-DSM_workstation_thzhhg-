#include <stdio.h>
#include <iostream>
#include <stdlib.h>                                                                 // include random number
#include <sstream>                                                                  // for precision of the number to string
#include <fstream>                                                                  // read and write file
#include <typeinfo>                                                                 // templete of cwise operation of vector
#include <string>	
#include <time.h>
#include <vector>
#include <thread>
#include <mutex>
#include <iomanip>
#include <limits>
#include <omp.h>
#include "my_func/include.hpp"
#include "my_mesh/include.hpp"
#include "my_plot/include.hpp"
#include "my_efield/include.hpp"

int main(int argc, char **argv){ 

//---------------------------------------------------------
//define input pulse parameters
//---------------------------------------------------------
omp_set_num_threads(omp_get_max_threads());
//omp_set_num_threads(28);
double T_arr[]= {77};
 for (double T_n : T_arr){
 double E_arr[]= {15e6};
for(double E_n:E_arr){


CONST C;
std::string name_add="_save_P";
C.name_tail=name_add;
mytime time_count;
time_count.start();


double lambda_1   = 1000e-9;
double lambda_2   = 2000e-9;
double tau_fwhm1  = 150e-15;                                                          //Input pulse fwhm
double tau1       = tau_fwhm1/(sqrt(2*log(2)));                                       //Converting from FWHM to 1/e^2 value. e^(-t^2/(2*tau^2))
double tau2       =tau1;
double f_max      = 30e12;                                                             //2*f_max =frequency range ,120
double df         = 0.05e12;           //0.05                                          //frequency step size
double f_c_thz=0e12;
double f_thz_3h=1e12;
double omega_1 =  2*C.pi*C.c/lambda_1;                                                 // 2w2-w1=thz
double omega_2 =  2*C.pi*C.c/lambda_2;
    
//---------------------------------------------------------
//define DSM parameters
//---------------------------------------------------------
//  double n_arr[]= {1};
// for(double n_iter:n_arr){
double L=0.55e-12;                                                                      // DSM length
double L_save=0.5e-12;
double dz=0.5e-12;
int N_z=1000;                                                                           //saving eff array size
MESH::spacial_z z(L,L_save,dz,N_z);

std::vector<double> v={1.28*1e6,1.3*1e6,0.33*1e6};                                      // [m/s] fermi velocity vx,vy,vz
double T=T_n;                                                                          // temperature           
//for (int n_iter=2;n_iter<11;n_iter++){

double E_fermi=0.45*C.e;                                                               // [J] fermi energy =chemical potential   
double gamma_e=1/(150e-15);                                                             // inter band decay rate
double gamma_i=gamma_e;                                                                 // intra band decay rate
//---------------------------------------------------------
//define time space mesh
//---------------------------------------------------------
int NN_f=2*(f_max/df);
int N_f=NN_f-(NN_f%4);                                                                       // N_f , N_f_thz even
int N_f_thz=N_f;
MESH::time_frequency tf(N_f, N_f_thz,df,omega_1,omega_2,f_c_thz);


tf.N_thz_3h=f_thz_3h/df;
tf.f_thz_3h=f_thz_3h;
//---------------------------------------------------------
// material class with chi1 chi3
//---------------------------------------------------------
cd3as2 chi(v,E_fermi,gamma_e,gamma_i,T);
chi.assign_chi1(tf);
chi.assign_sigma3(tf);
  

//---------------------------------------------------------
// define input electric fields
//---------------------------------------------------------

double E_01=E_n;        //5e6                                                             //V/m
double E_02=E_n;
EW::E_field E(tf,chi,tau1,tau2,E_01,E_02);

// ---------------------------------------------------------
// define Fourier transform
// ---------------------------------------------------------
myfft fftw(N_f,N_f_thz, 1,1);
fftw.assign();

// ---------------------------------------------------------
// initialize saving path 
// ---------------------------------------------------------
std::string save_path=make_string(E_01,E_02,T,tau_fwhm1,E_fermi,dz,C,gamma_e);
// --------------------------------------------------------------
// numerical method RK method initialization//lowest storage rk method
// ---------------------------------------------------------

//EW::low_storage_RK RK(tf,z);
EW::rk_4th_order RK(tf,z);

// --------------------------------------------------------------
// save output data and plot
// ---------------------------------------------------------
save_data my_data(save_path,tf,z);
python_plot mypython(save_path);

//  std::string name="sigma_iee2";
//  save_array(chi.sigma3_iee_b2,save_path+name+".txt");
//   name="sigma_spm1";
//  save_array(chi.sigma3_spm1,save_path+name+".txt");
// mypython.plot("thz_omega.txt",name+".txt",name);

EW::low_storage_RK  rk_low;   
for (int k=0;k<1;k++){
       //set new z and new inter_z each 
       z.z_p=dz*k;
        E.time_domain(fftw,z.z_p);
       J_func(tf,C,E,mypython,rk_low,chi,fftw,z.z_p,z.dz);
       
}


//my_data.calculate_eff_spec(E,tf,z,chi,(z.N_z-1));
//my_data.save_all(std::to_string((int) z.z_p*1e9)+"nm",fftw,mypython);


// -----------------------------------------------------------
// REMOVE all the input pointers from fftw and mypython
//     VERY IMPORTANT
// -----------------------------------------------------------
    
fftw.my_fftw_clean();
mypython.clean();
time_count.end();
// -----------------------------------------------------------
//for parameters scan
// -----------------------------------------------------------
//} //fermi
} //T
} //e0
return 0;
}
