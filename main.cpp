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
#include "DSM_THzHHG_LuWang/my_func/include.hpp"
#include "DSM_THzHHG_LuWang/my_mesh/include.hpp"
#include "DSM_THzHHG_LuWang/my_plot/include.hpp"
#include "DSM_THzHHG_LuWang/my_efield/include.hpp"

int main(int argc, char **argv){ 

//---------------------------------------------------------
//define input pulse parameters
//---------------------------------------------------------
//omp_set_num_threads(omp_get_max_threads());
omp_set_num_threads(56);
    
double T_arr[]= {77};                                                                 //set the temperature scan values
 for (double T_n : T_arr){  
 double E_arr[]= {5e5};                                                               //set E field strength scan values 
for(double E_n:E_arr){

CONST C;
std::string name_add="";                                                              //add any string you want in the saving name
C.name_tail=name_add;
mytime time_count;
time_count.start();


double f_max      = 30e12;                                                             //2*f_max =frequency range 
double df         = 0.05e12;                                                           //frequency step size
double f_thz_3h=1e12;                                                                 //fundamental thz center frequency
double E_01=E_n;                                                                      //V/m
double s=33;                                                                          // the  s thing you use in your paper
//---------------------------------------------------------
//define DSM parameters
//---------------------------------------------------------

double L=10.05e-9;                                                                     // DSM length
double L_save=10e-9;
double dz=0.5e-12;
int N_z=10;                                                                            //saving eff array size
MESH::spacial_z z(L,L_save,dz,N_z);

std::vector<double> v{1.28*1e6,1.3*1e6,0.33*1e6};                                       // [m/s] fermi velocity vx,vy,vz
double T=T_n;                                                                           // temperature           
double E_fermi=0.45*C.e;                                                                // [J] fermi energy =chemical potential   
double gamma_e=1/(150e-15);                                                             // inter band decay rate
double gamma_i=gamma_e;                                                                 // intra band decay rate
//---------------------------------------------------------
//define time space mesh
//---------------------------------------------------------
int NN_f=2*(f_max/df);
int N_f=NN_f-(NN_f%4);                                                                   // N_f , N_f_thz even
int N_f_thz=N_f;
MESH::time_frequency tf(N_f, N_f_thz,df,f_thz_3h);

//---------------------------------------------------------
// material class with chi1 chi3
//---------------------------------------------------------
cd3as2 chi(v,E_fermi,gamma_e,gamma_i,T);
chi.assign_chi1(tf);
chi.assign_sigma3(tf);

//---------------------------------------------------------
// define input electric fields
//---------------------------------------------------------

EW::E_field E(tf,chi,E_01,s);

// ---------------------------------------------------------
// define Fourier transform
// ---------------------------------------------------------
myfft fftw(N_f,N_f_thz, 1,1);
fftw.assign();

// ---------------------------------------------------------
// initialize saving path 
// ---------------------------------------------------------

std::string save_path=make_string(E_01,T,tf.tau_fwhm,E_fermi,dz,C,gamma_e);
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



for (int k=0;k<z.N_total;k++){
   
       z.z_p=dz*k;
       RK.method(E,fftw,chi,z.z_p,tf,mypython);
       my_data.save_loop(k,z,fftw,mypython,E,tf,chi);
       
}


my_data.calculate_eff_spec(E,tf,z,chi,(z.N_z-1));
my_data.save_all(std::to_string((int) z.z_p*1e9)+"nm",fftw,mypython);


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

} //T
} //e0
return 0;
}
