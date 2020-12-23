#include <sys/stat.h>
template<typename T1>
std::string make_string(double& E_01,double & E_02,double &T,double &tau_fwhm,double&E_fermi,double& dz,T1& C,double & gamma_e){

std::stringstream E_text2,E_text1,s_gamma;
std::string s_efield2,s_efield1;
if(E_02>1e5){
E_text2 << std::fixed << std::setprecision(2) <<E_02*1e-5 ;                                  // 2 digit significant number of the E_01
s_efield2 = E_text2.str()+"E5V_per_m";}else{
E_text2 << std::fixed << std::setprecision(2) <<E_02*1e-3 ;                                  // 2 digit significant number of the E_01
s_efield2 = E_text2.str()+"E3V_per_m";
}


E_text1 << std::fixed << std::setprecision(2) <<E_01*1e-6 ;                                  // 2 digit significant number of the E_01
s_efield1 = E_text1.str()+"MV_per_m";


s_gamma<<std::fixed << std::setprecision(3) <<gamma_e*150e-15;
std::string extra_name="E_01_"+s_efield1+"E_02_"+s_efield2+std::to_string((int)T)+"k_"+std::to_string((int)(tau_fwhm*1e15))+"fs";
extra_name+="_"+std::to_string((int)(E_fermi*1e3/C.e))+"em3_"+s_gamma.str()+"gamma_dz"+std::to_string((int)(dz*1e15))+"e15";
extra_name+=C.name_tail;
std::string save_path="/home/luwang/DSM_workstation/my_output/"+extra_name+"/";
if (mkdir(save_path.c_str(),0777)!= 0){mkdir(save_path.c_str(),0777);}


    return save_path;
};



template<typename T1>
std::string make_string(double & E_02,double &T,double &tau_fwhm,double&E_fermi,double& dz, T1& C,double & gamma_e){

std::stringstream E_text,s_gamma;
std::string s_efield2;
if(E_02>1e5){
E_text << std::fixed << std::setprecision(2) <<E_02*1e-5 ;                                  // 2 digit significant number of the E_01
s_efield2 = E_text.str()+"E5V_per_m";}else{
E_text << std::fixed << std::setprecision(2) <<E_02*1e-3 ;                                  // 2 digit significant number of the E_01
s_efield2 = E_text.str()+"E3V_per_m";
}




s_gamma<<std::fixed << std::setprecision(3) <<gamma_e*150e-15;
std::string extra_name="E_01_"+s_efield2+std::to_string((int)T)+"k_"+std::to_string((int)(tau_fwhm*1e15))+"fs";
extra_name+="_"+std::to_string((int)(E_fermi*1e3/C.e))+"em3_"+s_gamma.str()+"gamma_dz"+std::to_string((int)(dz*1e15))+"e15";
extra_name+=C.name_tail;
std::string save_path="/home/luwang/DSM_workstation/my_output/"+extra_name+"/";
if (mkdir(save_path.c_str(),0777)!= 0){mkdir(save_path.c_str(),0777);}


    return save_path;
}
