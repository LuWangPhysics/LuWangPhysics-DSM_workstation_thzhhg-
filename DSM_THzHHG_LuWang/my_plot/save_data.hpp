#ifndef SAVE_DATA_H
#define SAVE_DATA_H
#include"save_print.hpp"
//#include "my_func/include.hpp"


struct save_data{
    
      std::string save_path;
      Eigen::ArrayXd thz_spectrum,Energy_cons;
      Eigen::ArrayXd eff_3hg;
      bool cond2,cond3;
      Eigen::ArrayXcd e_thz;
      double Energy_thz,df;
      int N_f;
    
      template<typename T1,typename T2>
      save_data(std::string a,T1& tf,T2 &z){
          save_path=a;
          N_f=tf.N_f;
          df=tf.df;

          Energy_cons=Eigen::ArrayXd::Zero(z.N_z);
      

      //for thz 3hg
          eff_3hg=Eigen::ArrayXd::Zero(z.N_z);
          save_array(tf.omega_thz,save_path+"thz_omega.txt");
          save_array(z.z_0,save_path+"z.txt");
          save_array(tf.t,save_path+"t.txt");
          
    };
      ~save_data(){}
    template<typename T1,typename T2,typename T3,typename T4>
    void calculate_eff_spec(T1 &E,T2& tf,T3& z,T4& chi,int j);
    template<typename T1>
    void save_eff_spec(const std::string& b,T1&mypython);
    template<typename T1,typename T2>
    void save_efield(int N_f,const std::string& b,T1& fftw,T2& mypython);
    template<typename T1,typename T2>
    void save_all(const std::string& b,T1& fftw,T2& mypython);
    template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
    void save_loop(const int & k, T1& z,T2& fftw, T3& mypython,T4&E,T5& tf,T6& chi);
      
};

template<typename T1,typename T2,typename T3,typename T4, typename T5,typename T6>
void save_data::save_loop(const int & k, T1& z,T2& fftw, T3& mypython,T4& E,T5& tf,T6& chi){

       if((k%z.N_z_save)==0){
                    
                   // my_print(z.N_z_iter);
                    calculate_eff_spec(E,tf,z,chi,z.N_z_iter);
                  
                    if(z.N_z_iter>2){
                    
                            cond2=(eff_3hg[z.N_z_iter-2]<eff_3hg[z.N_z_iter-1]);
                            cond3=(eff_3hg[z.N_z_iter-1]>eff_3hg[z.N_z_iter]);
                            if(cond2&&cond3){
                            //my_data.save_eff_spec(std::to_string(z.z_0[k]*1e9)+"nm",mypython);
                            save_all(std::to_string((int)(z.z_0[z.N_z_iter]*1e9))+"nm",fftw,mypython);
                            }
   
                    }
                 
                    my_print(z.N_z_iter);
                    z.N_z_iter++;
        }
       //save at the desired saving location

       if((k%(int)(z.L_save/z.dz))==0){
           save_all(std::to_string((int)(z.z_0[z.N_z_iter]*1e9))+"nm",fftw,mypython);
       }
}
template<typename T1,typename T2,typename T3,typename T4>
void save_data::calculate_eff_spec(T1 &E,T2& tf,T3& z,T4& chi,int j){
CONST C;
   

    e_thz=E.uw_thz*exp(C.II*E.k_thz*z.z_0[j]);

    // use convention E(t)=E1(t)exp(-iwt)+c.c no 1/2 envolved
    // E(w)=E1(w)+E1(-w)

    // optical only contains half spectrum but the THz range contains
    //both the possitive and negative frequency parts.
    thz_spectrum=(chi.n_thz.real()*C.eps*C.c)*e_thz.cwiseAbs2();
    Energy_cons[j]= 2*thz_spectrum.tail(ceil(N_f/2)).sum()*tf.df;
  
    //check if there is numerical error
    // condition when energy conservation breaks ||(Energy_cons[j]>Energy_cons[0])
    if(isnan(eff_3hg[j])||isinf(eff_3hg[j])||(Energy_cons[j]>Energy_cons[0])){
        std::cout<<"numerical error found, sadly abort the program  :("<<std::endl;
        std::string b=std::to_string((int)(z.z_0[z.N_z_iter]*1e9))+"nm";

        save_array(e_thz.abs(),save_path+b+"thz_ef_abs.txt");
        save_array(eff_3hg,save_path+b+"eff.txt");
        save_array(Energy_cons,save_path+"energy_cons.txt");
       
        exit(0);
    }

    //sum up for THZ 3HG efficiency bandwidth from 2f_c_thz to 4 f_c_thz
    int N_start=N_f/2.0+3.0*tf.N_thz_3h-1.0*tf.N_thz_3h;
   
    Energy_thz=0;
    for (int i=0;i<2*tf.N_thz_3h;i++){
    Energy_thz+=2*thz_spectrum[N_start+i]*tf.df;
    }
           
    eff_3hg[j]=100*Energy_thz/Energy_cons[0];
       
        
    
};
template<typename T1>
void save_data::save_eff_spec(const std::string& b,T1&mypython){
    my_print("save eff and spectra to "+save_path+b);
   // save_array(Energy_cons,save_path+"energy_cons.txt");
    save_array(e_thz,save_path+b+"thz_e_w.txt");
    save_array(eff_3hg,save_path+b+"eff_3hg.txt");
    
    
   // mypython.plot_eff_spec(save_path+b);
    
    
};
template<typename T1,typename T2>
void save_data::save_efield(int N_f,const std::string& b,T1& fftw,T2& mypython){
   
    my_print("save efields to "+save_path+b);
    
    //save the entire electric field
    Eigen::ArrayXd E_t_final=Eigen::ArrayXd::Zero(N_f);
 
    for(int i=0;i<N_f;i++){
        E_t_final[i]=(fftw.E_t_thz[i])*df;
    }
    save_array(E_t_final,save_path+b+"et_thz_total.txt");
    
    
    
    //mypython.plot_efield(save_path+b);
}
template<typename T1,typename T2>
void save_data::save_all(const std::string& b,T1& fftw,T2& mypython){
    my_print("save all the files to "+save_path);
    save_eff_spec(b,mypython);
    save_efield(N_f,b,fftw,mypython);
}
#endif
