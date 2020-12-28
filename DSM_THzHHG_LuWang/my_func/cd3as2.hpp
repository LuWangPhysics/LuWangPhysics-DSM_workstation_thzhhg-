#ifndef CD3AS2_H
#define CD3AS2_H
#include <cmath>
#include <algorithm>
#include <functional>
#include <tuple>
#include "../my_func/include.hpp"
#include "../my_plot/include.hpp"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
CONST C;
struct cd3as2{
    double gamma_e;
    double gamma_i;
    double T;
    double e_fermi;
    double g;
    int N_f,N_e;
    int counter=0;
    std::vector<double> v;

    //refractive index with chi1 at different frequencies
    Eigen::ArrayXcd n_thz,n_p1,n_p2;
    Eigen::ArrayXcd a_thz,a_p2,a_p1;


    Eigen::ArrayXcd sigma3_3thz;

    

    Eigen::ArrayXd energy,f,f_d,f_dd;
    Eigen::ArrayXd w_cv;
    double d_energy;
    double omega_max;
    
    cd3as2( const std::vector<double>& v_,const double & e_fermi_,const double &gamma_e_,const double& gamma_i_,const double & T_) {
        v=v_;
        e_fermi=e_fermi_;
        gamma_e=gamma_e_;
        gamma_i=gamma_i_;
        T=T_;
        g=4.0;


     
    
    }
    ~cd3as2(){};
 

    std::complex<double> w_e( const double & w1){return (w1+C.II*gamma_e);};
    std::complex<double> w_i(const double & w1){return (w1+C.II*gamma_i);};
    Eigen::ArrayXcd w_e( const  Eigen::ArrayXd  & w1){return (w1+C.II*gamma_e);};
    Eigen::ArrayXcd w_i(const  Eigen::ArrayXd  & w1){return (w1+C.II*gamma_i);};


    std::tuple<Eigen::ArrayXcd, Eigen::ArrayXcd> n(Eigen::ArrayXd& w1);
    Eigen::ArrayXcd sigma3(Eigen::ArrayXd& w1,Eigen::ArrayXd& w2, Eigen::ArrayXd& w3);
    void integral_check(const Eigen::ArrayXcd& A);
    template <typename T1>
    void assign_chi1(T1 &tf);
    template <typename T1>
    void assign_sigma3(T1 & tf);
    
    Eigen::ArrayXcd sigma_i_f(const Eigen::ArrayXd & w1,const Eigen::ArrayXd &w2,const Eigen::ArrayXd &w3);
    Eigen::ArrayXcd sigma_e_f(const Eigen::ArrayXd& w1,const Eigen::ArrayXd& w2,const Eigen::ArrayXd& w3);
    Eigen::ArrayXcd sigma_ie_f(const Eigen::ArrayXd& w1,const Eigen::ArrayXd& w2,const Eigen::ArrayXd& w3);
    Eigen::ArrayXcd sigma_ei_f(const Eigen::ArrayXd& w1,const Eigen::ArrayXd& w2,const Eigen::ArrayXd& w3);
};

template <typename T1>
void cd3as2::assign_sigma3(T1 & tf){
  Eigen::ArrayXd w1,w2,w3;
//  make sure that w1+w2+w3 the first element is negative frequencies

   w1=tf.omega_thz;
   w2=tf.omega_thz;  
   w3=tf.omega_thz;
   sigma3_3thz=sigma3(w1,w2,w3);


  
};
template <typename T1>
void cd3as2::assign_chi1(T1 &tf){

 

   std::tie(n_thz,a_thz)=n(tf.omega_thz);

    
};


  
    //purly intra band-------------------------------------------------
Eigen::ArrayXcd  cd3as2::sigma_i_f(const Eigen::ArrayXd & w1,const Eigen::ArrayXd &w2,const Eigen::ArrayXd &w3)
     { 
         Eigen::ArrayXcd sigma_f=C.II*pow(C.e,4)*g*pow(v[0],3)/(5.0*pow(C.hbar,3)*pow(C.pi,2)*v[1]*v[2])/(w_i(w1)*w_i(w1+w2)*w_i(w1+w2+w3));
       
         return sigma_f;
     };

    //purly inter band  //---------------------------------------------- 
   
Eigen::ArrayXcd  cd3as2::sigma_e_f(const Eigen::ArrayXd& w1,const Eigen::ArrayXd& w2,const Eigen::ArrayXd& w3){
        
        Eigen::ArrayXcd coe_eee=-C.II*g*pow(C.e,4)*pow(v[0],3)/(15.0*pow(C.hbar,3)*pow(C.pi,2)*v[1]*v[2]*w_i(w1+w2));
        Eigen::ArrayXcd temp,sigma_e;
        sigma_e=Eigen::ArrayXd::Zero(N_f);

        #pragma omp parallel for private(temp)
      
        for(int i=0;i<N_f;i++){
            temp=f/energy/(w_cv-w_e(w1[i]))*(1/(w_cv-w_e(w1[i]+w2[i]+w3[i]))-1/(w_cv+w_e(w1[i]+w2[i]+w3[i])));
            integral_check(temp);
            sigma_e[i]=coe_eee[i]*temp.sum()*d_energy;
        }
        
        return sigma_e;

    };

    
        // inter-intra band ie //---------------------------------------------- 
Eigen::ArrayXcd cd3as2::sigma_ie_f(const Eigen::ArrayXd& w1,const Eigen::ArrayXd& w2,const Eigen::ArrayXd& w3){
      
        Eigen::ArrayXcd coe=-C.II*g*pow(C.e,4)*pow(v[0],3)/(30.0*pow(C.hbar,3)*pow(C.pi,2)*v[1]*v[2])/w_i(w1+w2+w3);
        Eigen::ArrayXcd temp,sigma_ie;
        sigma_ie=Eigen::ArrayXd::Zero(N_f);
        #pragma omp parallel for private(temp)
        for(int i=0;i<N_f;i++){
            temp=4.0*f/energy/(w_cv-w_e(w1[i]))/w_i(w1[i]+w2[i]);
            integral_check(temp);
            sigma_ie[i]=coe[i]*temp.sum()*d_energy;
         }


        Eigen::ArrayXcd diff,temp2;
        
        
        #pragma omp parallel for private (temp,diff,temp2)
        for(int i=0;i<N_f;i++){
            diff=Eigen::ArrayXcd::Zero(N_e);
            temp=f/(w_cv-w_e(w1[i]));
            diff.head(N_e-1)=(temp.tail(N_e-1)-temp.head(N_e-1))/d_energy;
            diff[N_e-1]=diff[N_e-2];
            temp2=coe[i]/(w_cv-w_e(w1[i]+w2[i]))*(diff-2.0*temp/energy-f_d/w_i(w1[i]));
            integral_check(temp2);
            sigma_ie[i]+=temp2.sum()*d_energy;
        }
      

        return sigma_ie;
    };

    

    
    //  more inter intra band ei--------------------------------------------
Eigen::ArrayXcd  cd3as2::sigma_ei_f(const Eigen::ArrayXd& w1,const Eigen::ArrayXd& w2,const Eigen::ArrayXd& w3){
    
        std::complex<double> coe=-C.II*g*pow(C.e,4)*pow(v[0],3)/(30.0*pow(C.hbar,3)*pow(C.pi,2)*v[1]*v[2]);
       Eigen::ArrayXcd temp=Eigen::ArrayXd::Zero(N_e);
       Eigen::ArrayXcd sigma_ei=Eigen::ArrayXd::Zero(N_f);
        #pragma omp parallel for private(temp)
        for(int i=0;i<N_f;i++){
            temp=(4.0*f_d+energy*f_dd)/(w_cv-w_e(w1[i]+w2[i]+w3[i]));
            integral_check(temp);
            sigma_ei[i]=coe*temp.sum()*d_energy/(w_i(w1[i])*w_i(w1[i]+w2[i]));
         }

         
    
        Eigen::ArrayXcd temp2;
        Eigen::ArrayXcd diff;
        
    
       #pragma omp parallel for private (temp,diff,temp2)

        for(int i=0;i<N_f;i++){
            diff=Eigen::ArrayXcd::Zero(N_e);
            temp=f/(w_cv-w_e(w1[i]));  
            diff.head(N_e-1)=(temp.tail(N_e-1)-temp.head(N_e-1))/d_energy;    
            diff[N_e-1]=diff[N_e-2];
       
            temp2=coe*(2.0*w_cv-w_e(w1[i]+w2[i]+w3[i]))*(diff-2.0*temp/energy-f_d/w_i(w1[i]))/(w_cv-w_e(w1[i]+w2[i]))/((w_cv-w_e(w1[i]+w2[i]+w3[i])).pow(2));
            integral_check(temp2);
            sigma_ei[i]+=temp2.sum()*d_energy;
        }
       

        return sigma_ei;
    };
    
    
//************************************************************************
    



#endif
