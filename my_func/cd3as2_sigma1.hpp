#pragma once
#include "../my_plot/include.hpp"
#include <tuple>
std::tuple<Eigen::ArrayXcd, Eigen::ArrayXcd> cd3as2::n(Eigen::ArrayXd & w1){
        Eigen::ArrayXcd par_1,sigma_intra,sigma_inter,chi,chi_x;
    
        par_1=C.II*(pow(C.e,2)*g/pow(2*C.pi*C.hbar,3))/w_i(w1);
        double par_2=4*C.pi*v[0]/(3*v[2]*v[1]);
        sigma_intra=par_1*par_2*(pow(C.kb*T*C.pi,2)/3+pow(e_fermi,2));
    
        //calculte interband transition numerically
      
        N_f=w1.size();
      
       
      
        sigma_inter=Eigen::ArrayXd::Zero(N_f);
        //the integration range has to be larger than the maximum of the third harmonics of optical
        N_e=8e4;
     
        double energy_h=3.0*omega_max*C.hbar/2.0+14.0*e_fermi;
        double energy_l=-(3*omega_max*C.hbar/2.0+9.0*e_fermi);
        energy=Eigen::ArrayXd::LinSpaced(N_e,energy_l,energy_h).transpose();

        
        d_energy=(energy[1]-energy[0]);
        double coef2=g*(v[0]/(v[2]*v[1]))*pow(C.e,2)/(24.0*pow(C.pi,2)*C.hbar);
        
        //define fermi districution and the derivatives
        f=1/(exp((energy-e_fermi)/(C.kb*T))+1)-1/(exp((-energy-e_fermi)/(C.kb*T))+1);

        double beta=C.kb*T;

        f_d=-exp((-energy.array()-e_fermi)/beta)/(beta*pow(1+exp((-energy.array()-e_fermi)/beta),2))\
        -exp((energy.array()-e_fermi)/beta)/(beta*pow(1+exp((energy.array()-e_fermi)/beta),2));

        f_dd=-2*exp(2*(-energy.array()-e_fermi)/beta)/(pow(beta,2)*pow(1+exp((-energy.array()-e_fermi)/beta),3))\
        +exp((-energy.array()-e_fermi)/beta)/(pow(beta,2)*pow(1+exp((-energy.array()-e_fermi)/beta),2))\
        +2*exp(2*(energy.array()-e_fermi)/beta)/(pow(beta,2)*pow(1+exp((energy.array()-e_fermi)/beta),3))\
        -exp((energy.array()-e_fermi)/beta)/(pow(beta,2)*pow(1+exp((energy.array()-e_fermi)/beta),2));

        //change the nan to zero due to too low T
        f_d = (isnan(f_d)).select(0, f_d);
        f_dd = (isnan(f_dd)).select(0, f_dd);      


        Eigen::ArrayXcd temp;
        #pragma omp parallel for private(temp)
        for(int i=0;i<N_f;i++){
        temp=f/(energy+C.hbar*w_e(w1[i])/2.0);
        sigma_inter[i]=C.II*coef2*w_e(w1[i])*temp.sum()*d_energy;

        }

        Eigen::ArrayXcd n,alpha;
        chi_x=(sigma_intra+sigma_inter)/(-C.II*C.eps*w1);
        n=sqrt(1+chi_x);
        alpha=0*w1*sqrt(1+chi_x).imag()/C.c;
        return {n,alpha};
};
