#pragma once
#include <iostream>
//#include "my_plot/include.hpp"


    
Eigen::ArrayXcd  cd3as2::sigma3( Eigen::ArrayXd& w1, Eigen::ArrayXd& w2, Eigen::ArrayXd& w3){
 
    //check the final frequency array is [f_small f_large]
    if((w1[N_f-1]+w2[N_f-1]+w3[N_f-1])<(w1[0]+w2[0]+w3[0])){
       my_print("array sequence error, f_min>f_max");
       exit(0);
    }

    
    w_cv=energy*2.0/C.hbar;
    
    Eigen::MatrixXd m(3,N_f);
    Eigen::ArrayXcd sigma_i=Eigen::ArrayXcd::Zero(N_f);
    Eigen::ArrayXcd sigma_e=Eigen::ArrayXcd::Zero(N_f);
    Eigen::ArrayXcd sigma_ie=Eigen::ArrayXcd::Zero(N_f);
    Eigen::ArrayXcd sigma_ei=Eigen::ArrayXcd::Zero(N_f);
    m.row(0)=w1;
    m.row(1)=w2;
    m.row(2)=w3;

    int arr[] = { 0, 1, 2 };
    std::sort(arr, arr + 3);
     
   do {
       //permute the order of arr

      sigma_i+=sigma_i_f(m.row(arr[0]),m.row(arr[1]),m.row(arr[2]));
      sigma_e+=sigma_e_f(m.row(arr[0]),m.row(arr[1]),m.row(arr[2]));
      sigma_ie+=sigma_ie_f(m.row(arr[0]),m.row(arr[1]),m.row(arr[2]));
      sigma_ei+=sigma_ei_f(m.row(arr[0]),m.row(arr[1]),m.row(arr[2]));

   } while (std::next_permutation(arr, arr + 3));

    Eigen::ArrayXcd sigma_3=(sigma_i+sigma_e+sigma_ei+sigma_ie)/6.0;
    
    Eigen::ArrayXd w_final=(w1+w2+w3);
 // if array extends to the negativ frequency, get the sigma(-w)=sigma(w)*  
    if ((w_final<0.0).any()){
        int p=floor(N_f/2);
        double w_c=w1[p]+w2[p]+w3[p]; 
        
        int N_c=floor(w_c/(w1[1]+w2[1]+w3[1]-w1[0]-w2[0]-w3[0]));

        for (int i=0;i<(p-N_c-1);i++){
        sigma_3[i]=conj(sigma_3[N_f-2*N_c-i-1]);
        }
    }

    

    return sigma_3;
 
};


void cd3as2::integral_check(const Eigen::ArrayXcd& A){
    double int_max=A.abs().maxCoeff();
    if (((std::abs(A[0])/int_max)>1.0e-2)||((std::abs(A[N_e-1])/int_max)>1.0e-2)){
    my_print("change integral domain \n");
   
    }
    
};
