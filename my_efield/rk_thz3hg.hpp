#ifndef RK_THZ3HG_H
#define RK_THZ3HG_H

 template<typename B2,typename B1,typename B3,typename B4,typename B5>
void low_storage_RK::polarization_thz3hg(B2& fftw, B1& E,B3& chi,B4 & tf,double z,B5& mypython){
    
            
       //------------------------------------
       //THz generation 3thz
       //------------------------------------     
             int N_half=tf.N_f/2;
              for(int i=0;i<N_f;i++)
               {
       
                tmp_r=std::pow((fftw.E_t_thz[i]),3)*std::pow(df,3);
                fftw.in_r[i]=tmp_r;
        
                }
             
               fftw_execute_dft_r2c(fftw.r2c,&fftw.in_r[0],&fftw.E_thzcr[0]);

              //generate the third harmonics
              //since e(iwt) is included in the thz spectrum, the 3w will automatically be shifted into the frequency array
               for(int i=0;i<N_half;i++)
               {// r2c the first array digit shouldn't be used
                   tmp.real(fftw.E_thzcr[i+1][0]);
                   tmp.imag(fftw.E_thzcr[i+1][1]);
                   P_thz[i+N_half]=-dt*chi.sigma3_3thz[i+N_half]*tmp*std::conj(std::exp(C.II*E.k_thz[i+N_half]*z))/(C.eps*2*C.c*chi.n_thz[i+N_half]);
                   
           
                }
                



    
};
#endif
