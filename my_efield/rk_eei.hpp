#ifndef RK_EEI_H
#define RK_EEI_H
 template<typename B2,typename B1,typename B3,typename B4,typename B5>
void low_storage_RK::polarization_eei(B2& fftw, B1& E,B3& chi,B4 & tf,double z,B5& mypython){
    
            

         #pragma omp parallel for private(tmp)
               for(int i=0;i<N_f;i++)
               {
                // ------------------------------------
                //THz back convert to optical 
                //w1=2w2-thz
                //------------------------------------
                tmp=std::pow(fftw.complex(fftw.E_t_ir2[i]),2)*std::conj(fftw.E_t_thz[i])*std::pow(df,3);
                fftw.in[i][0]=tmp.real();
                fftw.in[i][1]=tmp.imag();
                // ------------------------------------
                //spm w1
                //------------------------------------
                tmp=std::pow(fftw.complex(fftw.E_t_ir1[i]),2)*std::conj(fftw.complex(fftw.E_t_ir1[i]))*std::pow(df,3);
                fftw.in_spm[i][0]=tmp.real();
                fftw.in_spm[i][1]=tmp.imag();
                
                }
                

                fftw_execute_dft(fftw.t2f,&fftw.in[0],&fftw.out[0]);
                fftw.shift(&fftw.out[0],&fftw.out[N_half-1],&fftw.out[N_half]);
                
                fftw_execute_dft(fftw.t2f,&fftw.in_spm[0],&fftw.out_spm[0]);
                fftw.shift(&fftw.out_spm[0],&fftw.out_spm[N_half-1],&fftw.out_spm[N_half]);

                 #pragma omp parallel for private(tmp)
                for(int i=0;i<N_f;i++)
                {
                    tmp.real(fftw.out[i][0]);
                    tmp.imag(fftw.out[i][1]);
                    P_ir1[i]=-3.0*dt*chi.sigma3_iee_b1[i]*tmp*std::exp(-C.II*E.k_ir1[i]*z)/(C.eps*2*C.c*chi.n_p1[i]);
                    
                    tmp.real(fftw.out_spm[i][0]);
                    tmp.imag(fftw.out_spm[i][1]);
                    P_ir1[i]+=(-3.0*dt*chi.sigma3_spm1[i]*tmp*std::exp(-C.II*E.k_ir1[i]*z)/(C.eps*2*C.c*chi.n_p1[i]));
                    
                }
            

       //------------------------------------
       //THz back convert to optical w2=w1+thz-w2
       //------------------------------------
               #pragma omp parallel for private(tmp)
               for(int i=0;i<N_f;i++)
               {
           
                tmp=fftw.complex(fftw.E_t_ir1[i])*fftw.E_t_thz[i]*std::conj(fftw.complex(fftw.E_t_ir2[i]))*std::pow(df,3);
                fftw.in[i][0]=tmp.real();
                fftw.in[i][1]=tmp.imag();
                // ------------------------------------
                //spm w2
                //------------------------------------
                tmp=std::pow(fftw.complex(fftw.E_t_ir2[i]),2)*std::conj(fftw.complex(fftw.E_t_ir2[i]))*std::pow(df,3);
                fftw.in_spm[i][0]=tmp.real();
                fftw.in_spm[i][1]=tmp.imag();
  
                }
             
               fftw_execute_dft(fftw.t2f,&fftw.in[0],&fftw.out[0]);
               fftw.shift(&fftw.out[0],&fftw.out[N_half-1],&fftw.out[N_half]);
               
                fftw_execute_dft(fftw.t2f,&fftw.in_spm[0],&fftw.out_spm[0]);
                fftw.shift(&fftw.out_spm[0],&fftw.out_spm[N_half-1],&fftw.out_spm[N_half]);

                  #pragma omp parallel for private(tmp)
                for(int i=0;i<N_f;i++)
                {
                    tmp.real(fftw.out[i][0]);
                    tmp.imag(fftw.out[i][1]);
                    P_ir2[i]=-6.0*dt*chi.sigma3_iee_b2[i]*tmp*std::exp(-C.II*E.k_ir2[i]*z)/(C.eps*2*C.c*chi.n_p2[i]);
                    
                    tmp.real(fftw.out_spm[i][0]);
                    tmp.imag(fftw.out_spm[i][1]);
                    P_ir2[i]+=(-3.0*dt*chi.sigma3_spm2[i]*tmp*std::exp(C.II*E.k_ir2[i]*z)/(C.eps*2*C.c*chi.n_p2[i]));
                    
                }

       //------------------------------------
       //THz generation 2w2-w1=thz
       //------------------------------------     
               #pragma omp parallel for private(tmp)
              for(int i=0;i<N_f;i++)
               {
       
                tmp=std::conj(fftw.complex(fftw.E_t_ir1[i]))*std::pow(fftw.complex(fftw.E_t_ir2[i]),2)*std::pow(df,3);
                fftw.in[i][0]=tmp.real();
                fftw.in[i][1]=tmp.imag(); 
                
                 
                }
             
               fftw_execute_dft(fftw.t2f,&fftw.in[0],&fftw.out[0]);
               fftw.shift(&fftw.out[0],&fftw.out[N_half-1],&fftw.out[N_half]);

                #pragma omp parallel for private(tmp)
               for(int i=0;i<N_f_thz;i++)
               {
                   tmp.real(fftw.out[i][0]);
                   tmp.imag(fftw.out[i][1]);
                   P_thz[i]=-3.0*dt*chi.sigma3_iee[i]*tmp*std::exp(-C.II*E.k_thz[i]*z)/(C.eps*2*C.c*chi.n_thz[i]);
           
                }
                

    
};
#endif
