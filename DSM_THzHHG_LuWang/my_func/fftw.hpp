#ifndef MY_FFTW_H
#define MY_FFTW_H
#include <algorithm>                                                                   //include swap function
#include<complex.h>
#include<fftw3.h>

//**************************************************
//fftw initialization
// %matlab fft from time to frequency is int exp(-iwt) dt
// this is for exp(iwt-ikz), thus the frequency domian is the possitive
//frequency range
//************************************************
//fftw normalisation , overall operation(fft and ifft back) should be normalized by 1/N_f   (df*dt);
//fftw spectrum range [0 fs/2 -fs/2 0]
//in this code the propagation factor is  exp(-iwt+ikz)
//ifft exp(iwt) from time to frequency normalized by *df
//fft  exp(-iwt) from frequency to time normalized by *dt
// in this way the calculation falls in the possitive frequency domain.
//*******************************************************************

struct myfft{

fftw_plan f2t=nullptr;
fftw_plan t2f=nullptr;
fftw_plan c2r=nullptr;
fftw_plan r2c=nullptr;

int N_f;
int N_f_thz;
int n_thread;
int n_nested;
fftw_complex  *in,*out,*E_t_ir1,*E_ir1,*E_thz,*E_t_ir2,*E_ir2,*E_thzcr,*in_spm,*out_spm;
double *E_t_thz,*in_r;
std::complex<double> tmp;
myfft(const int& N_f_,const int&N_f_thz_,const int& n_thread_,const int& n_nested_): N_f(N_f_),N_f_thz(N_f_thz_),n_thread(n_thread_),n_nested(n_nested_){}
~myfft(){};
void my_fftw_clean();
void assign();
// template<class T1>
// void add_conj_flip(T1& E_in ,T1& E_revers,int& N_c);
template<class T1>
void shift(T1 first1, T1 last1,T1 first2);
std::complex<double> complex(fftw_complex A);
}; 


void myfft::assign(){

E_t_ir1=fftw_alloc_complex(N_f*n_thread);
E_t_ir2=fftw_alloc_complex(N_f*n_thread);
//E_t_thz=fftw_alloc_complex(N_f*n_thread);
E_t_thz=(double*)fftw_malloc(sizeof(double) * N_f);
in_r=(double*)fftw_malloc(sizeof(double) * N_f);
out= fftw_alloc_complex(N_f*n_thread);
in=fftw_alloc_complex(N_f*n_thread);
in_spm=fftw_alloc_complex(N_f*n_thread);
out_spm=fftw_alloc_complex(N_f*n_thread);
E_ir1=fftw_alloc_complex(N_f*n_thread);
E_ir2=fftw_alloc_complex(N_f*n_thread);
// zero pad one more element for complex to real Fourier transform
// E_thz negative frequency part is complex conjugate of itself, so E_t_thz always real;
//E_thz=fftw_alloc_complex((N_f_thz)*n_thread);
E_thzcr=fftw_alloc_complex((N_f/2+1)*n_thread);
    
    
t2f=fftw_plan_dft_1d(N_f,fftw_alloc_complex(N_f),fftw_alloc_complex(N_f),1, FFTW_MEASURE);
//the c to r change extend half the frequency array to full frequency array, thus the time 
//domian energy equals to frequency domian with both possitive and negative w
c2r= fftw_plan_dft_c2r_1d(N_f,fftw_alloc_complex(N_f_thz+1),fftw_alloc_real(N_f),FFTW_MEASURE);
r2c=fftw_plan_dft_r2c_1d(N_f,fftw_alloc_real(N_f),fftw_alloc_complex(N_f_thz+1),FFTW_MEASURE);
f2t=fftw_plan_dft_1d(N_f,fftw_alloc_complex(N_f),fftw_alloc_complex(N_f),-1,FFTW_MEASURE);

}



//*************************************************
//only for even input array
//**************************************************
template<class T1>
  void myfft::shift(T1 first1, T1 last1,T1 first2)
{
  while (first1<=last1) {
    std::swap(*first1, *first2);
    ++first1; ++first2;
  }
  
};

// template<class T1>
//   void myfft::add_conj_flip(T1& E_in ,T1& E_reverse,int& N_c)
// {
//     E_reverse=E_in.reverse().conjugate();
//     E_in.head(N_f-2*N_c)+=E_reverse.tail(N_f-2*N_c);
//   
// };
std::complex<double> myfft::complex(fftw_complex A){
    tmp.real(A[0]);
    tmp.imag(A[1]);
    return tmp;
};

void myfft::my_fftw_clean()
{
//fftw_cleanup_threads();
//fftw_mpi_cleanup();                                                           
fftw_destroy_plan(f2t);
fftw_destroy_plan(t2f);
fftw_destroy_plan(c2r);
fftw_destroy_plan(r2c);
fftw_free(E_ir1);
fftw_free(E_ir2);
fftw_free(in_spm);
fftw_free(out_spm);
fftw_free(in_r);
fftw_free(E_thzcr);
fftw_free(in);
fftw_free(E_t_thz);
fftw_free(out);
fftw_free(E_t_ir1);
fftw_free(E_t_ir2);
fftw_cleanup();
}




#endif
