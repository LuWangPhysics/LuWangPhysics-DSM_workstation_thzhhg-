#pragma once
struct  time_frequency
{   int N_f;
    int N_f_thz;
    int N_c;
    double f_max;
    double df;
    double omega_1; 
    double omega_2;
    double dt;
    double f_c_thz;
    double f_thz_3h;
    int N_thz_3h;
    Eigen::ArrayXd omega_thz;
    Eigen::ArrayXd omega_p1;

    Eigen::ArrayXd t;  
    time_frequency(int& N_f_,int& N_f_thz_,double& df_,double & f_thz_3h_): N_f(N_f_),N_f_thz(N_f_thz_),df(df_),f_thz_3h(f_thz_3h_)
        {   double pi = M_PI;
            dt=1/(df*N_f);
            N_c=(int) f_c_thz/df;
            f_c_thz=0e12;
            omega_thz=Eigen::ArrayXd::Zero(N_f);
            t=Eigen::ArrayXd::Zero(N_f,1);


            for(int i=0;i<N_f;i++)
            {

            // 1e-1 df is an offset at 0 frequency to avoid singularity
            omega_thz(i)=2*pi*df/2+2*pi*f_c_thz-2*pi*N_f*df/2+2*pi*df*i;
            t(i)= -(N_f-1)/(2*N_f*df)+i*dt; 
            }
            

        }
	~time_frequency(){} 

};
