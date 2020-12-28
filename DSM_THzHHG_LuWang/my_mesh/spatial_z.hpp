#pragma once

struct  spacial_z
{      

	double L,L_save;
        double dz;
        int N_z,N_total; 
        Eigen::ArrayXd z_0;
        double z_p; // the z at each iteration point
	int N_z_iter;
        int N_z_save;
	

	spacial_z(const double& L_,const double& Lsave_, const double& dz_, const int&N_z_): L(L_),L_save(Lsave_), dz(dz_),N_z(N_z_)
        {   
         
           z_0=Eigen::ArrayXd::LinSpaced(N_z+1,0,L); 
           N_z_iter=0;
           N_z_save=(L/dz/N_z)+((int)(L/dz)%N_z!=0);
           N_total=L/dz;

            
        };
	~spacial_z(){} ;
 	
};
