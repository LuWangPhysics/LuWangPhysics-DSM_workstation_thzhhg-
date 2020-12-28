#pragma once
#include <iostream>
#include <math.h>
#include <complex>
struct  CONST
{      
            const double pi=M_PI;   
            const double c=3e8;                                                                  //Speed of light  
            const double eps =8.85e-12;                                                         //Free space permittivity
            const double kb=1.380649e-23;                                                        //boltzmann constant J/k
            std::complex<double> II;                                                              //imaginary unite
            const double e=1.60217662e-19;                                                     //elementary charge
            const double hbar=1.0545718e-34;
            std::string name_tail;
	CONST(){
    using namespace std::complex_literals;
       II=1i;
    };
    
	~CONST(){} ;
 	
};
