//Parameters
#pragma once
//#ifndef parameters_h
//#define parameters_h
#include<cmath>



namespace Parameters
{

const int NN=8;										// Neutrons number
const int NP=0;										// Protons number
//const double mp= 1.6726219e-27;					// Proton mass
//const double mn= 1.6749273e-27;					// Neutron mass
const double mn= 939.565378;                        // Neutron mass in mnc^2 [MeV]
    const double me = 0.510;
const int A=NN+NP; 									// Mass number
const double R0= 1.27; 								// [fm]

// Wood-Saxon potential parameters
const double Rn = R0*pow(A,(1/3)); 					// Nuclear radius [fm]
const double a0 = 0.67;

// Spin-Orbit potential parameters
//const double k0=  ;
//const double r0=  ;

// HO potential parameters
//const double f = 2.417988e21;						// [Hz]
    //const double k = 1e-24;


// Eigenvalues generator parameters
/*double eigenvalue1=0.93250 ;
double eigenvalue2=1.36256 ;
double error=10e-8;*/

// Runge-Kutta parameters

const double x_in = -10*Rn;
const double x_fin = 10*Rn;
const double psi0= -0.0001;
//const double Psi_0=;
const double psiPrime0 = 0;

//unsigned long N_step=10;

// Other parameters
//const unsigned int energyLevel = 1;
//const int angularMomentum = 2;
    const double PI = 4*atan(1);
const double hbar = 6.58211928e-16; 				// Reduced constant Planck [eV*s]
const double hbarc= 197.3269788;					// [MeV*fm]
//const double x_min=0;								// Integration
//const double x_max=3*Rn;							// Interval
    
    const double hbar_omega=10;							// [MeV]

}

//#endif

