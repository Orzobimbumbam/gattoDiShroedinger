//Parameters
#pragma once
//#ifndef parameters_h
//#define parameters_h
#include<cmath>


double psi0(unsigned int l);
double psiPrime0(unsigned int l);


namespace Parameters
{

const int NN = 48;									// Neutrons number
const int NP = 0;									// Protons number
//const double mp= 1.6726219e-27;					// Proton mass [kg]
//const double mn= 1.6749273e-27;					// Neutron mass [kg]
const double mn = 939.565378;                       // Neutron mass in mnc^2 [MeV]
const double mp = 938.28;							// Proton mass in mpc^2 [MeV]
const int A = NN + NP; 								// Mass number
const double R0 = 1.27; 							// [fm]

// Wood-Saxon potential parameters
const double Rn = R0*pow(A,(1/3)); 					// Nuclear radius [fm]
const double a0 = 0.67;								// Nuclear surface thickness [fm]

// Spin-Orbit potential parameters
//const double k0=  ;
//const double r0=  ;

// HO potential parameters
//const double f = 2.417988e21;						// [Hz]
const double hbar_omega = 10;						// [MeV]
//const double k = 1e-24;


// Eigenvalues generator parameters
/*double eigenvalue1=0.93250 ;
double eigenvalue2=1.36256 ;
double error=10e-8;*/

// Runge-Kutta parameters
//const double x_in = -10*Rn;
const double x_in = 1e-12;
const double x_fin = 5*Rn;
//const double psi0= 0.001;
//const double psiPrime0 = 1;
//unsigned long N_step=10;

// Other parameters
//const unsigned int energyLevel = 1;
//const int angularMomentum = 2;
const double qe = 16021e-19;						// elementar charge [C]
const double PI = 4*atan(1);
const double hbar = 6.58211928e-16; 				// Reduced constant Planck [eV*s]
const double hbarc = 197.3269788;					// [MeV*fm]
const double rms = 3.460;							// root-mean-square radius of the charge distribution (<r^2>)^1/2 [fm]
const double rp = 0.8751;							// proton rms charge radius (<rp^2>)^1/2 [fm]
const double pregamma = 0.05;						// gamma prefactor for Kohn-Sham equations
//const double x_min=0;								// Integration
//const double x_max=3*Rn;							// Interval

}

//#endif

