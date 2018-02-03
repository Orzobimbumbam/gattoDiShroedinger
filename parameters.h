//Parameters
#pragma once
//#ifndef parameters_h
//#define parameters_h
#include<cmath>



namespace Parameters
{

const int NN=10;										// Neutrons number
const int NP=10;										// Protons number
const double mp= 1.6726219e-27;							// Proton mass
//const double mn= 1.6749273e-27;							// Neutron mass
const double mn= 1.;
// Wood-Saxon potential parameters
const double V0 = 25.025;
const double Rn = 1e-9;
const double a0 = 1;

// Spin-Orbit potential parameters
//const double k0=  ;
//const double r0=  ;

// HO potential parameters
const double m=mn   ;
const double f = 100;
//double r=   ;

// Eigenvalues generator parameters
/*double eigenvalue1=0.93250 ;
double eigenvalue2=1.36256 ;
double error=10e-8;*/

// Runge-Kutta parameters
const double x_in=0;
const double x_fin=Rn;
const double psi0=0;
const double Psi_0=0;
const double psiPrime0=1;
//unsigned long N_step=10;

// Other parameters
const unsigned int energyLevel = 1;
const int angularMomentum = 3;
const double PI = 4*atan(1);
//const double hbar = 1.0545718e-34;
const double hbar = 1;
const double x_min=0;									//Integration
const double x_max=3*Rn;								//Interval

}

//#endif


