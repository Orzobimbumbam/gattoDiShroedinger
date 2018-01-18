//Parameters
#pragma once
#include<cmath>

namespace Parameters
{

const double PI = 4*atan(1);

// Nucleons parameters
const unsigned int N=10;									// Neutrons number
const unsigned int P=1;									// Protons number
//const double mn=  ;							// Neutron mass
//const double mp=   ;						// Proton mass

// Wood-Saxon potential parameters
const double V0=25.025;			
//double Rn=  ;
//const double a0=   ;

// Spin-Orbit potential parameters
//const double k0=  ;
//const double r0=  ;

// HBO potential parameters
//const double m=   ;
double f = 100;
double omega=2*PI*f;
//double r=   ;
const  unsigned int n = 100;
const unsigned int nl=2;

// Eigenvalues generator parameters
//double E1= ;
//double E2= ;

// Other parameters
const double hbar = 1.0545718*1e-34;
}
