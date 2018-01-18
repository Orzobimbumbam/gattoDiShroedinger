//Class for initial potentials

#include <cmath>
#include "initpot.h" //remember to include header with interface


WSaxPot::WSaxPot(double V0, double Rn, double a0): m_V0(V0), m_Rn(Rn), m_a0(a0) {};

double WSaxPot::potential(double x) const
{
    double wspot = -m_V0/(1+exp ((x-m_Rn)/m_a0));
                   
    return wspot;
}



HBOPot::HBOPot(double m, double omega, double rhbo): m_m(m), m_omega(omega) {}

double HBOPot::potential(double x) const
{
    double hbopot=(1/2)*m_m*(m_omega*m_omega)*(x*x);
    
    return hbopot;
}



//=====================================================================
// Spin-orbit potential
//=====================================================================
/*
double SOPot(double k0, doubler0, double x, double hbar, double Rn, double a)
{
	double l, j;
	double LS=((hbar*hbar)/2)(j(j+1)-l(l+1)-3/4)
	
	double deriv 
	
	sopot=*/
