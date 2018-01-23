//Class for initial potentials

#include <cmath>
#include "initpot.h" //remember to include header with interface

InitialPot::~InitialPot() {}

WSaxPot::WSaxPot(double V0, double Rn, double a0): m_V0(V0), m_Rn(Rn), m_a0(a0) {};

double WSaxPot::potential(double x) const
{
    double wspot = -m_V0/(1+exp ((x-m_Rn)/m_a0));
                   
    return wspot;
}

InitialPot* WSaxPot::clone() const
{
    return new WSaxPot(*this); //return a derived class object through a base class pointer
}


HOPot::HOPot(double m, double omega): m_m(m), m_omega(omega) {}

double HOPot::potential(double x) const
{
    double hopot=(1/2)*m_m*(m_omega*m_omega)*(x*x);
    
    return hopot;
}
/*
double HOPot::getOmega() const
{
    return m_omega;
}*/

InitialPot* HOPot::clone() const
{
    return new HOPot(*this); //return a derived class object through a base class pointer
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
