//Class for initial potentials
#include <cmath>
#include "parameters.h"
#include "initpot.h"

InitialPot::~InitialPot() {}
InitialPot::InitialPot(double m, unsigned int l): m_m(m), m_anglmomentum(l) {};

/*===================================================================
 * Woods-Saxon potential
 *=================================================================*/

WSaxPot::WSaxPot(double V0, double Rn, double a0, double m, unsigned int l): InitialPot(m, l), m_V0(V0), m_Rn(Rn), m_a0(a0) {};

double WSaxPot::potential(double x) const
{
    double angularPart = 0;
    if (x != 0)
        angularPart = (Parameters::hbarc*Parameters::hbarc)*m_anglmomentum*(m_anglmomentum+1)/(2*m_m*x*x);
    
    double wspot = -m_V0/(1+exp ((x-m_Rn)/m_a0)) + angularPart;
    return wspot;
}

void WSaxPot::setL(unsigned int l)
{
    m_anglmomentum = l;
}

std::unique_ptr<InitialPot> WSaxPot::clone() const
{
    return std::make_unique<WSaxPot> (*this); //return a derived class object through a base class pointer
}


/*===================================================================
 * HO potential
 *=================================================================*/

HOPot::HOPot(double m, unsigned int l): InitialPot(m, l), m_m(m) {}

void HOPot::setL(unsigned int l)
{
    m_anglmomentum = l;
}


double HOPot::potential(double x) const
{
    double angularPart = 0;
    if (x != 0)
        angularPart = (Parameters::hbarc*Parameters::hbarc)*m_anglmomentum*(m_anglmomentum+1)/(2*m_m*x*x);
    
    const double c =(m_m*Parameters::hbar_omega*Parameters::hbar_omega)/(Parameters::hbarc*Parameters::hbarc);
    double hopot = angularPart + 0.5*c*x*x;

    return hopot;
}
/*
double HOPot::getOmega() const
{
    return m_omega;
}*/

std::unique_ptr<InitialPot> HOPot::clone() const
{
    return std::make_unique<HOPot> (*this); //return a derived class object through a base class pointer
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

/*=====================================================================
 * Kohn-Sham potential
 *===================================================================*/

potOut::potOut(KohnShamInverse outpot): m_outpot(outpot){}
double potOut::potential(double x) const
{
	std::map<double, double> potMap;
	m_outpot.getOutPot(potMap);

	return potMap[x];
}

std::unique_ptr<InitialPot> potOut::clone() const
{
	 return std::make_unique<potOut> (*this); //return a derived class object through a base class pointer
}








