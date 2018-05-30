//Class for initial potentials
#include <cmath>
#include <cstdlib>
#include "parameters.h"
#include "initpot.h"

InitialPot::~InitialPot() {}
InitialPot::InitialPot(double m, unsigned int l, double j): m_m(m), m_anglmomentum(l), m_spin(j) {};

/*===================================================================
 * Woods-Saxon potential
 *=================================================================*/

WSaxPot::WSaxPot(double Rn, double a0, double m): InitialPot(m, 0, 0), m_Rn(Rn), m_a0(a0) {};
WSaxPot::WSaxPot(double Rn, double a0, double m, unsigned int l): InitialPot(m, l, 0), m_Rn(Rn), m_a0(a0) {};
WSaxPot::WSaxPot(double Rn, double a0, double m, unsigned int l, double j): InitialPot(m, l, j), m_Rn(Rn), m_a0(a0) {};

double WSaxPot::potential(double x) const
{
    using namespace Parameters;
    /*int s = 2.;
    const double Ru = Parameters::Rn*sqrt((1 + (5*s*s)/(2*Parameters::Rn*Parameters::Rn))/
    		(1 + (3*s*s)/(4*Parameters::Rn*Parameters::Rn)));

    double copot;
    if (x <= Ru)
    	copot = 0.5*((Parameters::NP*Parameters::qe*Parameters::qe)/Ru)*(3 - (x/Ru)*(x/Ru));
    else
    	copot = (Parameters::NP*Parameters::qe*Parameters::qe)/x;*/

    int sign = -1;
    if (Parameters::nucleons == 0)
        sign = 1;
    
    double V0 = 51 + sign*33.0*(ElementConstants::NN() - ElementConstants::NP())/ElementConstants::A();
    double wspot = -V0/(1+exp ((x-m_Rn)/m_a0)) /*+ copot*/;
    return wspot;
}


std::unique_ptr<InitialPot> WSaxPot::clone() const
{
    return std::make_unique<WSaxPot> (*this); //return a derived class object through a base class pointer
}


/*===================================================================
 * HO potential
 *=================================================================*/

HOPot::HOPot(double m, unsigned int l, double j): InitialPot(m, l, j) {}
HOPot::HOPot(double m, unsigned int l): InitialPot(m, l, 0) {}
HOPot::HOPot(double m): InitialPot(m, 0, 0) {}

double HOPot::potential(double x) const
{
    using namespace Parameters;
    const double c =(m_m*ElementConstants::hBarOmega()*ElementConstants::hBarOmega())/(hbarc*hbarc);
    double hopot = 0.5*c*x*x;

    return hopot;
}


std::unique_ptr<InitialPot> HOPot::clone() const
{
    return std::make_unique<HOPot> (*this); //return a derived class object through a base class pointer
}

/*===================================================================
 * Coulomb potential
 *=================================================================*/

CoPot::CoPot(int Z, double m): InitialPot(m, 0, 0), m_Z(Z) {}
CoPot::CoPot(int Z, double m, unsigned int l): InitialPot(m, l, 0), m_Z(Z) {}
CoPot::CoPot(int Z, double m, unsigned int l, double j): InitialPot(m, l, j), m_Z(Z) {}

double CoPot::potential(double x) const
{
    using namespace Parameters;
    
    int s = 2.;
    const double Ru = ElementConstants::Rn()*sqrt((1 + (5*s*s)/(2*ElementConstants::Rn()*ElementConstants::Rn()))/
    		(1 + (3*s*s)/(4*ElementConstants::Rn()*ElementConstants::Rn())));

    double copot;
    if (x <= Ru)
    	copot = 0.5*((m_Z*Parameters::qe*Parameters::qe)/Ru)*(3 - (x/Ru)*(x/Ru));
    else
    	copot = (m_Z*Parameters::qe*Parameters::qe)/x;

    return copot;
}

std::unique_ptr<InitialPot> CoPot::clone() const
{
    return std::make_unique<CoPot> (*this); //return a derived class object through a base class pointer
}

/*=====================================================================
 * Kohn-Sham potential
 *===================================================================*/

PotOut::PotOut(const KohnShamInverse& outpot, double m, unsigned int l, double j): InitialPot(m, l, j), m_outpot(outpot.getKSPot()) {}
PotOut::PotOut(const KohnShamInverse& outpot, double m, unsigned int l): InitialPot(m, l, 0), m_outpot(outpot.getKSPot()) {}
PotOut::PotOut(const KohnShamInverse& outpot, double m): InitialPot(m, 0, 0), m_outpot(outpot.getKSPot()) {}
PotOut::PotOut(const KSPotential& inPot, double m): InitialPot(m, 0, 0), m_outpot(inPot) {}

double PotOut::potential(double x) const
{
    KSPotential ksp = m_outpot;
    if (ksp.find(x) == ksp.end())
        return interpolatedPotential(x);
    
    return m_outpot.at(x);
}

std::unique_ptr<InitialPot> PotOut::clone() const
{
	 return std::make_unique<PotOut> (*this); //return a derived class object through a base class pointer
}

/* Linear interpolator/extrapolator to enable integration with R-K method on a map (discrete values) as an input
 * rather than a continuous function.*/
double PotOut::interpolatedPotential(double x) const
{
    KSPotential ksp = m_outpot; 	// if x < first element in map,
    if (x < ksp.begin() -> first)			// set potential point value at x equal to the first element
        return ksp.begin() -> second; //lower extrapolation
    
    KSPotential::iterator it = ksp.begin();
    KSPotential::iterator p = it;
    ++it;
    
    for (; it != ksp.end(); ++it)					// if x is between two consecutive elements of
    {												// the map, add x value in this range
        if ( x >= p -> first && x < it -> first)	// using straight line equation
            return p -> second +
            (it -> second - p -> second)/(it -> first - p -> first)*(x - p -> first); //interpolation
        ++p;
    }
    
    return p -> second; //at this level p points to last map element, which implies x > lastMapVaulue; upper extrapolation
}

/*===================================================================
 * Test potential
 *=================================================================*/

TestPot::TestPot(double m, unsigned int l, double j): InitialPot(m, l, j) {}
TestPot::TestPot(double m, unsigned int l): InitialPot(m, l, 0) {}
TestPot::TestPot(double m): InitialPot(m, 0, 0) {}

double TestPot::potential(double x) const
{

	using namespace Parameters;
    const double c =(m_m*ElementConstants::hBarOmega()*ElementConstants::hBarOmega())/(hbarc*hbarc);
	const double verTraslation = -50.;
    //const double perturbativePart = 10*x;
    const double perturbativePart2 = (rand()/(static_cast<double>(RAND_MAX))) - 0.05;
    //double hopot = 0.05*x*x + perturbativePart2 /*+ verTraslation*/;
    double hopot = 0.05*c*x*x + perturbativePart2 + verTraslation;


    return hopot;
}


std::unique_ptr<InitialPot> TestPot::clone() const
{
    return std::make_unique<TestPot> (*this); //return a derived class object through a base class pointer
}





