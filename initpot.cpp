//Class for initial potentials
#include <cmath>
#include "parameters.h"
#include "initpot.h"

InitialPot::~InitialPot() {}
InitialPot::InitialPot(double m, unsigned int l): m_m(m), m_anglmomentum(l) {};

/*===================================================================
 * Woods-Saxon potential
 *=================================================================*/

WSaxPot::WSaxPot(double Rn, double a0, double m): InitialPot(m, 0), m_Rn(Rn), m_a0(a0) {};
WSaxPot::WSaxPot(double Rn, double a0, double m, unsigned int l): InitialPot(m, l), m_Rn(Rn), m_a0(a0) {};

double WSaxPot::potential(double x) const
{
    /*
    int s = 2.;
    const double Ru = Parameters::Rn*sqrt((1 + (5*s*s)/(2*Parameters::Rn*Parameters::Rn))/
    		(1 + (3*s*s)/(4*Parameters::Rn*Parameters::Rn)));

    double copot;
    if (x <= Ru)
    	copot = 0.5*((Parameters::NP*Parameters::qe*Parameters::qe)/Ru)*(3 - (x/Ru)*(x/Ru));
    else
    	copot = (Parameters::NP*Parameters::qe*Parameters::qe)/x;

    */int sign = -1;
    if (Parameters::NN == 0)
        sign *= -1;
    
    double V0 = 51 + sign*33.0*(Parameters::NN - Parameters::NP)/Parameters::A;
    double wspot = -V0/(1+exp ((x-m_Rn)/m_a0));// + copot;
    return wspot;
}


std::unique_ptr<InitialPot> WSaxPot::clone() const
{
    return std::make_unique<WSaxPot> (*this); //return a derived class object through a base class pointer
}


/*===================================================================
 * HO potential
 *=================================================================*/

HOPot::HOPot(double m, unsigned int l): InitialPot(m, l) {}
HOPot::HOPot(double m): InitialPot(m, 0) {}

double HOPot::potential(double x) const
{
    using namespace Parameters;
    double coulomb = 0;
    /*
    if (x != 0)
    {
        angularPart = (hbarc*hbarc)*m_anglmomentum*(m_anglmomentum+1)/(2*m_m*x*x);
        //coulomb = -NN*qe*qe/x;
    }*/

    const double c =(m_m*hbar_omega*hbar_omega)/(hbarc*hbarc);
    double hopot = 0.5*c*x*x + coulomb;

    return hopot;
}


std::unique_ptr<InitialPot> HOPot::clone() const
{
    return std::make_unique<HOPot> (*this); //return a derived class object through a base class pointer
}

/*===================================================================
 * Coulomb potential
 *=================================================================*/

CoPot::CoPot(int Z, double m): InitialPot(m, 0), m_Z(Z) {}
CoPot::CoPot(int Z, double m, unsigned int l): InitialPot(m, l), m_Z(Z) {}

double CoPot::potential(double x) const
{
    int s = 2.;
    const double Ru = Parameters::Rn*sqrt((1 + (5*s*s)/(2*Parameters::Rn*Parameters::Rn))/
    		(1 + (3*s*s)/(4*Parameters::Rn*Parameters::Rn)));

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

/*===================================================================
 * Total potential
 *=================================================================*/

/*TotPot& operator+(const InitialPot&);
{
	return TotPot()
}

TotPot::TotPot(double m): InitialPot(m, 0) {}
TotPot::TotPot(double m, unsigned int l): InitialPot(m, l) {}

double TotPot::potential(double x) const
{
	const WSaxPot wspot(Parameters::Rn, Parameters::a0, Parameters::mn);
	const CoPot copot(Parameters::NP, Parameters::mn);
	double totpot = wspot + copot;

    return totpot;
}

std::unique_ptr<InitialPot> TotPot::clone() const
{
    return std::make_unique<TotPot> (*this); //return a derived class object through a base class pointer
}*/

/*=====================================================================
 * Kohn-Sham potential
 *===================================================================*/

PotOut::PotOut(const KohnShamInverse& outpot, double m, unsigned int l): InitialPot(m, l), m_outpot(outpot) {}
PotOut::PotOut(const KohnShamInverse& outpot, double m): InitialPot(m, 0), m_outpot(outpot) {}

double PotOut::potential(double x) const
{
    KSPotential ksp = m_outpot.getKSPot();
    if (ksp.find(x) == ksp.end())
        return interpolatedPotential(x);
    
    return m_outpot.getKSPot().at(x);
}

std::unique_ptr<InitialPot> PotOut::clone() const
{
	 return std::make_unique<PotOut> (*this); //return a derived class object through a base class pointer
}

/* Linear interpolator/extrapolator for solve R-K method with a map (discrete values) as input
 * rather than a continues function.*/
double PotOut::interpolatedPotential(double x) const
{
    KSPotential ksp = m_outpot.getKSPot(); 	// if x < first element in map,
    if (x < ksp.begin() -> first)			// set potential point value at x equal to the first element
        return ksp.begin() -> second; //lower extrapolation
    
    KSPotential::iterator it = ksp.begin();
    KSPotential::iterator p = it;
    ++it;
    
    for (; it != ksp.end(); ++it)					// if x is between two consecutive elements of
    {												// the map, add a x value between this range
        if ( x >= p -> first && x < it -> first)	// using straight line equation
            return p -> second +
            (it -> second - p -> second)/(it -> first - p -> first)*(x - p -> first); //interpolation
        ++p;
    }
    
    return p -> second; //at this level p points to last map element, which implies x > lastMapVaulue; upper extrapolation
}

/*=====================================================================
 * Kohn-Sham potential from J-W method
 *===================================================================*/

PotOutJW::PotOutJW(const KohnShamInverseWithJW& outpot, double m, unsigned int l): InitialPot(m, l), m_outpot(outpot) {}
PotOutJW::PotOutJW(const KohnShamInverseWithJW& outpot, double m): InitialPot(m, 0), m_outpot(outpot) {}

double PotOutJW::potential(double x) const
{
    JWKSPotential ksp = m_outpot.getJWKSPot();
    if (ksp.find(x) == ksp.end())
        return interpolatedPotential(x);

    return m_outpot.getJWKSPot().at(x);
}

std::unique_ptr<InitialPot> PotOutJW::clone() const
{
	 return std::make_unique<PotOutJW> (*this); //return a derived class object through a base class pointer
}

/* Linear interpolator/extrapolator for solve R-K method with a map (discrete values) as input
 * rather than a continues function.*/
double PotOutJW::interpolatedPotential(double x) const
{
    JWKSPotential ksp = m_outpot.getJWKSPot(); 	// if x < first element in map,
    if (x < ksp.begin() -> first)			// set potential point value at x equal to the first element
        return ksp.begin() -> second; //lower extrapolation

    JWKSPotential::iterator it = ksp.begin();
    JWKSPotential::iterator p = it;
    ++it;

    for (; it != ksp.end(); ++it)					// if x is between two consecutive elements of
    {												// the map, add a x value between this range
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

TestPot::TestPot(double m, unsigned int l): InitialPot(m, l) {}
TestPot::TestPot(double m): InitialPot(m, 0) {}

double TestPot::potential(double x) const
{

	using namespace Parameters;
	const double c =(m_m*hbar_omega*hbar_omega)/(hbarc*hbarc);
    const double perturbativePart = 10*x;
    const int verTraslation = -50;
    double hopot = 0.5*c*x*x + perturbativePart /*+ verTraslation*/;


    return hopot;
}


std::unique_ptr<InitialPot> TestPot::clone() const
{
    return std::make_unique<TestPot> (*this); //return a derived class object through a base class pointer
}





