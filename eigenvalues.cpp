#include "eigenvalues.hpp"
#include "parameters.h"
//#include "schroddy.h"

Eigenvalues::~Eigenvalues() {}

/*=======================================================================
 * HO Eigenvalues generator
 *=====================================================================*/

HarmonicEigenvalues::HarmonicEigenvalues(unsigned int nr, int l): m_nr(nr), m_l(l) {}
double HarmonicEigenvalues::eigenvalue() const
{
    return Parameters::hbar_omega*(2*(m_nr-1)+m_l+(3/2));
}

std::unique_ptr<Eigenvalues> HarmonicEigenvalues::clone() const
{
    return std::make_unique<HarmonicEigenvalues>(*this);
}


/*========================================================================
 * Eigenvalues generator by shooting method
 *======================================================================*/


TrialEigenvalues* TrialEigenvalues::m_trialEigenvalusObj = nullptr;
TrialEigenvalues::TrialEigenvalues(): m_eigenval1(-50), m_eigenval2(100) {}

double TrialEigenvalues::getEigenval1()
{
    if (m_trialEigenvalusObj == nullptr)
        m_trialEigenvalusObj = new TrialEigenvalues; //lazy instantiation

    return m_trialEigenvalusObj -> m_eigenval1;
}

double TrialEigenvalues::getEigenval2()
{
    if (m_trialEigenvalusObj == nullptr)
        m_trialEigenvalusObj = new TrialEigenvalues; //lazy instantiation

    return m_trialEigenvalusObj -> m_eigenval2;
}

GenericEigenvalues::GenericEigenvalues (const Schroddy& sh, unsigned int nState, unsigned int lState):
m_sh(sh), m_nState(nState), m_lState(lState) {}

// Shooting method implementation
double GenericEigenvalues::_shootingMethod(double E1, double E2, unsigned int nState) const
{

    const double error = 1e-8;
    const unsigned int lState = m_lState;

    std::vector <double> psiArray;
    unsigned int nodes = 0;
    while (true)
    {

        m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(lState),
                                          psiPrime0(lState), E2 , psiArray);

        nodes = 0;
        for (unsigned long i = 0; i < psiArray.size()-1; ++i)
        {
            if (psiArray[i]*psiArray[i+1] < 0)
                ++ nodes;
        }

        if (nodes > nState)
            E2 -= 0.5;
        else if (nodes < nState)
            E2 += 0.5;
        else break;
    }

    int endPointSign = 1;
    if (*(psiArray.end() -1) < 0)
        endPointSign = -1;

    double midE = (E1 + E2)/2.;
    while (std::abs(E2 - E1) > error)
    {
        double sMid = m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(lState),
                                             psiPrime0(lState), midE , psiArray);

        if (endPointSign*sMid > 0)
            E2 = midE;
        else
            E1 = midE;

        midE = (E1 + E2)/2.;
    }

    return midE;
}

// level energy
double GenericEigenvalues::eigenvalue() const //this must return a double..
{
    m_ELev.clear();
    double E1 = TrialEigenvalues::getEigenval1();
    for (unsigned int i = 1; i <= m_nState; ++i)
    {
        E1 = _shootingMethod(E1, TrialEigenvalues::getEigenval2(), i);
        m_ELev[i] = E1;
    }
    return E1;

}

EnergyLevels GenericEigenvalues::getAllLevels() const
{
    return m_ELev;
}

std::unique_ptr<Eigenvalues> GenericEigenvalues::clone() const
{
    return std::make_unique<GenericEigenvalues>(*this);
}
