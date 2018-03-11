#include "eigenvalues.hpp"
#include "parameters.h"

Eigenvalues::~Eigenvalues() {}

/*=======================================================================
 * HO Eigenvalues generator
 *=====================================================================*/

HarmonicEigenvalues::HarmonicEigenvalues(unsigned int n, int l): m_n(n), m_l(l) {}
double HarmonicEigenvalues::eigenvalue() const
{
    return Parameters::hbar_omega*(2*(m_n-1)+m_l+(3/2));
}
Eigenvalues* HarmonicEigenvalues::clone() const
{
    return new HarmonicEigenvalues(*this);
}


/*========================================================================
 * Eigenvalues generator by shooting method
 *======================================================================*/


TrialEigenvalues* TrialEigenvalues::m_trialEigenvalusObj = nullptr;
TrialEigenvalues::TrialEigenvalues(): m_eigenval1(0), m_eigenval2(200) {}

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

double GenericEigenvalues::shootingMethod(double E1, double E2, unsigned int nState) const
{
    const double error = 1e-8;
    int parityFlag = 1;
    /*
     if (nState%2 == 0)
     parityFlag = -1;*/

    unsigned int lState = m_lState;


    std::vector <double> psiArray;
    unsigned int nodes = 0;
    while (true)
    {

        double s = m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(lState),
                                          psiPrime0(lState), E2 , psiArray);

        nodes = 0;
        for (unsigned long i = 0; i < psiArray.size()-1; ++i)
        {
            if (psiArray[i]*psiArray[i+1] < 0)
                ++ nodes;
        }

        if (nodes > nState)
            E2 *= 0.99;
        else if (nodes < nState)
            E2 *= 1.11;
        else break;
    }

    int endPointSign = 1;
    if (*(psiArray.end() -1) < 0)
        endPointSign = -1;


    double midE = (E1 + E2)/2.;
    while (std::abs(E2 - E1) > error)
    {
        //double s1 = m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(lState),
        //                                  psiPrime0(lState), E1 , psiArray);
        //double s2 = m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(lState),
        //                                  psiPrime0(lState), E2 , psiArray);


        //if (s1*s2 < 0)
        //{

        double sMid = m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(lState),
                                             psiPrime0(lState), midE , psiArray);

        if (/*s2*/endPointSign*sMid > 0)
            E2 = midE;
        else /*if (s1*sMid < 0)*/
            E1 = midE;
        //}
        /*
         else
         {
         midE = E2;
         break;
         }*/

        midE = (E1 + E2)/2.;
    }

    return midE;
}

double GenericEigenvalues::eigenvalue() const //this must return a double..
{
    double E1 = TrialEigenvalues::getEigenval1();
    for (unsigned int i = 0; i <= m_nState; ++i)
        E1 = shootingMethod(E1, TrialEigenvalues::getEigenval2(), i);

    return E1;

}



Eigenvalues* GenericEigenvalues::clone() const
{
    return new GenericEigenvalues(*this);
}




