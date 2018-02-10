//Class for Schrodinger equation solver
//# include <math.h>
//# include <stdlib.h>
//# include <stdio.h>
# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
//# include <cmath>
# include <vector>
# include "NLSolverClass.hpp"


Eigenvalues::~Eigenvalues() {}

/*=======================================================================
 * HO Eigenvalues generator
 *=====================================================================*/

HarmonicEigenvalues::HarmonicEigenvalues(double omega, unsigned int n, int l): m_omega(omega), m_n(n), m_l(l) {}
double HarmonicEigenvalues::eigenvalue() const
{
    return Parameters::hbar*m_omega*(2*m_n+m_l+(3/2));
}
Eigenvalues* HarmonicEigenvalues::clone() const
{
    return new HarmonicEigenvalues(*this);
}


/*========================================================================
 * Eigenvalues generator by shooting method
 *======================================================================*/


TrialEigenvalues* TrialEigenvalues::m_trialEigenvalusObj = nullptr;
TrialEigenvalues::TrialEigenvalues(): m_eigenval1(200), m_eigenval2(400) {}

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

GenericEigenvalues::GenericEigenvalues (const InitialPot& initPot, double H): m_initPot(initPot.clone()), m_h(H) {}

double GenericEigenvalues::eigenvalue() const //this must return a double..
{
	const double error = 10e-8;
	Schroddy sch(*m_initPot, m_h);

	NLSolver <schroddywrapper, &schroddywrapper::eigenfunction> sol(error);

	const schroddywrapper wrap(sch);
    return sol.solveByBisection(wrap, 0, TrialEigenvalues::getEigenval1() , TrialEigenvalues::getEigenval2());

}

Eigenvalues* GenericEigenvalues::clone() const
{
	return new GenericEigenvalues(*this);
}

/*=========================================================================
 * Schrodinger equation solver by Runge-Kutta method
 *=======================================================================*/

Schroddy::Schroddy(const InitialPot& pot, double H/*, const Eigenvalues& eigenval*/): m_pot(pot.clone()), m_h(H)/*, m_eigenval(eigenval.clone())*/ {}

Schroddy::Schroddy(const Schroddy& sourceSchroddy) //copy constructor
{
    if (&sourceSchroddy != this)
    {
    	m_pot = sourceSchroddy.m_pot -> clone();
    	m_h = sourceSchroddy.m_h;
    }
}

Schroddy& Schroddy::operator=(const Schroddy& rhsSchroddy) //copy assignment
{
    if (&rhsSchroddy != this)
    {
        delete m_pot;
        m_pot = rhsSchroddy.m_pot -> clone();
        m_h = rhsSchroddy.m_h;
    }
    return *this;
}

Schroddy::~Schroddy()
{
    delete m_pot;
}

double Schroddy::solveShroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E) const
{
    const double h = m_h;
    const unsigned long NStep = (x1 - x0)/h;
    const double factor = 2*Parameters::mn/(Parameters::hbarc*Parameters::hbarc);
    const double eigenvalue = E;//m_eigenval -> eigenvalue();

    double runningX = x0, runningPsi = psi0, runningPsiPrime = psiPrime0;
    std::vector<double> psiArray;
    psiArray.push_back(psi0);

    for (unsigned long i = 0; i < NStep ; ++i)
    {

        //compute decoupled RK factors
        const double k1 = runningPsiPrime;
        const double l1 = (m_pot -> potential(runningX) - eigenvalue)*runningPsi;
        const double k2 = runningPsiPrime + h/2.*l1;
        const double l2 = (m_pot -> potential(runningX + h/2.) - eigenvalue)*(runningPsi + h/2.*k1);
        const double k3 = runningPsiPrime + h/2.*l2;
        const double l3 = (m_pot -> potential(runningX + h/2.) - eigenvalue)*(runningPsi + h/2.*k2);
        const double k4 = runningPsiPrime + h*l3;
        const double l4 = (m_pot -> potential(runningX + h) - eigenvalue)*(runningPsi + h*k3);

        //advance running variables and store intermediate results
        runningPsi += h/6.*factor*(k1 + 2*k2 + 2*k3 + k4);
        runningPsiPrime += h/6.*factor*(l1 + 2*l2 + 2*l3 + l4);
        runningX += h;
        psiArray.push_back(runningPsi);

    }
    //return runningPsi/x1;

    //work out normalizatiion constant
    double psiSquared = 0;// normalPsi=0;
    for (std::vector<double>::iterator it = psiArray.begin(); it != psiArray.end(); ++it)
        psiSquared += (*it)*(*it); //derefence iterator at array's element and square (brackets are needed!)

    const double scalar = psiSquared*h;

    //return normalized eigenfunction (on interval [x0, x1]) and job done!
    const double normalPsi = (runningPsi)/sqrt(scalar);
    return normalPsi/x1;
}

/*========================================================================
 * Wrapper for resolve arguments inconsistency between schroddy and NLSolver classes
 *======================================================================*/

//schroddywrapper::schroddywrapper (const Schroddy& sh, double x0, double x1, double psi0, double psiPrime0, unsigned long NSteps): m_sh(sh), m_x0(x0), m_x1(x1), m_psi0(psi0), m_psiPrime0(psiPrime0), m_NSteps(NSteps) {}
schroddywrapper::schroddywrapper (const Schroddy& sh): m_sh(sh) {}

double schroddywrapper::eigenfunction(double E) const
{
    return m_sh.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, E);
}


/*
    double integral(double(*fun)(double), double xmin, double xmax, int n_int );
    // Eigenfunction normalization by trapezes method
    double h_t=0, scalarPsi=0, normalPsi=0;
    double integralPsi=runningPsi*runningPsi;
    h_t=(Parameters::x_min-Parameters::x_max)/NSteps;
    double integral(double(*fun)(double), double xmin, double xmax, int n_int )
    {
    	double hint=0, i=0, value=0;
    	hint=(xmin-xmax)/n_int;
    	for (i = xmin; i < xmax; i+=n_int)
    		value+=(((*fun)(i)+(*fun)(i+hint))*hint/2)
    	return value;
    }
    scalarPsi=integral(integralPsi, Parameters::x_min, Parameters::x_max, h_t)
    normalPsi=runningPsi/sqrt(scalarPsi);
}*/
