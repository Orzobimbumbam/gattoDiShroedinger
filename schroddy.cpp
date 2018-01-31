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
 * Wrapper for resolve arguments inconsistency between schroddy and NLSolver classes
 *======================================================================*/

//schroddywrapper::schroddywrapper (const Schroddy& sh, double x0, double x1, double psi0, double psiPrime0, unsigned long NSteps): m_sh(sh), m_x0(x0), m_x1(x1), m_psi0(psi0), m_psiPrime0(psiPrime0), m_NSteps(NSteps) {}
schroddywrapper::schroddywrapper (const Schroddy& sh): m_sh(sh) {}

double schroddywrapper::eigenfunction(double E) const
{
	return m_sh.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, E, 10);
}

/*========================================================================
 * Eigenvalues generator by shooting method
 *======================================================================*/


TrialEigenvalues* TrialEigenvalues::m_trialEigenvalusObj = nullptr;
TrialEigenvalues::TrialEigenvalues(): m_eigenval1(0.93250), m_eigenval2(1.36256) {}

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

GenericEigenvalues::GenericEigenvalues (const InitialPot& initPot): m_initPot(initPot.clone()) {}

double GenericEigenvalues::eigenvalue() const //this must return a double..
{
	const double error = 10e-8;
	Schroddy sch(*m_initPot);

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

Schroddy::Schroddy(const InitialPot& pot/*, const Eigenvalues& eigenval*/): m_pot(pot.clone())/*, m_eigenval(eigenval.clone())*/ {}

Schroddy::~Schroddy()
{
    //delete m_pot; this should be fixed by implementing copy constructor
    //delete m_eigenval;
}

double Schroddy::solveShroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, unsigned long NSteps) const
{
    const double h = (x1 - x0)/NSteps;
    const double factor = 2*Parameters::mn/(Parameters::hbar*Parameters::hbar);
    const double eigenvalue = E;//m_eigenval -> eigenvalue();

    double runningX = x0, runningPsi = psi0, runningPsiPrime = psiPrime0;
    std::vector<double> psiArray;
    psiArray.push_back(psi0);

    for (unsigned long i = 0; i < NSteps; ++i)
    {

        //compute decoupled RK factors
        const double k1 = runningPsiPrime;
        const double l1 = (m_pot -> potential(runningX) - eigenvalue)*runningPsi;
        const double k2 = runningPsiPrime + h/2*l1;
        const double l2 = (m_pot -> potential(runningX + h/2) - eigenvalue)*(runningPsi + h/2*k1);
        const double k3 = runningPsiPrime + h/2*l2;
        const double l3 = (m_pot -> potential(runningX + h/2) - eigenvalue)*(runningPsi + h/2*k2);
        const double k4 = runningPsiPrime + h*l3;
        const double l4 = (m_pot -> potential(runningX + h) - eigenvalue)*(runningPsi + h*k3);

        //advance running variables and store intermediate results
        runningPsi += h/6*factor*(k1 + 2*k2 + 2*k3 + k4);
        runningPsiPrime =+ h/6*factor*(l1 + 2*l2 + 2*l3 + l4);
        runningX += h;
        psiArray.push_back(runningPsi);

    }
    //work out normalizatiion constant
    double psiSquared = 0, normalPsi=0;
    for (std::vector<double>::iterator it = psiArray.begin(); it != psiArray.end(); ++it)
        psiSquared += (*it)*(*it); //derefence iterator at array's element and square (brackets are needed!)

    const double scalar = psiSquared*h;

    //return normalized eigenfunction (on interval [x0, x1]) and job done!
    return normalPsi=runningPsi/sqrt(scalar);
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
