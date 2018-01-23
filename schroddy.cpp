//Class for Schrodinger equation solver

# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include <cmath>


Eigenvalues::~Eigenvalues() {}

HarmonicEigenvalues::HarmonicEigenvalues(double omega, unsigned int n, int l): m_omega(omega), m_n(n), m_l(l) {}

double HarmonicEigenvalues::eigenvalue() const
{
    return Parameters::hbar*m_omega*(2*m_n+m_l+(3/2));
}

Eigenvalues* HarmonicEigenvalues::clone() const
{
    return new HarmonicEigenvalues(*this);
}

Schroddy::Schroddy(const InitialPot& pot, const Eigenvalues& eigenval): m_pot(pot.clone()), m_eigenval(eigenval.clone()) {}

Schroddy::~Schroddy()
{
    delete m_pot;
    delete m_eigenval;
}

double Schroddy::solveShroddyByRK(double x0, double x1, double psi0, double psiPrime0, unsigned long NSteps) const
{
    //implement runge-kutta here..
    const double h = (x1 - x0)/NSteps;
    const double factor = 2*Parameters::mn/(Parameters::hbar*Parameters::hbar);
    const double eigenvalue = m_eigenval -> eigenvalue();
    
    double runningX = x0, runningPsi = psi0, runningPsiPrime = psiPrime0;
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
        
        //advance running variables
        runningPsi += h/6*factor*(k1 + 2*k2 + 2*k3 + k4);
        runningPsiPrime =+ h/6*factor*(l1 + 2*l2 + 2*l3 + l4);
        runningX += h;
        
    }
    return runningPsi;
}
	
	

