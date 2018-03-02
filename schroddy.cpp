//Class for Schrodinger equation solver

# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include <cmath>
# include <vector>
# include <fstream>
# include "NLSolverClass.hpp"
//# include "gnuplot_i.hpp"


Eigenvalues::~Eigenvalues() {}

/*=======================================================================
 * HO Eigenvalues generator
 *=====================================================================*/

HarmonicEigenvalues::HarmonicEigenvalues(unsigned int n, int l): m_n(n), m_l(l) {}
double HarmonicEigenvalues::eigenvalue() const
{
    return Parameters::hbar_omega*(2*m_n+m_l+(3/2));
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
    
    std::vector <double> psiArray;
    unsigned int nodes = 0;
    while (true)
    {
        
        double s = m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(m_lState),
                                          psiPrime0(m_lState), E2 , psiArray);
        
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
    
    double midE = (E1 + E2)/2.;
    
    while (std::abs(E2 - E1) > error)
    {
        double s = m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(m_lState),
                                          psiPrime0(m_lState), midE , psiArray);
        
        if (parityFlag*s > 0)
            E1 = midE;
        else
            E2 = midE;
        
        midE = (E1 + E2)/2.;
    }
    
    return midE;
}

double GenericEigenvalues::eigenvalue() const //this must return a double..
{
    double E1 = TrialEigenvalues::getEigenval1();
    
    for (unsigned int i = 1; i <= m_nState; ++i)
    {
        E1 = shootingMethod(E1, TrialEigenvalues::getEigenval2(), i);
    }
    
    return E1;

}



Eigenvalues* GenericEigenvalues::clone() const
{
	return new GenericEigenvalues(*this);
}

/*=========================================================================
 * Schrodinger equation solver by Runge-Kutta method
 *=======================================================================*/

Schroddy::Schroddy(const InitialPot& pot): m_pot(pot.clone()), m_h(0.001) {}//default precision
Schroddy::Schroddy(const InitialPot& pot, double H): m_pot(pot.clone()), m_h(H) {}

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

double Schroddy::solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, std::vector<double>& psiArray) const
{
    psiArray.clear();
    const unsigned long NSteps = std::abs(x1 - x0)/m_h;
    const double factor = 2*Parameters::mn/(Parameters::hbarc*Parameters::hbarc);
    const double eigenvalue = E;//m_eigenval -> eigenvalue();
    //std::ofstream file("RKout.txt");

    double runningX = x0, runningPsi = psi0, runningPsiPrime = psiPrime0;
    psiArray.push_back(psi0);
    std::vector<double> x;
    x.push_back(runningX);

    for (unsigned long i = 0; i < NSteps ; ++i)
    {

        //compute decoupled RK factors
        const double k1 = m_h*runningPsiPrime;
        const double l1 = m_h*factor*(m_pot -> potential(runningX) - eigenvalue)*runningPsi;
        const double k2 = m_h*(runningPsiPrime + 1./2.*l1);
        const double l2 = m_h*factor*(m_pot -> potential(runningX + m_h/2.) - eigenvalue)*(runningPsi + 1./2.*k1);
        const double k3 = m_h*(runningPsiPrime + 1/2.*l2);
        const double l3 = m_h*factor*(m_pot -> potential(runningX + m_h/2.) - eigenvalue)*(runningPsi + 1./2.*k2);
        const double k4 = m_h*(runningPsiPrime + 1*l3);
        const double l4 = m_h*factor*(m_pot -> potential(runningX + m_h) - eigenvalue)*(runningPsi + 1*k3);

        //advance running variables and store intermediate results
        runningPsi += 1./6.*(k1 + 2*k2 + 2*k3 + k4);
        runningPsiPrime += 1./6.*(l1 + 2*l2 + 2*l3 + l4);
        runningX += m_h;
        psiArray.push_back(runningPsi);
        x.push_back(runningX);
    }
    /*
    std::vector <double>::iterator walk = psiArray.begin();
    while (walk != psiArray.end())
=======
        psiArray.push_back(runningPsi/runningX);
        psiRadius.push_back(runningX);

    }

    std::vector<double>::iterator walk1 = psiRadius.begin();
    std::vector<double>::iterator walk2 = psiArray.begin();
    while (walk1 != psiRadius.end() && walk2 != psiArray.end())

    {
    	file << *walk1 << "\t"<< *walk2 << std::endl;
    	//std::cout << *walk << std::endl;
    	walk1++, walk2++;
    }


    file.close();


    const std::string f = "RKout.txt";
    Gnuplot gp("lines");
    gp.set_title("Plotfile\\nNew Line");
    //gp.plotfile_xy(&f,1,2,'Psi');
    gp.plotfile_xy("RKout.txt",1,2,"Psi");
    gp.unset_title();*/



    //work out normalizatiion constant
    double psiSquared = 0;// normalPsi=0;
    for (auto it = psiArray.begin(); it != psiArray.end(); ++it)
        psiSquared += (*it)*(*it); //derefence iterator at array's element and square (brackets are needed!)

    const double scalar = psiSquared*m_h;
    for (auto it = psiArray.begin(), itt = x.begin(); it != psiArray.end() && itt != x.end(); ++it, ++itt)
        *it = *it/(*itt*sqrt(scalar));

    //return normalized eigenfunction (on interval [x0, x1]) and job done!
    const double normalPsi = runningPsi/sqrt(scalar);
    return normalPsi;
}

/*========================================================================
 * Wrapper for resolve arguments inconsistency between schroddy and NLSolver classes
 *======================================================================*/

//schroddywrapper::schroddywrapper (const Schroddy& sh, double x0, double x1, double psi0, double psiPrime0, unsigned long NSteps): m_sh(sh), m_x0(x0), m_x1(x1), m_psi0(psi0), m_psiPrime0(psiPrime0), m_NSteps(NSteps) {}
schroddywrapper::schroddywrapper (const Schroddy& sh): m_sh(sh) {}
/*
double schroddywrapper::eigenfunction(double E) const
{
	std::vector <double> psiArray;
    return m_sh.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0,
    		Parameters::psiPrime0, E, psiArray );
}*/


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
