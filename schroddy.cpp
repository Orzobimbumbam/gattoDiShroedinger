//Class for Schrodinger equation solver

# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include <cmath>
# include <vector>
# include <fstream>
# include "NLSolverClass.hpp"
//# include "gnuplot_i.hpp"


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
        //no need to release pointer as it is done automatically
        m_pot = rhsSchroddy.m_pot -> clone();
        m_h = rhsSchroddy.m_h;
    }
    return *this;
}

double Schroddy::solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, std::vector<double>& psiArray) const
{
    psiArray.clear();
    const unsigned long NSteps = std::abs(x1 - x0)/m_h;
    const double factor = 2*Parameters::mn/(Parameters::hbarc*Parameters::hbarc);
    const double eigenvalue = E;
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
    double psiSquared = 0;
    for (const auto& it : psiArray)
        psiSquared += (it)*(it);

    const double scalar = psiSquared*m_h;
    
    //normalize radial eigenfunction element-wise -> [u(r)/r]*scalar^-1
    for (auto it = psiArray.begin(), itt = x.begin(); it != psiArray.end() && itt != x.end(); ++it, ++itt)
        *it = *it/(*itt*sqrt(scalar));

    //return terminal value of normalized radial eigenfunction (on interval [x0, x1]) and job done!
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
