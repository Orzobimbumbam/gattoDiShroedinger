# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include <cmath>
# include <vector>
# include <fstream>


/*=========================================================================
 * Schrodinger equation solver by Runge-Kutta method
 *=======================================================================*/

Schroddy::Schroddy(const InitialPot& pot) : m_pot(pot.clone()), m_h(0.001) {} //default precision
Schroddy::Schroddy(const InitialPot& pot, double H) : m_pot(pot.clone()), m_h(H) {}

Schroddy::Schroddy(const Schroddy& sourceSchroddy) : m_pot(sourceSchroddy.m_pot -> clone()), m_h(sourceSchroddy.m_h) {} //copy constructor

void Schroddy::swap(Schroddy& sourceSh)
{
    std::swap(m_h, sourceSh.m_h);
    std::swap(m_pot, sourceSh.m_pot);
}

Schroddy& Schroddy::operator=(const Schroddy& rhsSchroddy) //copy assignment
{
    Schroddy tempSh(rhsSchroddy); //cloning takes place here in the copy ctor
    swap(tempSh);
    
    return *this; //memory automatically released here by smart pointer as method goes out of scope
}

double Schroddy::solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, std::vector<double>& psiArray) const
{
    psiArray.clear();
    const unsigned long NSteps = static_cast<unsigned long>(std::abs(x1 - x0)/m_h);
    const double factor = 2*Parameters::mn/(Parameters::hbarc*Parameters::hbarc);
    const double eigenvalue = E;

    double runningX = x0, runningPsi = psi0, runningPsiPrime = psiPrime0;
    psiArray.push_back(psi0);

    for (unsigned long i = 0; i < NSteps ; ++i)
    {
        //lambda expression for decoupled angular term
        auto angularPart = [this](double x) -> const double {
            if (x == 0)
                return 0;
            return (Parameters::hbarc*Parameters::hbarc)*m_pot -> getL()*(m_pot -> getL() + 1)/(2*Parameters::mn*x*x);};
        
        //compute decoupled RK factors
        const double k1 = m_h*runningPsiPrime;
        const double l1 = m_h*factor*(m_pot -> potential(runningX) + angularPart(runningX) - eigenvalue)*runningPsi;
        const double k2 = m_h*(runningPsiPrime + 1./2.*l1);
        const double l2 = m_h*factor*(m_pot -> potential(runningX + m_h/2.) + angularPart(runningX + m_h/2.)
                                      - eigenvalue)*(runningPsi + 1./2.*k1);
        const double k3 = m_h*(runningPsiPrime + 1/2.*l2);
        const double l3 = m_h*factor*(m_pot -> potential(runningX + m_h/2.) + angularPart(runningX + m_h/2.)
                                      - eigenvalue)*(runningPsi + 1./2.*k2);
        const double k4 = m_h*(runningPsiPrime + 1*l3);
        const double l4 = m_h*factor*(m_pot -> potential(runningX + m_h) + angularPart(runningX + m_h) - eigenvalue)*(runningPsi + 1*k3);

        //advance running variables and store intermediate results
        runningPsi += 1./6.*(k1 + 2*k2 + 2*k3 + k4);
        runningPsiPrime += 1./6.*(l1 + 2*l2 + 2*l3 + l4);
        runningX += m_h;
        psiArray.push_back(runningPsi);
    }
    
    //work out normalizatiion constant
    double psiSquared = 0;
    for (const auto& it : psiArray)
        psiSquared += (it)*(it);

    const double scalar = psiSquared*m_h;
    
    //normalize radial eigenfunction element-wise -> [u(r)]*scalar^-1
    for (auto& it : psiArray)
        it = it/(sqrt(scalar));

    //return terminal value of normalized radial eigenfunction (on interval [x0, x1]) and job done!
    const double normalPsi = *(psiArray.end() - 1);
    return normalPsi;
}

// Calculate eigenfunctions and load a map defined in eigenfuntion.cpp
const Eigenfunction Schroddy::solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E) const
{
    Eigenfunction psi;
    std::vector<double> psiArray;
    
    double x = x0;
    solveSchroddyByRK(x0, x1, psi0, psiPrime0, E, psiArray);
    for (const auto& it : psiArray)
    {
        psi(x) = it;
        x += m_h;
    }
    
    return psi; //return an object with implemented move semantic
}




