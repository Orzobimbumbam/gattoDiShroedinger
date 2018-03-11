#pragma once
#include <vector>
#include <memory>
#include "initpot.h"


/*=======================================================================
 * Shroedinger class + solver
 *=====================================================================*/

class Schroddy
{
public:
    Schroddy(const InitialPot& pot);
    Schroddy(const InitialPot& pot, double H);
    Schroddy(const Schroddy& sourceSh);
    Schroddy& operator=(const Schroddy& rhsSchroddy);
    
    double getH() const {return m_h;}
    void setH(double H) {m_h = H;}
    std::unique_ptr<InitialPot>& getInitialPotPtr() {return m_pot;};

    double solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, std::vector<double>& psiArray) const; //psiPrime0 is boundary condition on first derivative of eigenfunction


private:
    Schroddy();
    std::unique_ptr<InitialPot> m_pot;
    mutable double m_h;
};

/*=======================================================================
 * Wrapper for resolve arguments inconsistency between schroddy and NLSolver classes
 *=====================================================================*/

class schroddywrapper
{
public:
	schroddywrapper (const Schroddy& sh);
	double eigenfunction (double E) const;

private:
	Schroddy m_sh;
};









