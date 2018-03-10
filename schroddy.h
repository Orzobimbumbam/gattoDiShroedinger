//double HBOEigen ( double omega, double l, double hbar, int n);
//double schrodinger (double r, double l, double vks, double u);
#pragma once
#include <vector>

class InitialPot;
class Schroddy;

/*=======================================================================
 * Shroedinger class + solver
 *=====================================================================*/

class Schroddy
{
public:
    Schroddy(const InitialPot& pot/*, double H*/);
    Schroddy(const InitialPot& pot, double H);
    Schroddy(const Schroddy& sourceSh);
    Schroddy& operator=(const Schroddy& rhsSchroddy);
    
    double getH() const {return m_h;}
    void setH(double H) {m_h = H;}
    InitialPot* getInitialPotPtr() const {return m_pot;};

    double solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, std::vector<double>& psiArray ) const; //psiPrime0 is boundary condition on first derivative of eigenfunction

    ~Schroddy();

private:
    Schroddy();
    InitialPot* m_pot;
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









