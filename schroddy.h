//double HBOEigen ( double omega, double l, double hbar, int n);
//double schrodinger (double r, double l, double vks, double u);
#pragma once
#include <vector>

class InitialPot;
class Schroddy;

/*==================================================================
 * Pure virtual class for eigenvalues
 *================================================================*/

class Eigenvalues
{
public:
    virtual double eigenvalue() const = 0;
    virtual ~Eigenvalues();

    virtual Eigenvalues* clone() const = 0;

private:

};

/*===================================================================
 * HO eigenvalues generator class
 *=================================================================*/

class HarmonicEigenvalues: public Eigenvalues
{
public:
    HarmonicEigenvalues(unsigned int n, int l);
    double eigenvalue() const override;
    //setEnergyLevel(unsigned int n);
    //setAngularMomentum(int l);

    Eigenvalues* clone() const override;

private:
    HarmonicEigenvalues();

    //if you want to be able to change quantum numbers within same object, implement set methods and remove const keyword
    const unsigned int m_n;
    const int m_l;
};


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

/*=====================================================================
 * Eigenvalues generator (shooting method) class
 *===================================================================*/

//singleton class for trial eigenvalues
class TrialEigenvalues
{
public:
    static double getEigenval1();
    static double getEigenval2();
    
    
private:
    static TrialEigenvalues* m_trialEigenvalusObj;
    const double m_eigenval1, m_eigenval2;
    TrialEigenvalues();
};


class GenericEigenvalues: public Eigenvalues
{
public:
    GenericEigenvalues (const Schroddy& sh, unsigned int nState, unsigned int lState);
    double eigenvalue() const override;
    
    Eigenvalues* clone() const override;
    
private:
    GenericEigenvalues();
    double shootingMethod(double E1, double E2, unsigned int parity) const;
    
    const Schroddy m_sh;
    unsigned int m_nState, m_lState;
    
};








