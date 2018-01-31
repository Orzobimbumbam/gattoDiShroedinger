//double HBOEigen ( double omega, double l, double hbar, int n);
//double schrodinger (double r, double l, double vks, double u);
#pragma once

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
    HarmonicEigenvalues(double omega, unsigned int n, int l);
    double eigenvalue() const override;
    //setEnergyLevel(unsigned int n);
    //setAngularMomentum(int l);

    Eigenvalues* clone() const override;

private:
    HarmonicEigenvalues();
    const double m_omega;

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
    Schroddy(const InitialPot& pot/*, const Eigenvalues& eigenval*/);
    Schroddy(const Schroddy& sourceSh);
    Schroddy& operator=(const Schroddy& rhsSchroddy);
    
    double solveShroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, unsigned long NSteps) const; //psiPrime0 is boundary condition on first derivative of eigenfunction

    ~Schroddy();

private:
    InitialPot* m_pot;
    //const Eigenvalues* const m_eigenval;
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
	//double m_x0, m_x1, m_psi0, m_psiPrime0, m_NSteps;
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
    GenericEigenvalues (const InitialPot& initPot);
    double eigenvalue() const override;
    
    Eigenvalues* clone() const override;
    
private:
    GenericEigenvalues();
    InitialPot* m_initPot;
};








