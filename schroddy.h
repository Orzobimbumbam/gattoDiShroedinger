//double HBOEigen ( double omega, double l, double hbar, int n);
//double schrodinger (double r, double l, double vks, double u);
#pragma once

class InitialPot;

//hierarchy for eigenvalue objects

class Eigenvalues
{
public:
    virtual double eigenvalue() const = 0;
    virtual ~Eigenvalues();
    
    virtual Eigenvalues* clone() const = 0;

private:
    
};


class HarmonicEigenvalues: public Eigenvalues
{
public:
    HarmonicEigenvalues(double omega, unsigned int n, int l);
    double eigenvalue() const override;
    
    Eigenvalues* clone() const override;
    
private:
    HarmonicEigenvalues();
    const double m_omega;
    unsigned int m_n;
    int m_l;
};


//Shroedinger class + solver


class Schroddy
{
public:
    Schroddy(const InitialPot& pot, const Eigenvalues& eigenval); //if you have more parameters to pass in, change this
    double solveShroddyByRK(double x0, double x1, double psi0, double psiPrime0, unsigned long NSteps) const; //psiPrime0 is boundary condition on first derivative of eigenfunction
    
    ~Schroddy();
    
private:
    const InitialPot* const m_pot;
    const Eigenvalues* const m_eigenval;
    
};
