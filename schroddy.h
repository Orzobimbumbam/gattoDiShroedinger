//double HBOEigen ( double omega, double l, double hbar, int n);
//double schrodinger (double r, double l, double vks, double u);
#pragma once

class InitialPot;

//hierarchy for eigenvalue objects

class Eigenvalues
{
public:
    virtual double eigenvalue(unsigned int n, int l) const = 0;
    virtual ~Eigenvalues();

private:
    
};



class HarmonicEigenvalues: public Eigenvalues
{
public:
    HarmonicEigenvalues(double omega);
    double eigenvalue(unsigned int n, int l) const override;
    
private:
    HarmonicEigenvalues();
    const double m_omega;
};


//Shroedinger class + solver


class Schroddy
{
public:
    Schroddy(); //if you have more parameters to pass in, change this
    double solveShroddyByRK(const InitialPot& pot, const Eigenvalues& eigenval, double x0, double x1, double step) const;
    
private:
    
};
