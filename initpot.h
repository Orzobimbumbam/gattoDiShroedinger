//class hierarchy interface goes in header
#pragma once
#include <memory>
#include "kohn-sham.h"

/*====================================================================
 * Abstract base class for in/out potentials
 *===================================================================*/

class InitialPot
{
public:
    virtual double potential(double x) const = 0;
    
    virtual void setL(unsigned int l) {m_anglmomentum = l;};
    virtual unsigned int getL() const {return m_anglmomentum;};
    
    virtual void setJ(int j) {m_spin = j;};
    virtual int getJ() const {return m_spin;};
    
    virtual ~InitialPot();
    
    virtual std::unique_ptr<InitialPot> clone() const = 0;
    
protected:
    const double m_m;
    unsigned int m_anglmomentum;
    double m_spin;
    InitialPot(double m, unsigned int l, double j);

};

/*=====================================================================
 * Woods-Saxon potential class
 *===================================================================*/

class WSaxPot: public InitialPot //derived class
{

public:
    WSaxPot(double Rn, double a0, double m);
    WSaxPot(double Rn, double a0, double m, unsigned int l);
    WSaxPot(double Rn, double a0, double m, unsigned int l, double j);
    double potential(double x) const override;
    
    std::unique_ptr<InitialPot> clone() const override;

private:
    WSaxPot(); //make defaul ctor private as we must initialize parameters
    const double m_Rn, m_a0;

};


/*=====================================================================
 * HO potential class
 *===================================================================*/

class HOPot: public InitialPot //derived class
{
public:
    HOPot(double m); //l is assumed ground state
    HOPot(double m, unsigned int l);
    HOPot(double m, unsigned int l, double j);
    double potential(double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

private:
    HOPot(); //same as before
    //const double m_m;
};

/*=====================================================================
 * Coulomb potential class
 *===================================================================*/

class CoPot: public InitialPot //derived class
{
public:
    CoPot(int Z, double m);
    CoPot(int Z, double m, unsigned int l);
    CoPot(int Z, double m, unsigned int l, double j);
    double potential(double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

private:
    CoPot(); //same as before
    int m_Z;
};

/*======================================================================
 * Kohn-Sham potential class
 *====================================================================*/

class PotOut: public InitialPot
{
public:
    PotOut(const KohnShamInverse& outpot, double m, unsigned int l, double j);
	PotOut(const KohnShamInverse& outpot, double m, unsigned int l);
    PotOut(const KohnShamInverse& outpot, double m); //l is assumed ground state; j is zero
    
    PotOut(const KSPotential& inPot, double m);
	double potential (double x) const override;
    
    std::unique_ptr<InitialPot> clone() const override;

protected:
    PotOut();
	double interpolatedPotential(double x) const; //protected only for testing purposes, fixture classes to be derived
    
private:
    //KohnShamInverse m_outpot;
    KSPotential m_outpot;
};

/*======================================================================
 * HO Test potential class
 *====================================================================*/

class TestPot: public InitialPot
{
public:
	TestPot(double m); //l is assumed ground state
	TestPot(double m, unsigned int l);
    TestPot(double m, unsigned int l, double j);
    double potential(double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

private:
    TestPot(); //same as before
};








