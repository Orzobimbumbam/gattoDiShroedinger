//class hierarchy interface goes in header
#pragma once
#include <memory>
#include "kohn-sham.h"
#include "kohn-shamJW.hpp"

/*====================================================================
 * Pure virtual base class for in/out potentials
 *===================================================================*/

class InitialPot
{
public:
    virtual double potential(double x) const = 0;
    virtual void setL(unsigned int l) {m_anglmomentum = l;};
    virtual double getL() const {return m_anglmomentum;};
    virtual ~InitialPot();
    
    virtual std::unique_ptr<InitialPot> clone() const = 0;
    
protected:
    const double m_m;
    unsigned int m_anglmomentum;
    InitialPot(double m, unsigned int l);

};

/*=====================================================================
 * Woods-Saxon potential class
 *===================================================================*/

class WSaxPot: public InitialPot //derived class
{

public:
    WSaxPot(double Rn, double a0, double m);
    WSaxPot(double Rn, double a0, double m, unsigned int l);
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
    double potential(double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

private:
    CoPot(); //same as before
    int m_Z;
};

/*=====================================================================
 * Total potential class
 *===================================================================*/

/*class TotPot: public InitialPot //derived class
{
public:
    TotPot(double m);
    TotPot(double m, unsigned int l);
    double potential(double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

private:
    TotPot(); //same as before
};

friend TotPot& operator+(const InitialPot& rhsPotential);*/

/*=====================================================================
 * Spin-Orbit potential class
 *===================================================================*/
/*
class SOPot: public InitialPot //derived class
{
public:

    SOPot(double m, double omega);
    double potential(double x) const override;

    //double getOmega() const;
    InitialPot* clone() const override;

private:
    SOPot(); //same as before
    const double m_m, m_omega;
};
*/

/*======================================================================
 * Kohn-Sham potential class
 *====================================================================*/

class PotOut: public InitialPot
{
public:
	PotOut(const KohnShamInverse& outpot, double m, unsigned int l);
    PotOut(const KohnShamInverse& outpot, double m); //l is assumed ground state
	double potential (double x) const override;
    
    std::unique_ptr<InitialPot> clone() const override;

protected:
    PotOut();
	double interpolatedPotential(double x) const; //protected only for testing purposes, fixture classes to be derived
    
private:
    KohnShamInverse m_outpot;
};

/*======================================================================
 * Kohn-Sham potential class from JW
 *====================================================================*/

class PotOutJW: public InitialPot
{
public:
	PotOutJW(const KohnShamInverseWithJW& outpot, double m, unsigned int l);
    PotOutJW(const KohnShamInverseWithJW& outpot, double m); //l is assumed ground state
	double potential (double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

protected:
    PotOutJW();
	double interpolatedPotential(double x) const; //protected only for testing purposes, fixture classes to be derived

private:
    KohnShamInverseWithJW m_outpot;
};

/*======================================================================
 * Test potential class
 *====================================================================*/

class TestPot: public InitialPot
{
public:
	TestPot(double m); //l is assumed ground state
	TestPot(double m, unsigned int l);
    double potential(double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

private:
    TestPot(); //same as before
};








