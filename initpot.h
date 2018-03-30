//class hierarchy interface goes in header
#pragma once
#include <memory>
#include "kohn-sham.h"

/*====================================================================
 * Pure virtual base class for in/out potentials
 *===================================================================*/

class InitialPot
{
public:
    virtual double potential(double x) const = 0;
    virtual void setL(unsigned int l) {m_anglmomentum = l;};
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
    WSaxPot(double V0, double Rn, double a0, double m);
    WSaxPot(double V0, double Rn, double a0, double m, unsigned int l);
    double potential(double x) const override;
    //void setL(unsigned int l) override;
    
    std::unique_ptr<InitialPot> clone() const override;

private:
    WSaxPot(); //make defaul ctor private as we must initialize parameters
    const double m_V0, m_Rn, m_a0;

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
	double potential (double x) const override;

    std::unique_ptr<InitialPot> clone() const override;

protected:
    PotOut();
	double interpolatedPotential(double x) const; //protected only for testing purposes
    
private:
    KohnShamInverse m_outpot;

};








