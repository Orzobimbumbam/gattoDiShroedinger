//class hierarchy interface goes in header
#pragma once

class InitialPot //pure virtual base class
{
public:
    virtual double potential(double x) const = 0;
    virtual ~InitialPot();

    virtual InitialPot* clone() const = 0; //virtual constructor

private:

};

/*
class WSaxPot: public InitialPot //derived class
{
    //=====================================================================
    // Woods-Saxon potential
    //=====================================================================
public:
    WSaxPot(double V0, double Rn, double a0);
    double potential(double x) const override;

    InitialPot* clone() const override;

private:
    WSaxPot(); //make defaul ctor private as we must initialize parameters
    const double m_V0, m_Rn, m_a0;

};
*/

class HOPot: public InitialPot //derived class
{
public:
    //=====================================================================
    // HO potential
    //=====================================================================
    HOPot(double m, double omega, double hbar, double angular);
    double potential(double x) const override;

    //double getOmega() const;
    InitialPot* clone() const override;

private:
    HOPot(); //same as before
    const double m_m, m_omega, m_hbar, m_angular;
};






