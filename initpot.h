//class hierarchy interface goes in header

class InitialPot //pure virtual base class
{
public:
    virtual double potential(double x) const = 0;
    virtual ~InitialPot();
    
private:
    
};


class WSaxPot: public InitialPot //derived class
{
    //=====================================================================
    // Woods-Saxon potential
    //=====================================================================
public:
    WSaxPot(double V0, double Rn, double a0);
    double potential(double x) const override;
    
private:
    WSaxPot(); //make defaul ctor private as we must initialize parameters
    const double m_V0, m_Rn, m_a0;
    
};


class HBOPot: public InitialPot //derived class
{
public:
    //=====================================================================
    // HBO potenzial
    //=====================================================================
    HBOPot(double m, double omega, double rhbo);
    double potential(double x) const override;
    
private:
    HBOPot(); //same as before
    const double m_m, m_omega;
};


