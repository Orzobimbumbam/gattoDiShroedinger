
#pragma once

//How many theoretical densities are you planning to have? If this is the only one, no need to define a virtual interface.
//Only use virtual and abstract classes when you want a class hierarchy (i.e. derived classes)
//Abstract classes MUST be derived and CANNOT be instantiated

class Densities
{
public:
    virtual double density() const = 0;
    virtual ~Densities();

    virtual Densities* clone() const = 0;

private:

};


class Theoreticaldensity: public Densities
{
public:
	Theoreticaldensity(double eigenfunc, double x);
    double density() const override;

    Densities* clone() const override;

private:
    Theoreticaldensity();
    const double m_x;
    double m_eigenfunc;
};


/*
class Theoreticaldensity
{
public:
    virtual double theodensity(double x) const = 0; //this line makes class ABSTRACT
    virtual ~Theoreticaldensity();

    virtual Theoreticaldensity* clone() const = 0; //this line makes class ABSTRACT

private:
    Theoreticaldensity(double normalPsi);
    const double m_normalPsi;
};
*/

//SUGGESTION: don't use virtual functions here. A normal class same as Schroddy will be perfectly fine.


