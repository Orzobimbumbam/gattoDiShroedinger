
#pragma once

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
	Theoreticaldensity(double eigenfunc, int degen, double x);
    double density() const override;

    Densities* clone() const override;

private:
    Theoreticaldensity();
    const double m_x;
    double m_eigenfunc;
    double m_degen;
};




