#pragma once
#include <memory>
#include <vector>

class Densities
{
public:
    virtual double density() const = 0;
    virtual ~Densities();

    virtual std::unique_ptr<Densities> clone() const = 0;

private:

};


class Theoreticaldensity: public Densities
{
public:
	Theoreticaldensity(char namefile);
    double density(double x0, double x1, std::vector<double> thDensArray, std::vector<double> xArray) const override;

    double getH() const {return m_h;}
    void setH(double H) {m_h = H;}

    virtual std::unique_ptr<Densities> clone() const override;

private:
    Theoreticaldensity();
    char m_namefile;
    mutable double m_h;

};




