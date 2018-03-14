#pragma once
#include <memory>
#include <string>
#include <vector>

class Theoreticaldensity
{
public:
	Theoreticaldensity();
    double density(double psi, int degen, double x) const;

private:
    //Theoreticaldensity();
};


/*class Densities
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
	Theoreticaldensity(char filename);
    double density(double x0, double x1, std::vector<double> thDensArray, std::vector<double> xArray) const override;

    double getH() const {return m_h;}
    void setH(double H) {m_h = H;}

    virtual std::unique_ptr<Densities> clone() const override;

private:
    Theoreticaldensity();
    char m_filename;
    mutable double m_h;

};*/

