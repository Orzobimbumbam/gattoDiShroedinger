#pragma once
#include <memory>
#include <string>
#include <vector>

class Theoreticaldensity
{
public:
	Theoreticaldensity();
    void density(const std::vector<double>& psi, std::vector<double>& thDensArray, unsigned int degen, double step) const;

private:
    //Theoreticaldensity();
};


class SOGdensity
{
public:
	SOGdensity();
	void sogDensity (const std::vector<std::vector<double> >& QRparameters, std::vector<double>& sogdensity, double h)const;

private:

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
