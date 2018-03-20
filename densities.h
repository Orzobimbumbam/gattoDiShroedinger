#pragma once
#include <memory>
#include <string>
#include <vector>

class Theoreticaldensity
{
public:
	Theoreticaldensity();
    void density(const std::vector<double>& psi, std::vector<double>& thDensArray, unsigned int degen, double step) const;
    bool convergence (const std::vector<double>& empidensity, const std::vector<double>& thdensity) const;

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
