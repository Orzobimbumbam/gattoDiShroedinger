#pragma once
#include <memory>
#include <string>
#include <vector>
#include <map>

class Eigenfunction;
typedef std::map<double, double> ThDensity;

class Theoreticaldensity
{
public:
    void density(const Eigenfunction& psi, unsigned int degen, double step);
    bool hasConverged (const std::map<double, double>& empidensity) const; //empidensity should be wrapped up in a singleton object (static)

private:
    ThDensity m_thDensity;
};


class SOGdensity
{
public:
	SOGdensity();
	void sogDensity (const std::vector<std::vector<double> >& QRparameters, std::vector<double>& sogdensity, double h)const;

private:

};
