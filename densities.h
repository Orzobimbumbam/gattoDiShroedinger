#pragma once
#include <memory>
#include <string>
#include <vector>
#include <map>

#include "IOUtils.hpp"

class Eigenfunction;
//class
typedef std::map<double, double> ThDensity;

class Theoreticaldensity
{
public:
    void density(const Eigenfunction& psi, unsigned int degen);
    bool hasConverged (const std::map<double, double>& empidensity) const; //empidensity should be wrapped up in a singleton object (static)

    std::ostream& operator<<(std::ostream& wStream) const;
    
private:
    friend std::ostream& writeMap(const ThDensity&, std::ostream& wStream);
    
    ThDensity m_thDensity;
};


class SOGdensity
{
public:
	SOGdensity();
	void sogDensity (const std::vector<std::vector<double> >& QRparameters, std::vector<double>& sogdensity, double h)const;

private:

};
