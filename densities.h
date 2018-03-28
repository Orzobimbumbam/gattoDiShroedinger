#pragma once
#include <memory>
#include <string>
#include <vector>
#include <map>

#include "IOUtils.hpp"
#include "element.hpp"

class Eigenfunction;
class Element;
typedef std::map<double, double> Density;

class Theoreticaldensity
{
public:
    //void density(const Element& element);
    void density(const ElementEigenfunctions& psi, const OrderedLevelDegeneration& degen);
    bool hasConverged (const std::map<double, double>& empidensity) const; 

    friend std::ostream& operator<<(std::ostream& wStream, const Theoreticaldensity& thDensity);
    
private:
    friend std::ostream& writeMap(const Density&, std::ostream& wStream, bool header);
    
    Density m_thDensity;
};

std::ostream& operator<<(std::ostream& wStream, const Theoreticaldensity& thDensity);

class SOGdensity
{
public:
	SOGdensity();
	void sogDensity (const std::vector<std::vector<double> >& QRparameters, std::vector<double>& sogdensity, double h)const;

private:

};
