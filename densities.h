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

class NuclearDensity
{
public:
    void theoreticalDensity(const ElementEigenfunctions& psi, const OrderedLevelDegeneration& degen);
    void sogDensity (const std::vector<std::vector<double>>& QRparameters, double h);
    bool hasConverged () const;
    
    Density getTheoreticalDensity() const;
    Density getSOGDensity() const;

    //friend std::ostream& operator<<(std::ostream& wStream, const NuclearDensity& thDensity);
    
private:
    
    Density m_thDensity;
    Density m_sogDensity;
};

std::ostream& operator<<(std::ostream& wStream, const Density& density);

