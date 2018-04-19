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
typedef std::pair<std::string, NuclearDensity> NuclearDensityOutputQuery;
typedef std::vector<double> DensityVector;

class NuclearDensity
{
public:
    void theoreticalDensity(const ElementEigenfunctions& psi, const OrderedLevelDegeneration& degen);
    void sogDensity (const std::vector<std::vector<double>>& QRparameters, double h);
    void mcDensity (std::ifstream& inStream);
    bool hasConverged () const;
    bool hasConvergedMC () const;
    
    Density getTheoreticalDensity() const;
    Density getSOGDensity() const;
    Density getMCDensity() const;
    
    double distanceToConvergence() const;
    double epsilon() const;
    double distanceToConvergenceMC() const;
    double epsilonMC() const;

    friend std::ostream& operator<<(std::ostream& wStream, const NuclearDensityOutputQuery& ouputQuery);
    
private:
    
    Density m_thDensity;
    Density m_sogDensity;
    Density m_mcDensity;
    mutable DensityVector m_theoDensity;
    mutable DensityVector m_MCDensity;
    mutable double m_distanceToConvergenge;
    mutable double m_epsilon;
    mutable double m_distanceToConvergengeMC;
    mutable double m_epsilonMC;

};

std::ostream& operator<<(std::ostream& wStream, const Density& density);
std::ostream& operator<<(std::ostream& wStream, const NuclearDensityOutputQuery& outputQuery);



