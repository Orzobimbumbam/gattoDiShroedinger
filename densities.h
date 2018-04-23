#pragma once
#include <memory>
#include <string>
#include <vector>
#include <map>

#include "IOUtils.hpp"
//#include "element.hpp"
#include "idensity.h"


typedef std::vector<double> DensityVector;

class NuclearDensityWithSOG : public NuclearDensity
{
public:
    NuclearDensityWithSOG() : NuclearDensity() {};
    //void theoreticalDensity(const ElementEigenfunctions& psi, const OrderedLevelDegeneration& degen);
    void benchmarkDensity (const std::vector<std::vector<double>>& QRparameters, double h = 0) override;
    //void mcDensity (std::ifstream& inStream);
    //bool hasConverged () const override;
    //bool hasConvergedMC () const;
    
    //Density getTheoreticalDensity() const;
    Density getBenchmarkDensity() const override;
    //Density getMCDensity() const;
    
    //double distanceToConvergence() const;
    //double epsilon() const;
    //double distanceToConvergenceMC() const;
    //double epsilonMC() const;

    friend std::ostream& operator<<(std::ostream& wStream, const NuclearDensityOutputQuery& ouputQuery);
    
private:
    
    //Density m_thDensity;
    //Density m_sogDensity;
    //Density m_mcDensity;
    //mutable DensityVector m_theoDensity;
    //mutable DensityVector m_MCDensity;
    //mutable double m_distanceToConvergenge;
    //mutable double m_epsilon;
    //mutable double m_distanceToConvergengeMC;
    //mutable double m_epsilonMC;

};

std::ostream& operator<<(std::ostream& wStream, const Density& density);
std::ostream& operator<<(std::ostream& wStream, const NuclearDensityOutputQuery& outputQuery);



