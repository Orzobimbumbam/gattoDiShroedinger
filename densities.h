#pragma once
#include <memory>
#include <string>
#include <vector>
#include <map>

#include "IOUtils.hpp"
#include "idensity.h"


typedef std::vector<double> DensityVector;

class NuclearDensityWithSOG : public NuclearDensity
{
public:
    NuclearDensityWithSOG() : NuclearDensity() {};
    void benchmarkDensity (const std::vector<std::vector<double>>& QRparameters, double h = 0) override;
   
    Density getBenchmarkDensity() const override;
    friend std::ostream& operator<<(std::ostream& wStream, const NuclearDensityOutputQuery& ouputQuery);
    
private:

};

std::ostream& operator<<(std::ostream& wStream, const Density& density);
std::ostream& operator<<(std::ostream& wStream, const NuclearDensityOutputQuery& outputQuery);



