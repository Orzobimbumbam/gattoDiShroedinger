#ifndef energytot_hpp
#define energytot_hpp

#include "element.hpp"

class EnergyTOT
{
public:
    
    double energyTot(const ElementEigenfunctions& elEigf, const ElementEigenvalues& elEigV) const;

private:
    double _laplacian(std::map<double, double>::const_iterator itt, double h) const;

};

#endif /* energytot_hpp */


