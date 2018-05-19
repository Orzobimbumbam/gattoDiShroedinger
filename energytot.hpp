#ifndef energytot_hpp
#define energytot_hpp

#include "element.hpp"

class EnergyTOT
{
public:
    void energyTot (const ElementEigenfunctions& elEigf, const ElementEigenvalues& elEigV, double h) const;

private:

};

#endif /* energytot_hpp */


