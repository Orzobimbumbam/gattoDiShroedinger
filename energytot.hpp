#ifndef energytot_hpp
#define energytot_hpp

#include "element.hpp"

class EnergyTOT
{
public:
    void energyTot (const ElementEigenfunctions& elEigf, const ElementEigenValues& elEigV, double h) override;

private:

};

#endif /* energytot_hpp */


