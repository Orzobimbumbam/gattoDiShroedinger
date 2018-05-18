#ifndef energytot_hpp
#define energytot_hpp

#include "element.hpp"

class EnergyTOT
{
public:
    void energyTot (const ElementEigenfunctions& elEigf, ElementEigenValues& eigenValMatrix) override;

private:

};

#endif /* energytot_hpp */


