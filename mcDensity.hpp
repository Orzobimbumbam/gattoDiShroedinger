//
//  mcDensity.hpp
//  Codice
//
//  Created by Alberto Campi on 21/04/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#ifndef mcDensity_hpp
#define mcDensity_hpp

#include "idensity.h"

typedef std::vector<double> KeysArray;

class NuclearDensityWithMC : public NuclearDensity
{
public:
    NuclearDensityWithMC() : NuclearDensity() {}
    
    void benchmarkDensity (const std::vector<std::vector<double>>& mcDensity, double h = 0) override; //data are read and stored into a matrix in the client code (i.e. main), then passed in here
    Density getBenchmarkDensity() const override;
    
private:
    KeysArray getMatchingKeys() const;
};

#endif /* mcDensity_hpp */
