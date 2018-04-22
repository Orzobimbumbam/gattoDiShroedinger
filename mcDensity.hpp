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

class NuclearDensityWithMC : public NuclearDensity
{
public:
    void benchmarkDensity (const std::vector<std::vector<double>>& mcDensity, double h = 0) override;
    Density getBenchmarkDensity() const override;
};

#endif /* mcDensity_hpp */
