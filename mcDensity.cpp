//
//  mcDensity.cpp
//  Codice
//
//  Created by Alberto Campi on 21/04/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#include "mcDensity.hpp"

void NuclearDensityWithMC::benchmarkDensity(const std::vector<std::vector<double>>& mcDensity, double h)//h is optional here
{
    KeysArray arrayOfKeys = getMatchingKeys();
    const long NRows = mcDensity.size();
    if (arrayOfKeys.size() == NRows)
    {
        m_benchmarkDensity.clear();
        for (unsigned long i = 0; i < NRows; ++i)
            m_benchmarkDensity.insert(std::make_pair(arrayOfKeys[i], mcDensity[i][1])); //align mcDensity into a map structure
    }
    else
        throw std::out_of_range("NuclearDensityWithMC::benchmarkDensity : mismatching key-value size."); //this behaviour can be changed once more info is available about mc data (should we interpolate?)
}

Density NuclearDensityWithMC::getBenchmarkDensity() const
{
    return m_benchmarkDensity;
}

KeysArray NuclearDensityWithMC::getMatchingKeys() const
{
    Density thDens = getTheoreticalDensity();
    KeysArray arrayOfKeys;
    for (const auto& it : thDens)
        arrayOfKeys.push_back(it.first);
    
    return arrayOfKeys;
}
