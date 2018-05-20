#include "mcDensity.hpp"

void NuclearDensityWithMC::benchmarkDensity(const std::vector<std::vector<double>>& mcDensity, double h)//h is optional here
{
    KeysArray arrayOfKeys = _getMatchingKeys();
    const long NRows = mcDensity.size();
    if (arrayOfKeys.size() == NRows)
    {
        m_benchmarkDensity.clear();
        for (unsigned long i = 0; i < NRows; ++i)
            m_benchmarkDensity.insert(std::make_pair(arrayOfKeys[i], mcDensity[i][1])); //align mcDensity to a map structure
    }
    else
        throw std::out_of_range("NuclearDensityWithMC::benchmarkDensity : mismatching key-value size."); //this behaviour can be changed once more info is available about mc data (should we interpolate?)
}

Density NuclearDensityWithMC::getBenchmarkDensity() const
{
    return m_benchmarkDensity;
}

// Load keys vector (radius values) from theoretical density map
KeysArray NuclearDensityWithMC::_getMatchingKeys() const
{
    Density thDens = getTheoreticalDensity();
    KeysArray arrayOfKeys;
    for (const auto& it : thDens)
        arrayOfKeys.push_back(it.first);
    
    return arrayOfKeys;
}
