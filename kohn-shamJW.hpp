#pragma once

#include <string>
#include <vector>
#include <map>
#include "idensity.h"
#include "densities.h"

class NuclearDensity;
class InitialPot;
typedef std::map<double, double> JWKSPotential;
typedef std::vector<std::vector<double>> KSMatrix;

class KohnShamInverseWithJW
{
public:
	KohnShamInverseWithJW();
	KohnShamInverseWithJW(const InitialPot& iPot, double h);
    void KSinverseWithJW(const NuclearDensity& density, const KohnShamInverseWithJW& inKSPot);


    JWKSPotential getJWKSPot() const;

protected:
    JWKSPotential m_KSOutPot;
};


