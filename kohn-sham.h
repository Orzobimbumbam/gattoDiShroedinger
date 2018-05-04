#pragma once
#include <string>
#include <vector>
#include <map>
#include "idensity.h"
#include "densities.h"

class NuclearDensity;
class InitialPot;
typedef std::map<double, double> KSPotential;
typedef std::vector<std::vector<double>> KSMatrix;

class KohnShamInverse
{
public:
	KohnShamInverse();
    KohnShamInverse(const InitialPot& iPot, double h);
    void KSinverseWithLB(const NuclearDensity& density, const KohnShamInverse& inKSPot);
    //void KSinverseWithJW(const NuclearDensity& density, const KohnShamInverse& inKSPot);


    KSPotential getKSPot() const;

protected:
    KSPotential m_KSOutPot;
};



