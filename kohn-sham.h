#pragma once
#include <string>
#include <vector>
#include <map>

class NuclearDensity;
class InitialPot;
typedef std::map<double, double> KSPotential;

class KohnShamInverse
{
public:
	KohnShamInverse();
    KohnShamInverse(const InitialPot& iPot, double h);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot);

    KSPotential getKSPot() const;

protected:
    KSPotential m_KSOutPot;
};



