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
    /*bool hasConverged (const KohnShamInverse& inKSPot) const;

    double distanceToConvergence() const;*/

    KSPotential getKSPot() const;

protected:
    KSPotential m_KSOutPot;
    //mutable double m_distanceToConvergenge;
};



