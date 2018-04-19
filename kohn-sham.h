#pragma once
#include <string>
#include <vector>
#include <map>

class NuclearDensity;
class InitialPot;
typedef std::map<double, double> KSPotential;
typedef std::vector<std::vector<double>> KSMatrix;

class KohnShamInverse
{
public:
	KohnShamInverse();
    KohnShamInverse(const InitialPot& iPot, double h);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot);
    void KSinverseMC(const NuclearDensity& density, const KohnShamInverse& inKSPot);


    KSPotential getKSPot() const;

protected:
    KSPotential m_KSOutPot;
    KSMatrix m_theoDensity;
    KSMatrix m_mcDensity;
    KSMatrix m_inKSPot;
};



