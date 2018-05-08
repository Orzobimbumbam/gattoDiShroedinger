#pragma once
#include <string>
#include <vector>
#include <map>
#include "idensity.h"

class InitialPot;
typedef std::map<double, double> KSPotential;

class KohnShamInverse
{
public:
    KohnShamInverse();
    KohnShamInverse(const InitialPot& iPot, double h);
    virtual void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) = 0;
    
    virtual ~KohnShamInverse() = default;
    virtual const KSPotential getKSPot() const;
    
protected:
    KSPotential m_KSOutPot;
};



class KohnShamInverseWithLB : public KohnShamInverse
{
public:
	KohnShamInverseWithLB();
    KohnShamInverseWithLB(const InitialPot& iPot, double h);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) override;

};



