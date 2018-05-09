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

class KohnShamInverseWithJW : public KohnShamInverse
{
public:
	KohnShamInverseWithJW();
    KohnShamInverseWithJW(const InitialPot& iPot, double h);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) override;

};

class KohnShamInverseWithWP : public KohnShamInverse
{
public:
	KohnShamInverseWithWP();
    KohnShamInverseWithWP(const InitialPot& iPot, double h, const ElementEigenValues& eVal, const ElementEigenfunctions& inKSPsi);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) override;

};



