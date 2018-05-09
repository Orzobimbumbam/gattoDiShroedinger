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

private:
    ElementEigenValues m_eVal;
    ElementEigenfunctions m_inKSPsi;
    
    //DESIGN QUESTION: do you plan to create one object within the main scope and re-use it (same as we've been doing) or instantiate a new object after each loop? It is important to decide because the private fields are initialised in the constructor, which is called only once.
    
    //Solution 1: - Re-usable object -> create set methods for private fields, to be called after every loop
    //Solution 2: - Single-use object -> leave everything as it is and mark private fields as const
};



