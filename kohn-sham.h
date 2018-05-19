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
    KohnShamInverse(const InitialPot& iPot, double h);
    virtual void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) = 0;

    virtual ~KohnShamInverse() = default;
    virtual const KSPotential getKSPot() const;
    
protected:
    KohnShamInverse() = default;
    KSPotential m_KSOutPot;
};

class KohnShamInverseWithLB : public KohnShamInverse
{
public:
    KohnShamInverseWithLB(const InitialPot& iPot, double h);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) override;
    
private:
    KohnShamInverseWithLB();

};

class KohnShamInverseWithJW : public KohnShamInverse
{
public:
    KohnShamInverseWithJW(const InitialPot& iPot, double h);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) override;
    
private:
    KohnShamInverseWithJW();

};

class KohnShamInverseWithWP : public KohnShamInverse
{
public:
    KohnShamInverseWithWP(const InitialPot& iPot, double h, const ElementEigenvalues& eVal, const ElementEigenfunctions& inKSPsi);
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) override;
    
    void setElementEigenvalues(const ElementEigenvalues& eVal);
    void setElementEigenfunctions(const ElementEigenfunctions& inKSPsi);

private:
    KohnShamInverseWithWP();
    ElementEigenvalues m_eVal;
    ElementEigenfunctions m_inKSPsi;
    
    //[Orzobimbumbam] :
    //DESIGN QUESTION: do you plan to create one object within the main scope and re-use it (same as we've been doing) or instantiate a new object after each loop? It is important to decide because the private fields are initialised in the constructor, which is called only once.
    
    //Solution 1: - Re-usable object -> create set methods for private fields, to be called after every loop
    //Solution 2: - Single-use object -> leave everything as it is and mark private fields as const
};



