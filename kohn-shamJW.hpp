#pragma once

#include <string>
#include <vector>
#include <map>
#include "idensity.h"
#include "kohn-sham.h"

class KohnShamInverseWithJW : public KohnShamInverse
{
public:
	KohnShamInverseWithJW();
	KohnShamInverseWithJW(const InitialPot& iPot, double h);
    
    void KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot) override;

};


