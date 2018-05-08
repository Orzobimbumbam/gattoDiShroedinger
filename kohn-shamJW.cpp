# include "parameters.h"
# include "kohn-shamJW.hpp"
# include "initpot.h"
# include <vector>
# include <cmath>

KohnShamInverseWithJW::KohnShamInverseWithJW() : KohnShamInverse() {}
KohnShamInverseWithJW::KohnShamInverseWithJW(const InitialPot& iPot, double h) : KohnShamInverse(iPot, h) {}


/*==============================================================================================
 * Kohm-Sham inverse equations Jensen-Wasserman method
 *============================================================================================*/

void KohnShamInverseWithJW::KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot)
{
	//I'm assuming here that the two maps have the same keys, so I can use one iterator only
    double ratio1, ratio2, newPot;
    for (const auto& it : density.getTheoreticalDensity())
    {
        const double alpha = 5;
        if ( inKSPot.getKSPot().at(it.first) < 0)
        {
        	ratio1 = (density.getBenchmarkDensity().at(it.first) - it.second)/density.getBenchmarkDensity().at(it.first);
            newPot = inKSPot.getKSPot().at(it.first) + alpha*ratio1;
        }
        else
        {
            ratio2 = (it.second - density.getBenchmarkDensity().at(it.first))/density.getBenchmarkDensity().at(it.first);
            newPot = inKSPot.getKSPot().at(it.first) + alpha*ratio2;
        }

        m_KSOutPot[it.first] = newPot;
    }
	return;
}




