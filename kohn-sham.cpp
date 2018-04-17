# include "parameters.h"
# include "schroddy.h"
# include "kohn-sham.h"
# include "densities.h"
#include "initpot.h"
# include <vector>
# include <cmath>

/*===============================================================================
 * Kohm-Sham inverse equations
 *=============================================================================*/

KohnShamInverse::KohnShamInverse (): m_KSOutPot() {}
KohnShamInverse::KohnShamInverse (const InitialPot& iPot, double h)
{
    using namespace Parameters;
    const unsigned int NSteps = static_cast<unsigned int>(std::abs(x_fin - x_in)/h);
    
    double x = x_in;
    m_KSOutPot.insert(std::make_pair(x, iPot.potential(x)));
    for (unsigned int i = 0; i < NSteps; i++)
    {
        x += h;
        m_KSOutPot.insert(std::make_pair(x, iPot.potential(x)));
    }
}

void KohnShamInverse::KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot)
{
	//I'm assuming here that the two maps have the same keys, so I can use one iterator only
    //double ratio1, ratio2, newPot;
    for (const auto& it : density.getTheoreticalDensity())
    {
        double newPot;
        const double alpha = 1;
        if (inKSPot.getKSPot().at(it.first) < 0)
        {
            double ratio1 = 2 - it.second/density.getSOGDensity().at(it.first);
            if (alpha*ratio1 < 1 - Parameters::pregamma)
                ratio1 = 1 - Parameters::pregamma;
            
            else if (alpha*ratio1 > 1 + Parameters::pregamma)
                ratio1 = 1 + Parameters::pregamma;
            
            newPot = alpha*ratio1*inKSPot.getKSPot().at(it.first);
        }
        else
        {
            double ratio2 = it.second/density.getSOGDensity().at(it.first);
            /*
            if (alpha*ratio2 < 1 - Parameters::pregamma)
                ratio2 = 1 - Parameters::pregamma;
            
            else if (alpha*ratio2 > 1 + Parameters::pregamma)
                ratio2 = 1 + Parameters::pregamma;
            */
            newPot = ratio2*inKSPot.getKSPot().at(it.first);
        }
        
        m_KSOutPot[it.first] = newPot;
    }
	return;
}

KSPotential KohnShamInverse::getKSPot() const
{
	return m_KSOutPot;
}






