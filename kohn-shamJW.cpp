# include "parameters.h"
# include "kohn-shamJW.hpp"
# include "initpot.h"
# include <vector>
# include <cmath>

/*==========================================================================================
 * Kohm-Sham inverse equations
 *========================================================================================*/

KohnShamInverseWithJW::KohnShamInverseWithJW (): m_KSOutPot() {}
KohnShamInverseWithJW::KohnShamInverseWithJW (const InitialPot& iPot, double h)
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


/*==============================================================================================
 * Kohm-Sham inverse equations Jensen-Wasserman method
 *============================================================================================*/

void KohnShamInverseWithJW::KSinverseWithJW(const NuclearDensity& density, const KohnShamInverseWithJW& inKSPot)
{
	//I'm assuming here that the two maps have the same keys, so I can use one iterator only
    double ratio1, ratio2, newPot;
    for (const auto& it : density.getTheoreticalDensity())
    {
        const double alpha = 5;
        if ( inKSPot.getJWKSPot().at(it.first) < 0)
        {
        	ratio1 = (density.getBenchmarkDensity().at(it.first) - it.second)/density.getBenchmarkDensity().at(it.first);
            newPot = inKSPot.getJWKSPot().at(it.first) + alpha*ratio1;
        }
        else
        {
            ratio2 = (it.second - density.getBenchmarkDensity().at(it.first))/density.getBenchmarkDensity().at(it.first);
            newPot = inKSPot.getJWKSPot().at(it.first) + alpha*ratio2;
        }

        m_KSOutPot[it.first] = newPot;
    }
	return;
}

KSPotential KohnShamInverseWithJW::getJWKSPot() const
{
    return m_KSOutPot;
}




