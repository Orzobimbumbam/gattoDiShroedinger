# include "parameters.h"
# include "schroddy.h"
# include "kohn-sham.h"
# include "densities.h"
# include "initpot.h"
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
    double ratio1, ratio2, newPot;
    for (const auto& it : density.getTheoreticalDensity())
    {
        const double alpha = 1.;//ratio1;
        if ( inKSPot.getKSPot().at(it.first) < 0)
        {
        	ratio1 = 2 - it.second/density.getMCDensity().at(it.first);
        	//ratio1 = 2 - it.second/density.getSOGDensity().at(it.first);
        	if (ratio1 < 1 - Parameters::pregamma)
                ratio1 = 1 - Parameters::pregamma;
            // Prefactor
        	else if (ratio1 > 1 + Parameters::pregamma)
                ratio1 = 1 + Parameters::pregamma;	// condition

            newPot = alpha*ratio1*inKSPot.getKSPot().at(it.first);
        }
        else
        {
            ratio2 = it.second/density.getMCDensity().at(it.first);
            //ratio2 = it.second/density.getSOGDensity().at(it.first);
			if (ratio2 < 1 - Parameters::pregamma)
                ratio2 = 1 - Parameters::pregamma;		// Prefactor
            
			else if (ratio2 > 1 + Parameters::pregamma)
                ratio2 = 1 + Parameters::pregamma;	// condition
            
            newPot = alpha*ratio2*inKSPot.getKSPot().at(it.first);
        }
        
        m_KSOutPot[it.first] = newPot;
    }
	return;
}

KSPotential KohnShamInverse::getKSPot() const
{
	return m_KSOutPot;
}

/*bool KohnShamInverse::hasConverged (const KohnShamInverse& inKSPot) const
{
	double maxDiff = std::abs(m_KSOutPot.begin() -> second - inKSPot.getKSPot().begin() -> second);
	double xMax = m_KSOutPot.begin() -> first;
    for (const auto& it : m_KSOutPot)
	{
		if(std::abs(it.second - inKSPot.getKSPot().at(it.first)) > maxDiff) //access only, throw exception if key is not found
		{
			maxDiff = std::abs(it.second - inKSPot.getKSPot().at(it.first));
			xMax = it.first;
		}
	}

	const double epsilon = inKSPot.getKSPot().at(xMax)*50;
    m_distanceToConvergenge = maxDiff - epsilon;
	return maxDiff < epsilon; // convergence condition
}

double KohnShamInverse::distanceToConvergence() const
{
    return m_distanceToConvergenge;
}*/





