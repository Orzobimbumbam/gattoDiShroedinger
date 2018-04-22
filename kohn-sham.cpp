# include "parameters.h"
# include "kohn-sham.h"
//# include "densities.h"
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

/*===============================================================================
 * Kohm-Sham inverse equations method for SOG Densities
 *=============================================================================*/

void KohnShamInverse::KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot)
{
	//I'm assuming here that the two maps have the same keys, so I can use one iterator only
    double ratio1, ratio2, newPot;
    for (const auto& it : density.getTheoreticalDensity())
    {
        const double alpha = 1.;//ratio1;
        if ( inKSPot.getKSPot().at(it.first) < 0)
        {
        	ratio1 = 2 - it.second/density.getBenchmarkDensity().at(it.first);
        	if (ratio1 < 1 - Parameters::pregamma)
                ratio1 = 1 - Parameters::pregamma;
            // Prefactor
        	else if (ratio1 > 1 + Parameters::pregamma)
                ratio1 = 1 + Parameters::pregamma;	// condition

            newPot = alpha*ratio1*inKSPot.getKSPot().at(it.first);
        }
        else
        {
            ratio2 = it.second/density.getBenchmarkDensity().at(it.first);
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


/*===============================================================================
 * Kohm-Sham inverse equations method for Monte-Carlo Densities
 *=============================================================================*/
/*
void KohnShamInverse::KSinverseMC(const NuclearDensity& density, const KohnShamInverse& inKSPot)
{
	// Load matrices from maps to avoid keys incompatibility
	KSMatrix temp(density.getTheoreticalDensity().size(), std::vector<double>(2));
	m_theoDensity = temp;
	m_mcDensity = temp;
	m_inKSPot = temp;

	int i = 0;
	for (const auto& it : density.getTheoreticalDensity())
	{
		m_theoDensity[i][0] = m_inKSPot[i][0] = it.first;
		m_theoDensity[i][1] = it.second;
		m_inKSPot[i][1] = inKSPot.getKSPot().at(it.first);
		++i;
	}

	int j = 0;
	for (const auto& it : density.getMCDensity())
	{
		m_mcDensity[j][0] = it.first;
		m_mcDensity[j][1] = density.getMCDensity().at(it.first);
		++j;
	}

	// New potential by Kohm-Sham inverse equations
	m_KSOutPot.clear();
	double ratio1, ratio2, newPot;
    for (int i = 0; i < m_theoDensity.size(); ++i)
    {
    	const double alpha = 1.;//ratio1;
    	if ( m_inKSPot[i][1] < 0)
    	{
    		ratio1 = 2 - ((m_theoDensity[i][1])/(m_mcDensity[i][1]));
    		if (ratio1 < 1 - Parameters::pregamma)
    			ratio1 = 1 - Parameters::pregamma;
    		// Prefactor
    		else if (ratio1 > 1 + Parameters::pregamma)
    			ratio1 = 1 + Parameters::pregamma;	// condition

    		newPot = alpha*ratio1*m_inKSPot[i][1];
    	}
    	else
    	{
    		ratio2 = (m_theoDensity[i][1])/(m_mcDensity[i][1]);
    		if (ratio2 < 1 - Parameters::pregamma)
    			ratio2 = 1 - Parameters::pregamma;		// Prefactor

    		else if (ratio2 > 1 + Parameters::pregamma)
    			ratio2 = 1 + Parameters::pregamma;	// condition

    		newPot = alpha*ratio2*m_inKSPot[i][1];
    	}

    	std::pair<double, double> mapElement = {m_theoDensity[i][0], newPot};
    	m_KSOutPot.insert(mapElement);
    }
	return;
}
*/





