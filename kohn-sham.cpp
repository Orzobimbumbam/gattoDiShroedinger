# include "parameters.h"
# include "kohn-sham.h"
# include "initpot.h"
# include "eigenfunction.hpp"
# include <vector>
# include <cmath>

/*==========================================================================================
 * Kohm-Sham inverse equations
 *========================================================================================*/

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

const KSPotential KohnShamInverse::getKSPot() const
{
    return m_KSOutPot;
}

KohnShamInverseWithLB::KohnShamInverseWithLB() : KohnShamInverse() {}
KohnShamInverseWithLB::KohnShamInverseWithLB(const InitialPot& iPot, double h) : KohnShamInverse(iPot, h) {}
KohnShamInverseWithJW::KohnShamInverseWithJW() : KohnShamInverse() {}
KohnShamInverseWithJW::KohnShamInverseWithJW(const InitialPot& iPot, double h) : KohnShamInverse(iPot, h) {}
KohnShamInverseWithWP::KohnShamInverseWithWP() : KohnShamInverse() {}
KohnShamInverseWithWP::KohnShamInverseWithWP(const InitialPot& iPot, double h, const ElementEigenvalues& eVal, const ElementEigenfunctions& inKSPsi) :
		KohnShamInverse(iPot, h), m_eVal(eVal), m_inKSPsi(inKSPsi) {}

/*===========================================================================================
 * Kohm-Sham inverse equations van Leeuwen-Baerends method
 *=========================================================================================*/

void KohnShamInverseWithLB::KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot)
{
	//I'm assuming here that the two maps have the same keys, so I can use one iterator only
    double ratio1, ratio2, newPot;
    for (const auto& it : density.getTheoreticalDensity())
    {
        const double alpha = 1.;//ratio1;
        if (inKSPot.getKSPot().at(it.first) < 0)
        {
        	ratio1 = 2 - alpha*it.second/density.getBenchmarkDensity().at(it.first);
        	if (ratio1 < 1 - Parameters::pregamma)
                ratio1 = 1 - Parameters::pregamma;
            											// Prefactor
        	else if (ratio1 > 1 + Parameters::pregamma)	// condition
                ratio1 = 1 + Parameters::pregamma;

            newPot = alpha*ratio1*inKSPot.getKSPot().at(it.first);
        }
        else
        {
            ratio2 = alpha*it.second/density.getBenchmarkDensity().at(it.first);
			if (ratio2 < 1 - Parameters::pregamma)
                ratio2 = 1 - Parameters::pregamma;		// Prefactor
														// condition
			else if (ratio2 > 1 + Parameters::pregamma)
                ratio2 = 1 + Parameters::pregamma;
            
            newPot = alpha*ratio2*inKSPot.getKSPot().at(it.first);
        }
        
        m_KSOutPot[it.first] = newPot;
    }
	return;
}


/*==============================================================================================
 * Kohm-Sham inverse equations Jensen-Wasserman method
 *============================================================================================*/

void KohnShamInverseWithJW::KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot)
{
	//I'm assuming here that the two maps have the same keys, so I can use one iterator only
    double ratio1, ratio2, newPot;
    for (const auto& it : density.getTheoreticalDensity())
    {
        const double alpha = 1;
        if (inKSPot.getKSPot().at(it.first) < 0)
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


/*==============================================================================================
 * Kohm-Sham inverse equations Wang-Parr method
 *============================================================================================*/

void KohnShamInverseWithWP::KSinverse(const NuclearDensity& density, const KohnShamInverse& inKSPot)
{
	//I'm assuming here that the two maps have the same keys, so I can use one iterator only
    for (const auto& it : density.getTheoreticalDensity())
    {
    	double sumPsi = 0;
        for (unsigned int i = 0; i < m_inKSPsi.size(); ++i)
        {
            const unsigned long energyIndex = m_eVal[i].size() - 1; //retrieve last column element index per each row
            sumPsi += m_inKSPsi[i].get().at(it.first)*m_inKSPsi[i].get().at(it.first)/m_eVal[i][energyIndex];
        }

    	const double newPot = ((density.getBenchmarkDensity().at(it.first) - it.second)/sumPsi) + inKSPot.getKSPot().at(it.first);
        m_KSOutPot[it.first] = newPot;
    }
	return;
}

void KohnShamInverseWithWP::setElementEigenfunctions(const ElementEigenfunctions &elEigf)
{
    m_inKSPsi = elEigf;
}

void KohnShamInverseWithWP::setElementEigenvalues(const ElementEigenvalues &elEigV)
{
    m_eVal = elEigV;
}









