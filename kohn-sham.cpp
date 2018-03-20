# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include "kohn-sham.h"
# include <vector>
# include <cmath>

/*===============================================================================
 * Kohm-Sham inverse equations
 *=============================================================================*/

KohnShamInverse::KohnShamInverse () {}

void KohnShamInverse::KSinverse(const std::vector<double>& inTheoDensity, const std::vector<double>& empiDensity,
    	const std::vector<double>& inPot)
{
	m_outPot.clear();
	unsigned long N = inPot.size();
	const double H = (Parameters::x_fin - Parameters::x_in)/N;
	double ratio1, ratio2, newPot, x = Parameters::x_in;
	for (int i = 0; i < N; ++i)
	{

		if (inPot[i] < 0)
		{
			ratio1 = (2 - (inTheoDensity[i]/empiDensity[i]));
			if (ratio1 < 1 - Parameters::pregamma) ratio1 = Parameters::pregamma;
			else if (ratio1 > 1 + Parameters::pregamma) ratio1 = Parameters::pregamma;
			newPot = ratio1*inPot[i];
		}
		else
		{
			ratio2 = (inTheoDensity[i]/empiDensity[i]);
			if (ratio2 < 1 - Parameters::pregamma) ratio2 = Parameters::pregamma;
			else if (ratio2 > 1 + Parameters::pregamma) ratio2 = Parameters::pregamma;
			newPot = ratio2*inPot[i];
		}

		//outPot.push_back(newPot);
		m_outPot[x] = newPot;

		x += H;
	}
	return;
}

void KohnShamInverse::getOutPot (std::map<double, double>& outPot) const
{
	outPot = m_outPot;
}






