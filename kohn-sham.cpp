# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include "kohn-sham.h"
# include <cmath>

/*=============================================
 * Kohm-Sham inverse equations
 *===========================================*/

KohnShamInverse::KohnShamInverse () {}

void KSinverse(const std::vector<double>& inTheoDensity, const std::vector<double>& empiDensity,
    	const std::vector<double>& inPot, std::vector<double>& outPot) const
{
	outPot.clear();
	for (int i = 0; i < inPot.size(); ++i)
	{
		if (inPot[i] < 0)
			double newPot = (2 - (inTheoDensity[i]/empiDensity[i]))*inPot[i];
		else
			double newPot = (inTheoDensity[i]/empiDensity[i])*inPot[i];

		outPot.push_back(newPot);
	}
	return;
}






