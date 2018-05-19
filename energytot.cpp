# include "parameters.h"
# include "energytot.hpp"
# include "eigenfunction.hpp"
# include "element.hpp"

/*=================================================================================================
 * Total Energy
 *===============================================================================================*/

void EnergyTOT::energyTot (const ElementEigenfunctions& elEigf, const ElementEigenvalues& elEigV, double h) const
{
	double factor = (-1)*(Parameters::hbarc*Parameters::hbarc)/(2*Parameters::mn);

	double nrgtot = 0;
	for (unsigned int i = 0; i < elEigf.size(); ++i)
	{
        std::map<double, double>::iterator it = elEigf[i].get().begin();
		double tjj = 0;
		std::map<double, double>::iterator p = it;
		++it;

		for (; it != elEigf[i].get().end(); ++it)
		{
			tjj += std::abs(((elEigf[i].get().at(p -> first)*p -> first*p -> first) +
                             (elEigf[i].get().at(it -> first)*it -> first*it -> first)))*h/2;
		}

		nrgtot += 0.5*(tjj + elEigV[i][3]);
	}

    return;
}




