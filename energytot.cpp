# include "parameters.h"
# include "energytot.hpp"
# include "eigenfunction.hpp"
# include "element.hpp"
# include <map>
# include <vector>

/*=================================================================================================
 * Total Energy
 *===============================================================================================*/

double EnergyTOT::energyTot(const ElementEigenfunctions& elEigf, const ElementEigenvalues& elEigV) const
{
	const double factor = -Parameters::hbarc*Parameters::hbarc/(2*Parameters::mn);
	double nrgtot = 0;

	for (unsigned int i = 0; i < elEigf.size(); ++i)
	{
		std::map<double, double>::const_iterator it = ++elEigf[i].get().begin();
		std::map<double, double>::const_iterator p = it;
		std::map<double, double>::const_iterator end = --elEigf[i].get().end();
		double tjj = 0;

		++it;

		for (; it != end; ++it)
		{
            const double h = std::abs(it -> first - p -> first);
			tjj += p-> second*_laplacian(p, h) + it -> second*_laplacian(it, h)*h/2;
			++p;
		}

		nrgtot += 0.5*(factor*tjj + elEigV[i][3]);

	}

	return nrgtot;
}

double EnergyTOT::_laplacian(std::map<double, double>::const_iterator itt, double h) const
{
	std::map<double, double>::const_iterator minusItt = --itt;
	std::map<double, double>::const_iterator plusItt = ++itt;
	//const double h = itt -> first - minusItt -> first;
	const double laplacian = ((minusItt -> second) + (plusItt -> second) -2*(itt -> second))/(h*h);

	return laplacian;
}





