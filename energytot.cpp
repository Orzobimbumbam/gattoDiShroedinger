# include "parameters.h"
# include "energytot.hpp"
# include "eigenfunction.hpp"
# include "element.hpp"
# include <map>
# include <vector>

/*=================================================================================================
 * Total Energy
 *===============================================================================================*/

/*void EnergyTOT::energyTot (const Laplacians& LaplacVec, const ElementEigenvalues& elEigV, double h) const
{
	double factor = (-1)*(Parameters::hbarc*Parameters::hbarc)/(2*Parameters::mn);

	double nrgtot = 0;
	for (unsigned int i = 0; i < LaplacVec.size(); ++i)
	{
        std::map<double, double>::iterator it = LaplacVec[i].begin();
		double tjj = 0;
		std::map<double, double>::iterator p = it;
		++it;

		for (; it != LaplacVec[i].end(); ++it)
		{
			tjj += (LaplacVec[i].at(p -> first)*p -> first*p -> first) +
                             (LaplacVec[i].at(it -> first)*it -> first*it -> first)*h/2;
		}

		nrgtot += 0.5*(factor*tjj + elEigV[i][3]);
	}

    return;
}*/

double EnergyTOT::energyTot(const ElementEigenfunctions& elEigf, const ElementEigenvalues& elEigV, double h) const
{
	const double factor = -Parameters::hbarc*Parameters::hbarc/(2*Parameters::mn);
	double nrgtot = 0;

	for (unsigned int i = 0; i < elEigf.size(); ++i)
	{
		std::map<double, double>::const_iterator it = ++elEigf[i].get().begin();
		//++it;
		std::map<double, double>::const_iterator p = it;
		std::map<double, double>::const_iterator end = --elEigf[i].get().end();
		double tjj = 0;

		//--end;
		++it;

		for (; it != end; ++it)
		{
			tjj += p-> second*_laplacian(p) + it -> second*_laplacian(it)*h/2;
			++p;
		}

		nrgtot += 0.5*(factor*tjj + elEigV[i][3]);

	}

	return nrgtot;
}

double EnergyTOT::_laplacian(std::map<double, double>::const_iterator itt) const
{
	std::map<double, double>::const_iterator minusItt = --itt;
	std::map<double, double>::const_iterator plusItt = ++itt;
	const double h = itt -> first - minusItt -> first;
	const double laplacian = ((minusItt -> second) + (plusItt -> second) -2*(itt -> second))/(h*h);

	return laplacian;
}



// laplacians for each eigenfunctions
/*void laplacian(const ElementEigenfunctions& elEigf, double h) const
{
	Laplacians LaplacVec;
	Laplacian LaplacMap;
	LaplacVec.clear();
	for (unsigned int i = 0; i < elEigf.size(); ++i)
	{
		LaplacMap.clear();
		double abs1 = 0, abs2 = 0, abs3 = 0, laplac = 0;
		std::map<double, double>::iterator it = elEigf[i].get().begin();
		for (; it != elEigf[i].get().end(); ++it)
		{
			abs1 = std::abs(elEigf[i].get().at(it -> first))*std::abs(elEigf[i].get().at(it -> first));

			if(it -> first + h > elEigf[i].get().end() -> first)
				abs2 = std::abs(elEigf[i].get().end() -> second)*std::abs(elEigf[i].get().end() -> second);
			else
				abs2 = std::abs(elEigf[i].get().at(it -> first + h))*std::abs(elEigf[i].get().at(it -> first + h));

			if(it -> first - h < elEigf[i].get().begin() -> first)
				abs3 = std::abs(elEigf[i].get().begin() -> second)*std::abs(elEigf[i].get().begin() -> second);
			else
				abs3 = std::abs(elEigf[i].get().at(it -> first - h))*std::abs(elEigf[i].get().at(it -> first - h));

			laplac = (abs2 + abs3 -2*abs1)/(h*h);
			std::pair<double, double> elem = {it -> first, laplac};
			LaplacMap.insert(elem);

		}

		LaplacVec.push_back(LaplacMap);

	}
	return;
}*/

/*Laplacian EnergyTOT::get() const
{
    return m_LaplacMap;
}*/




