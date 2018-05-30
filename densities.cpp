# include "parameters.h"
# include "densities.h"
# include "eigenfunction.hpp"
# include <fstream>
# include <algorithm>


/*=============================================================
 * Class operators (friend, global)
 *===========================================================*/

std::ostream& operator<<(std::ostream& wStream, const Density& density)
{
    return writeMap(density, wStream, false);
}

std::ostream& operator<<(std::ostream& wStream, const NuclearDensityOutputQuery& outputQuery)
{
    std::string oqs = outputQuery.first;
    std::transform(oqs.begin(), oqs.end(), oqs.begin(), ::tolower); //transform string to lower case
    
    if (oqs == "thdensity" || oqs == "theoreticaldensity" || oqs == "theoretical")
        return writeMap(outputQuery.second -> getTheoreticalDensity(), wStream, false);
    
    if (oqs == "sogdensity" || oqs == "sog")
        return writeMap(outputQuery.second -> getBenchmarkDensity(), wStream, false);
    
    std::string usageMessage = "Possible candidates are: thdensity, sogdensity";
    return std::cerr << "NuclearDensity::operator<< : No matching results for " << outputQuery.first << std::endl
        << usageMessage << std::endl;
}

/*=========================================================================================================================================
 * SOG density
 *=======================================================================================================================================*/

void NuclearDensityWithSOG::benchmarkDensity(const std::vector<std::vector<double>>& QRparameters, double h)
{
    using namespace Parameters;
    
    const double alpha = sqrt(2./3.)*Parameters::rp;
	const double gamma = sqrt(2./3.)*ElementConstants::rms();
	const double beta = sqrt((gamma*gamma)-(alpha*alpha));

	const unsigned long NSteps = static_cast<unsigned long>(std::abs(IntegrationParameters::x1() - IntegrationParameters::x0())/h);
	double radiusx = IntegrationParameters::x0();
	for (unsigned int r = 0 ; r < NSteps + 1; ++r)
	{
        const double c1 = 1./(2*pow(Parameters::PI, 3./2.)*radiusx);

		double c3 = 0;
		for (unsigned int i = 0; i < QRparameters.size(); ++i)
		{
            const double Ai = (ElementConstants::NP()*QRparameters[i][1])/(1 + (2.*QRparameters[i][0]*QRparameters[i][0])/(gamma*gamma));
            const double exp1 = exp((-1)*((radiusx - QRparameters[i][0])/beta)*((radiusx - QRparameters[i][0])/beta));
            const double exp2 = exp((-1)*((radiusx + QRparameters[i][0])/beta)*((radiusx + QRparameters[i][0])/beta));
            const double c2p1 = ((radiusx - QRparameters[i][0])/(beta*beta*beta)) + (QRparameters[i][0]/(beta*gamma*gamma));
            const double c2p2 = ((radiusx + QRparameters[i][0])/(beta*beta*beta)) - (QRparameters[i][0]/(beta*gamma*gamma));
            const double c2 = exp1*c2p1 + exp2*c2p2;
            
            c3 += Ai*c2;
		}

		double ratio = 1;
		if(Parameters::nucleons == 1)
			ratio = ElementConstants::NN()/ElementConstants::NP();

		const double sogdens = ratio*c1*c3;
		m_benchmarkDensity[radiusx] = sogdens;
		radiusx += h;
	}

	// Normalization by trapezoid method integration
	/*Density::iterator it = m_benchmarkDensity.begin();
	double scalar = 0;
	Density::iterator p = it;
	++it;
    
    for (; it != m_benchmarkDensity.end(); ++it)
    {
        scalar += ((p -> second*p -> first*p -> first) + (it -> second*it -> first*it -> first))*h/2;
        ++p;
    }

	double norm = ElementConstants::NP()/scalar/4/Parameters::PI;

    for (auto& it : m_benchmarkDensity)
    {
        it.second *= norm;

    }*/

	// Control on SOG density normalization: if normalized the control returns the number of nucleons
	/*Density::iterator itt = m_benchmarkDensity.begin();
    double normControl = 0;
	Density::iterator pp = itt;
	++itt;
    for (; itt != m_benchmarkDensity.end(); ++itt)
    {
    	normControl += 4*Parameters::PI*((pp -> second*pp -> first*pp -> first) + (itt -> second*itt -> first*itt -> first))*h/2;
        ++pp;
    }
    std::cout << "SOG normalization control - Number of nucleons: " << normControl << std::endl;*/
    
	return;
}

Density NuclearDensityWithSOG::getBenchmarkDensity() const
{
    return m_benchmarkDensity;
}

/*=========================================================================================================================================
 * Neutrons of 208Pb SOG density by Zenihiro
 *=======================================================================================================================================*/

void NuclearDensityWithNeutronsSOG::benchmarkDensity(const std::vector<std::vector<double>>& QRparameters, double h)
{
	using namespace Parameters;

	const double gamma = sqrt(2./3.)*ElementConstants::rms();

	const unsigned long NSteps = static_cast<unsigned long>(std::abs(IntegrationParameters::x1() - IntegrationParameters::x0())/h);
	double radiusx = IntegrationParameters::x0();
	for (unsigned int r = 0 ; r < NSteps + 1; ++r)
	{
		const double c1 = ElementConstants::NN()/(2*pow(Parameters::PI, 3./2.)*gamma*gamma*gamma);

		double c3 = 0;
		for (unsigned int i = 0; i < QRparameters.size(); ++i)
		{
            const double exp1 = exp((-1)*(radiusx - QRparameters[i][0])*(radiusx - QRparameters[i][0])/(gamma*gamma));
            const double exp2 = exp((-1)*(radiusx + QRparameters[i][0])*(radiusx + QRparameters[i][0])/(gamma*gamma));
            const double c2 = (QRparameters[i][1])/(1 + (2.*QRparameters[i][0]*QRparameters[i][0])/(gamma*gamma));

            c3 += c2*(exp1 + exp2);
		}

		const double sogdens = c1*c3;
		m_benchmarkDensity[radiusx] = sogdens;
		radiusx += h;
	}

	// Control on SOG density normalization: if normalized the control returns the number of nucleons
	/*Density::iterator itt = m_benchmarkDensity.begin();
    double normControl = 0;
	Density::iterator pp = itt;
	++itt;
    for (; itt != m_benchmarkDensity.end(); ++itt)
    {
    	normControl += 4*Parameters::PI*((pp -> second*pp -> first*pp -> first) + (itt -> second*itt -> first*itt -> first))*h/2;
        ++pp;
    }
    std::cout << "SOG normalization control - Number of nucleons: " << normControl << std::endl;*/

	return;
}

Density NuclearDensityWithNeutronsSOG::getBenchmarkDensity() const
{
    return m_benchmarkDensity;
}


