# include "parameters.h"
# include "densities.h"
# include "eigenfunction.hpp"
# include <fstream>
# include <vector>
# include <algorithm>


/*========================================================
 * Theoretical density
 *======================================================*/

void NuclearDensity::theoreticalDensity(const ElementEigenfunctions& psi, const OrderedLevelDegeneration& degen)
{
    m_thDensity.clear();
    std::vector<Eigenfunction>::const_iterator el = psi.begin();
    std::vector<unsigned int>::const_iterator d = degen.begin();
    
    for (; el != psi.end() && d != degen.end(); ++el, ++d)
    {
        for (const auto& it : el -> get())
        {
            const double thdensity = 1./(4*Parameters::PI*(it.first*it.first))*(*d)*(it.second*it.second);
            m_thDensity[it.first] += thdensity;
        }
    }
	return;
}

/*==========================================================
 * Convergence condition for SOG densities
 *========================================================*/

bool NuclearDensity::hasConverged () const
{
	double maxDiff = std::abs(m_thDensity.begin() -> second - m_sogDensity.begin() -> second);
	double xMax = m_thDensity.begin() -> first;
    for (const auto& it : m_thDensity)
	{
		if(std::abs(it.second - m_sogDensity.at(it.first)) > maxDiff) //access only, throw exception if key is not found
		{
			maxDiff = std::abs(it.second - m_sogDensity.at(it.first));
			xMax = it.first;
		}
	}

	//const double epsilon = m_sogDensity.at(xMax)*0.05;
    //const double epsilon = 0.1*0.04;
    //m_distanceToConvergenge = maxDiff - epsilon;
    m_epsilon = m_sogDensity.at(xMax)*0.05;
    m_distanceToConvergenge = maxDiff - m_epsilon;

	return maxDiff < m_epsilon; // convergence condition
}

double NuclearDensity::distanceToConvergence() const
{
    return m_distanceToConvergenge;
}

double NuclearDensity::epsilon() const
{
    return m_epsilon;
}

Density NuclearDensity::getTheoreticalDensity() const
{
    return m_thDensity;
}

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
        return writeMap(outputQuery.second.getTheoreticalDensity(), wStream, false);
    
    if (oqs == "sogdensity" || oqs == "sog")
        return writeMap(outputQuery.second.getSOGDensity(), wStream, false);
    
    std::string usageMessage = "Possible candidates are: thdensity, sogdensity";
    return std::cerr << "NuclearDensity::operator<< : No matching results for " << outputQuery.first << std::endl
        << usageMessage << std::endl;
}

/*=================================================================
 * SOG density
 *===============================================================*/

void NuclearDensity::sogDensity(const std::vector<std::vector<double>>& QRparameters, double h)
{
    const double alpha = sqrt(2./3.)*Parameters::rp;
	const double gamma = sqrt(2./3.)*Parameters::rms;
	const double beta = sqrt((gamma*gamma)-(alpha*alpha));

	const unsigned long NSteps = static_cast<unsigned long>(std::abs(Parameters::x_fin - Parameters::x_in)/h);
	double radiusx = Parameters::x_in;
	for (unsigned int r = 0 ; r < NSteps + 1; ++r)
	{
        const double c1 = 1./(2*pow(Parameters::PI, 3./2.)*radiusx);

		double c3 = 0;
		for (unsigned int i = 0; i < QRparameters.size(); ++i)
		{
            const double Ai = (Parameters::NN*QRparameters[i][1])/(1 + (2.*QRparameters[i][0]*QRparameters[i][0])/(gamma*gamma));
            const double exp1 = exp((-1)*((radiusx - QRparameters[i][0])/beta)*((radiusx - QRparameters[i][0])/beta));
            const double exp2 = exp((-1)*((radiusx + QRparameters[i][0])/beta)*((radiusx + QRparameters[i][0])/beta));
            const double c2p1 = ((radiusx + QRparameters[i][0])/(beta*beta*beta)) - (QRparameters[i][0]/(beta*gamma*gamma));
            const double c2p2 = ((radiusx - QRparameters[i][0])/(beta*beta*beta)) + (QRparameters[i][0]/(beta*gamma*gamma));
            const double c2 = exp1*c2p1 + exp2*c2p2;
            
            c3 += Ai*c2;
		}
		const double sogdens = c1*c3;
		m_sogDensity[radiusx] = sogdens;
		radiusx += h;
	}

	// Normalization by trapezes method integration
	Density::iterator it = m_sogDensity.begin();
	//Density::iterator lastEval = m_sogDensity.end();
	//--lastEval;

	//integral end-point evaluations
	double scalar = 0;//( it -> second*it -> first*it -> first + lastEval -> second * lastEval -> first * lastEval -> first)*h/2.;
	//++it;
	Density::iterator p = it;
	++it;
    for (; it != m_sogDensity.end(); ++it)
    {
        scalar += ((p -> second*p -> first*p -> first) + (it -> second*it -> first*it -> first))*h/2;
        ++p;
    }

	double norm = Parameters::NN/scalar/4/Parameters::PI;

    for (auto& it : m_sogDensity)
    {
        it.second *= norm;
    }
    
	return;
}

Density NuclearDensity::getSOGDensity() const
{
    return m_sogDensity;
}

/*=========================================================================
 * Monte - Carlo Density
 *=======================================================================*/

void NuclearDensity::mcDensity(std::ifstream& inStream)
{
	readMap(m_mcDensity, inStream, false);
}

Density NuclearDensity::getMCDensity() const
{
    return m_mcDensity;
}

/*==========================================================================
 * Convergence condition for Monte-Carlo densities
 *========================================================================*/

bool NuclearDensity::hasConvergedMC () const
{
	// Load vectors from maps to avoid keys incompatibility
	m_theoDensity.clear();
	for (const auto& it : m_thDensity)
	{
		m_theoDensity.push_back(it.second);
	}

	m_MCDensity.clear();
	for (const auto& it : m_mcDensity)
	{
		m_MCDensity.push_back(it.second);
	}

	// Convergence condition implementation
	int i = 0;
	double maxDiff = std::abs(m_theoDensity[i] - m_MCDensity[i]);
	int iMax = i;
    for ( int i = 1; i < m_theoDensity.size(); ++i)
	{
		if(std::abs(m_theoDensity[i] - m_MCDensity[i]) > maxDiff)
		{
			maxDiff = std::abs(m_theoDensity[i] - m_MCDensity[i]);
			iMax = i;
		}
	}

    m_epsilonMC = m_MCDensity[iMax]*0.05;
    m_distanceToConvergengeMC = maxDiff - m_epsilonMC;

	return maxDiff < m_epsilonMC; // convergence condition
}

double NuclearDensity::distanceToConvergenceMC() const
{
    return m_distanceToConvergengeMC;
}

double NuclearDensity::epsilonMC() const
{
    return m_epsilonMC;
}

