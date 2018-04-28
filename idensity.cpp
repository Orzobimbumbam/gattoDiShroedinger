#include "idensity.h"
#include "parameters.h"
#include "eigenfunction.hpp"

NuclearDensity::NuclearDensity(): m_thDensity(), m_benchmarkDensity(), m_epsilon(0), m_distanceToConvergenge(0) {}

/*==========================================================================================================================
 * Theoretical density
 *========================================================================================================================*/

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

/*=========================================================================================================================
 * Convergence criteria
 *=======================================================================================================================*/

bool NuclearDensity::hasConverged () const
{
    double maxDiff = std::abs(m_thDensity.begin() -> second - m_benchmarkDensity.begin() -> second);
    double xMax = m_thDensity.begin() -> first;
    for (const auto& it : m_thDensity)
    {
        if(std::abs(it.second - m_benchmarkDensity.at(it.first)) > maxDiff) //access only, throw exception if key is not found
        {
            maxDiff = std::abs(it.second - m_benchmarkDensity.at(it.first));
            xMax = it.first;
        }
    }
    
    /*std::ofstream fOut("Outputs/distance.txt");
    double maxDiff = std::abs(1 - (m_thDensity.begin() -> second/m_benchmarkDensity.begin() -> second));
    double xMax = m_thDensity.begin() -> first;
    m_epsilon = (m_benchmarkDensity.begin() -> second)*0.05;
    for (const auto& it : m_thDensity)
	{
		if(std::abs(1 - (it.second/m_benchmarkDensity.at(it.first))) > maxDiff) //access only, throw exception if key is not found
		{
			maxDiff = std::abs(1 - (it.second/m_benchmarkDensity.at(it.first)));
			xMax = it.first;
			m_epsilon = (m_benchmarkDensity.at(it.first))*0.05;

		}

		m_distanceToConvergenge = maxDiff - m_epsilon;
		fOut << it.first << "\t\t" << it.second << "\t\t" << m_benchmarkDensity.at(it.first) << "\t\t" << maxDiff << "\t\t" <<
				m_epsilon << "\t\t" << m_distanceToConvergenge << std::endl;
	}
    fOut.close();*/

    //const double epsilon = m_sogDensity.at(xMax)*0.05;
    //const double epsilon = 0.1*0.04;
    //m_distanceToConvergenge = maxDiff - epsilon;
    m_epsilon = m_benchmarkDensity.at(xMax)*0.05;
    m_distanceToConvergenge = maxDiff - m_epsilon;
    
    return maxDiff < m_epsilon; // convergence condition
}

Density NuclearDensity::getTheoreticalDensity() const
{
    return m_thDensity;
}

double NuclearDensity::distanceToConvergence() const
{
    return m_distanceToConvergenge;
}

double NuclearDensity::epsilon() const
{
    return m_epsilon;
}








