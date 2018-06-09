# include "parameters.h"
# include "recap.hpp"
# include <fstream>

std::string RecapFile::m_recapFilePath = "Outputs/recap.txt";
RecapFile::RecapFile(double h, double loops, double time): m_h(h), m_loops(loops), m_time(time){}

void RecapFile::paramRecap()
{
	std::ofstream fOut(m_recapFilePath);
	fOut << "Element name: " << Parameters::ElementConstants::elementName() << std::endl;

	std::string nclType = "Protons";
	unsigned int numberType = Parameters::ElementConstants::NP();
	if(Parameters::NucleonType::isNeutron())
	{
		nclType = "Neutrons";
		numberType = Parameters::ElementConstants::NN();
	}

	fOut << "Nucleon type: " << nclType << std::endl;
	fOut << "Number of nucleons: " << numberType << std::endl;
    fOut << "Integration step size: " << m_h << " [fm]" << std::endl;
    fOut << "Initial radius: " << Parameters::IntegrationParameters::x0() << " [fm]" << std::endl;
    fOut << "Final radius: " << Parameters::IntegrationParameters::x1() << " [fm]" << std::endl;
    fOut << "Percentage of convergence " << Parameters::convperc << std::endl;
    fOut << "Convergence achieved after: " << m_loops << " iterations" << std::endl;
    fOut << "Convergence achieved in: " << m_time << " seconds" << std::endl;
    fOut.close();

    return;
}




