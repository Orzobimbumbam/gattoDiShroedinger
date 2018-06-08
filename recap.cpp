# include "parameters.h"
# include "recap.hpp"
# include <fstream>


RecapFile::RecapFile(double h, double xin, double xfin): m_h(h), m_xin(xin), m_xfin(xfin){}

void RecapFile::paramRecap()
{
	std::ofstream fOut("/Outputs/recap.txt", "a");
	fOut << "Element name: " << Parameters::ElementConstants::elementName() << std::endl;

	unsigned int nuclType = Parameters::ElementConstants::NP();
	if(Parameters::NucleonType::isNeutron())
		nuclType = Parameters::ElementConstants::NN();

	fOut << "Nucleons type: " << nuclType << std::endl;
    fOut << "Integration pass: " << m_h << " [fm]" << std::endl;
    fOut << "Initial radius: " << m_xin << " [fm]" << std::endl;
    fOut << "Final radius: " << m_xfin << " [fm]" << std::endl;
    fOut << "Percentage of convergence " << Parameters::convperc << std::endl;



}




