# include "parameters.h"
# include "recap.hpp"
# include <time.h>
# include <fstream>

std::string RecapFile::m_recapFilePath = "Outputs/recap.txt";
RecapFile::RecapFile(double h, double loops, double time, std::string fSoG,
		std::string fPot, std::string fDens,  std::string fQNumb, int pType, int edType): m_h(h), m_loops(loops),
				m_time(time), m_fSoG(fSoG), m_fPot(fPot), m_fDens(fDens), m_fQNumb(fQNumb), m_pType(pType), m_edType(edType){}

void RecapFile::paramRecap()
{
	std::ofstream fOut(m_recapFilePath);

	time_t t = time(NULL);
	tm *time;
	time = localtime(&t);

	int day = time -> tm_mday;
	int month = time -> tm_mon; //month - 1
	int year = time -> tm_year;	//year - 1900
	int hour = time -> tm_hour;
	int min = time -> tm_min;
	int sec = time -> tm_sec;

	fOut << "Date/Time: " << day << "." << month + 1 << "." << year - 100 << " - " << hour << ":" << min << ":" << sec << std::endl;
	fOut << "Element name: " << Parameters::ElementConstants::elementName() << std::endl;

	std::string nclType = "Protons";
	unsigned int numberType = Parameters::ElementConstants::NP();
	if(Parameters::NucleonType::isNeutron())
	{
		nclType = "Neutrons";
		numberType = Parameters::ElementConstants::NN();
	}

	std::string potentialType, filepot;
	switch(m_pType)
	{
		case 0:
			potentialType = "Total";
			filepot = "none";
		break;
		case 1:
			potentialType = "Woods-Saxon";
			filepot = "none";
		break;
		case 2:
			potentialType = "Harmonic oscillator";
			filepot = "none";
		break;
		case 3:
			potentialType = "Test";
			filepot = "none";
		break;
		case 4:
			potentialType = "KS/Other";
			filepot = m_fPot;
		break;
	}

	std::string densityType, filedens;
	switch(m_edType)
	{
		case 0:
			densityType = "Protons SoG";
			filedens = m_fSoG;
		break;
		case 1:
			densityType = "Neutrons SoG";
			filedens = m_fSoG;
		break;
		case 2:
			densityType = "MC/Other";
			filedens = m_fDens;
	}

	fOut << "Nucleon type: " << nclType << std::endl;
	fOut << "Number of nucleons: " << numberType << std::endl;
	fOut << "Initial potential type: " << potentialType << std::endl;
	fOut << "Potential from file: " << filepot << std::endl;
	fOut << "Empirical density type: " << densityType << std::endl;
	fOut << "Empirical density from file: " << filedens << std::endl;
	fOut << "Orbitals from file: " << m_fQNumb << std::endl;
    fOut << "Integration step size: " << m_h << " [fm]" << std::endl;
    fOut << "Initial radius: " << Parameters::IntegrationParameters::x0() << " [fm]" << std::endl;
    fOut << "Final radius: " << Parameters::IntegrationParameters::x1() << " [fm]" << std::endl;
    fOut << "Percentage of convergence " << Parameters::convperc << std::endl;
    fOut << "Convergence achieved after: " << m_loops << " iterations" << std::endl;
    fOut << "Convergence achieved in: " << m_time << " seconds" << std::endl;
    fOut.close();

    return;
}




