#include "parameters.h"
#include <sstream>
#include <iostream>

using namespace Parameters;

ElementConstants* ElementConstants::m_instancePtr = nullptr;

/*=====================================================================================================
 * Load element name, neutron number, proton number and rms value from filename
 *===================================================================================================*/

ElementConstants::ElementConstants(const std::string& fileName, const char delimiter)
{
    _initialiseParameters(fileName, delimiter);
}

void ElementConstants::_initialiseParameters(const std::string& fileName, const char delimiter)
{
    std::stringstream ss(_extractFileNameKey(fileName));
    std::string key;
    int count = 0;
    const unsigned int fileNameKeys = 4;
    while (std::getline(ss, key, delimiter))
    {
        ++count;
        if (count == 1) //all mandatory values, key must be complete
        	m_elementName = key;
        if (count == 2)
        	m_NN = std::stoi(key);
        if (count == 3)
        	m_NP = std::stoi(key);
        if (count == 4)
            m_rms = std::stod(key);
    }
    
    if (count > fileNameKeys)
        std::cerr << "initialiseParameters : Warning : file name may contain more keys than needed. " << std::endl;
    
    if (count < fileNameKeys)
        throw std::invalid_argument("initialiseParameters : invalid or incomplete file name keys.");
    
    // Parameters
    const unsigned int A = m_NN + m_NP;	// Mass number
    m_Rn = R0*pow(A,(1./3.));			// Nucleus radius
    m_hBarOmega = 41./pow(A,(1./3.));	// [MeV*fm]
}

/*=====================================================================================================
 * Parameters extractor from filename
 *===================================================================================================*/

std::string ElementConstants::_extractFileNameKey(const std::string& fileName)
{
    const char delimiter = '.';
    std::stringstream ss(fileName);
    std::string fileNameKey;
    do
    {
        std::string token;
        std::getline(ss, token, delimiter);
        
        if(ss.peek() != EOF)
            fileNameKey += (token + delimiter);
        else
        {
            fileNameKey.pop_back();
            break;
        }
        
    } while (true);
    
    //std::getline(ss, fileNameKey, '.');
    
    return fileNameKey;
}

void ElementConstants::initialiseElementConstants(const std::string& fileName, const char delimiter)
{
    if (ElementConstants::m_instancePtr == nullptr)
    {
        ElementConstants::m_instancePtr = new ElementConstants(fileName, delimiter);
    }
}

/*=====================================================================================================
 * Integration parameters
 *===================================================================================================*/

IntegrationParameters* IntegrationParameters::m_instancePtr = nullptr;
IntegrationParameters::IntegrationParameters(double x0, double x1) : m_x0(x0), m_x1(x1) {}

void IntegrationParameters::initialiseIntegrationParameters(double x0, double x1)
{
    if (IntegrationParameters::m_instancePtr == nullptr)
        IntegrationParameters::m_instancePtr = new IntegrationParameters(x0, x1);
}

/*=====================================================================================================
 * Nucleon type
 *===================================================================================================*/

NucleonType* NucleonType::m_instancePtr = nullptr;
NucleonType::NucleonType(bool isNeutron) : m_isNeutron(isNeutron) {}

void NucleonType::initialiseNucleonType(bool isNeutron)
{
    if (NucleonType::m_instancePtr == nullptr)
        m_instancePtr = new NucleonType(isNeutron);
}

/*=====================================================================================================
 * Spin-Orbit ON/OFF
 *===================================================================================================*/

SpinOrbit* SpinOrbit::m_instancePtr = nullptr;
SpinOrbit::SpinOrbit(bool SpinOrbitOn) : m_SpinOrbitOn(SpinOrbitOn) {}

void SpinOrbit::initialiseSpinOrbit(bool SpinOrbitOn)
{
    if (SpinOrbit::m_instancePtr == nullptr)
        m_instancePtr = new SpinOrbit(SpinOrbitOn);
}

/*=====================================================================================================
 * Boundary conditions for Schrodinger solver
 *===================================================================================================*/

// eigenfunction condition in r=0
double psi0(unsigned int l)
{
    return pow(Parameters::IntegrationParameters::x0(), l + 1);
}

// first derivative eigenfunction in r=0
double psiPrime0(unsigned int l)
{
	/*if (l == 0) return 0.;
	else*/ return (l + 1)*pow(Parameters::IntegrationParameters::x0(), l);
}

