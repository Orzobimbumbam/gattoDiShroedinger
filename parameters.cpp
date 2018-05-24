#include "parameters.h"
#include <sstream>
#include <iostream>

using namespace Parameters;

ElementConstants* ElementConstants::m_instancePtr = nullptr;

ElementConstants::ElementConstants(const std::string& fileName, const char delimiter)
{
    _initialiseParameters(fileName, delimiter);
}

void ElementConstants::_initialiseParameters(const std::string& fileName, const char delimiter)
{
    m_NP = 0; //optional value
    m_elementName = ""; //optional value
    std::stringstream ss(_extractFileNameKey(fileName));
    std::string key;
    int count = 0;
    while (std::getline(ss, key, delimiter))
    {
        ++count;
        if (count == 1)
            m_NN = std::stoi(key); //mandatory value
        if (count == 2)
            m_NP = std::stoi(key);
        if (count == 3)
            m_elementName = key;
    }
    
    if (count > 3)
        std::cerr << "initialiseParameters : Warning : file name may contain more keys than needed. " << std::endl;
    
    if (count == 0)
        throw std::runtime_error("initialiseParameters : invalid file name keys.");
    
    const unsigned int A = m_NN + m_NP;
    m_R0 = R0*pow(A,(1./3.));
    m_hBarOmega = 41./pow(A,(1./3.));
}

std::string ElementConstants::_extractFileNameKey(const std::string& fileName)
{
    std::stringstream ss(fileName);
    std::string fileNameKey;
    std::getline(ss, fileNameKey, '.');
    
    return fileNameKey;
}

void ElementConstants::initialiseElementConstants(const std::string& fileName, const char delimiter)
{
    if (ElementConstants::m_instancePtr == nullptr)
    {
        //ElementConstants::m_instance = true;
        ElementConstants::m_instancePtr = new ElementConstants(fileName, delimiter);
    }
}
/*
unsigned int Parameters::NN = _NN();
unsigned int Parameters::NP = _NP();
unsigned int Parameters::A = _A();*/

double psi0(unsigned int l)
{
    return pow(Parameters::x_in, l + 1);
}

double psiPrime0(unsigned int l)
{
	/*if (l == 0) return 0.;
	else*/ return (l + 1)*pow(Parameters::x_in, l);
}

