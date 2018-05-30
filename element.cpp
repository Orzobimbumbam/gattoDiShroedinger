#include "element.hpp"
#include "schroddy.h"
#include "parameters.h"
#include "eigenvalues.hpp"

#include <map>
#include <iomanip>
#include <sstream>

const unsigned long Element::maximumAlphaSetSize = 3; // nr, l, j
const unsigned long Element::minimumAlphaSetSize = 2; // nr, l

// Filling the orbitals
Element::Element(const OrderedOrbitalMatrix& orbitalMatrix)
{
    int nuclType = Parameters::ElementConstants::NP();
    if(Parameters::NucleonType::isNeutron())
        nuclType = Parameters::ElementConstants::NN();

    unsigned long nuclNum = 0;
    m_orbitalMatrixRows = 0;
    const unsigned long quantumAlphaSetSize = _getQuantumAlphaSetSize(orbitalMatrix, 0);
    
    for (unsigned int i = 0; i < orbitalMatrix.size(); ++i) //loop through matrix rows
    {
        const unsigned long degen = _getLevelDegeneration(orbitalMatrix, i);
        nuclNum += degen;
        ++m_orbitalMatrixRows;
        if (nuclNum >= nuclType)
        {
            const unsigned long outerShellDegen = nuclType - (nuclNum - degen); // number of nucleons in outer shell
            m_levelDegen.push_back(outerShellDegen);
            break;
        }
        else
            m_levelDegen.push_back(degen);
    }

    m_eigenvalMatrix.resize(m_orbitalMatrixRows, std::vector<double>(quantumAlphaSetSize + 1)); 
}

unsigned long Element::_getQuantumAlphaSetSize(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex) const
{
    if (!orbitalMatrix.empty())
        return orbitalMatrix[0].size();
    else
        return 0;
}

unsigned long Element::_getLevelDegeneration(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex) const
{
    unsigned long degeneration = 0;
    const unsigned long quantumAlphaSetSize = _getQuantumAlphaSetSize(orbitalMatrix, levelIndex);
    
    if (quantumAlphaSetSize < Element::minimumAlphaSetSize || quantumAlphaSetSize > Element::maximumAlphaSetSize)
        throw std::runtime_error ("Element::_getLevelDegeneration : invalid quantum number set size.");
    
    else if (quantumAlphaSetSize == Element::minimumAlphaSetSize)
        degeneration = 2*(2*orbitalMatrix[levelIndex][1] + 1);
    
    else if (quantumAlphaSetSize == Element::maximumAlphaSetSize)
        degeneration = 2*orbitalMatrix[levelIndex][2] + 1;
    /*
    else
        degeneration = 2*orbitalMatrix[levelIndex][1]*orbitalMatrix[levelIndex][1]; //[Orzobimbumbam] : can be excluded, however it may result in undefined behaviour; should we throw?*/
    
    return degeneration;
}

OrderedLevelDegeneration Element::getLevelDegeneration() const
{
    return m_levelDegen;
}

// Compute eigenfunctions for each shell
ElementEigenfunctions Element::orbitalEigenfunction(const Schroddy& sh, const OrderedOrbitalMatrix& orbitalMatrix) const
{
    using namespace Parameters;
    ElementEigenfunctions elEigf;
    
    const std::unique_ptr<InitialPot> ptPot = sh.getInitialPotPtr() -> clone();
    const double h = sh.getH();
    
    for (unsigned int i = 0; i < m_orbitalMatrixRows; ++i)
    {
        const double j = _checkAndGetJ(orbitalMatrix, i);
        const unsigned int l = orbitalMatrix[i][1];
        const unsigned int nr = orbitalMatrix[i][0];
        
        ptPot -> setL(l);
        ptPot -> setJ(j);
        
        const Schroddy tempSh(*ptPot, h);
        const GenericEigenvalues genEig(tempSh, nr, l);
        const double E = genEig.eigenvalue();
        const Eigenfunction eigf = tempSh.solveSchroddyByRK(IntegrationParameters::x0(), IntegrationParameters::x1(), psi0(l), psiPrime0(l), E);
        elEigf.push_back(eigf);
        
        const LevelTuple eigenvalMatrixRow = std::make_tuple(nr, l, j, E); // {nr, l, j, E};
        /*
    	m_eigenvalMatrix[i][0] = nr;
    	m_eigenvalMatrix[i][1] = l;
    	m_eigenvalMatrix[i][2] = j;
    	m_eigenvalMatrix[i][3] = E;*/
        _addTuple(eigenvalMatrixRow, i);
    }
    return elEigf;
}

ElementEigenvalues Element::getLevelEigenvalue() const
{
    return m_eigenvalMatrix;
}

double Element::_checkAndGetJ(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex) const
{
    if (_getQuantumAlphaSetSize(orbitalMatrix, levelIndex) == Element::maximumAlphaSetSize) //check if j is in the orbitals matrix
        return orbitalMatrix[levelIndex][2];
    
    else
        return 0.0;
}

void Element::_addTuple(const LevelTuple& row, unsigned long levelIndex) const
{
    const unsigned long nColumns = m_eigenvalMatrix[levelIndex].size();
    const double j = std::get<2>(row);
    const double E = std::get<3>(row);
    
    m_eigenvalMatrix[levelIndex][nColumns - 1] = E;
    m_eigenvalMatrix[levelIndex][0] = std::get<0>(row);
    m_eigenvalMatrix[levelIndex][1] = std::get<1>(row);
    
    if (j != 0)
        m_eigenvalMatrix[levelIndex][2] = j;
}

// Write all eigenfunctions to a single file
void writeElementEigenfunctions(const ElementEigenfunctions& elEigf, std::ostream& outStream)
{
    PsiArrayKVP kvp = elEigf[0].keyValues();
    for (int i = 0; i < elEigf[0].get().size(); ++i)
    {
        const short conversionPrecision = 10;
        const double key = kvp[i].first;
        std::string rowEigf = std::to_string(key);
        for (const auto& it : elEigf)
        {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(conversionPrecision) << it.get()[key];
            rowEigf += "\t" + ss.str(); //std::to_string(it.get()[key]);
        }
        outStream << rowEigf << std::endl;
    }
}



