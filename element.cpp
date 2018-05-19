#include "element.hpp"
#include "schroddy.h"
#include "parameters.h"
#include "eigenvalues.hpp"
#include <map>
#include <iomanip>
#include <sstream>

// Filling the orbitals
Element::Element(const OrderedOrbitalMatrix& orbitalMatrix)
{
    unsigned long nuclNum = 0;
    m_orbitalMatrixRows = 0;
    const unsigned long quantumAlphaSetSize = _getQuantumAlphaSetSize(orbitalMatrix);
    for (unsigned int i = 0; i < orbitalMatrix.size(); ++i) //loop through matrix rows
    {
        /*const unsigned int quantL = orbitalMatrix[i][1];
        const unsigned int degen = 2*(2*quantL + 1);
    	const unsigned int quantJ = orbitalMatrix[i][2];
    	const unsigned int degen = 2*quantJ + 1;*/
        const unsigned long degen = _getLevelDegeneration(orbitalMatrix, i, quantumAlphaSetSize);
        nuclNum += degen;
        ++m_orbitalMatrixRows;
        if (nuclNum >= Parameters::NN)
        {
            const unsigned long outerShellDegen = Parameters::NN - (nuclNum - degen); // number of nucleons in outer shell
            m_levelDegen.push_back(outerShellDegen);
            break;
        }
        else
            m_levelDegen.push_back(degen);
    }

    m_eigenValMatrix.resize(m_orbitalMatrixRows, std::vector<double>(quantumAlphaSetSize + 1));
}

unsigned long Element::_getQuantumAlphaSetSize(const OrderedOrbitalMatrix& orbitalMatrix) const
{
    if (!orbitalMatrix.empty())
        return orbitalMatrix[0].size();
    else
        return 0; // [Orzobimbumbam] : maybe should throw here?
}

unsigned long Element::_getLevelDegeneration(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex, unsigned long quantumAlphaSetSize) const
{
    const unsigned long maximumAlphaSetSize = 3; //[Orzobimbumbam] : nr, l, j
    unsigned long degeneration = 0;
    
    if (quantumAlphaSetSize == 0 || quantumAlphaSetSize > maximumAlphaSetSize)
        throw std::runtime_error ("Element::_getLevelDegeneration : invalid quantum number set size.");
    
    else if (quantumAlphaSetSize == 2)
        degeneration = 2*(2*orbitalMatrix[levelIndex][1] + 1);
    
    else if (quantumAlphaSetSize == 3)
        degeneration = 2*orbitalMatrix[levelIndex][2] + 1;
    
    else
        degeneration = 2*orbitalMatrix[levelIndex][1]*orbitalMatrix[levelIndex][1]; //[Orzobimbumbam] : double-check this
    
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
    	const double j = orbitalMatrix[i][2];
        const unsigned int l = orbitalMatrix[i][1];
        const unsigned int nr = orbitalMatrix[i][0];
        
        ptPot -> setL(l);
        ptPot -> setJ(j);
        
        const Schroddy tempSh(*ptPot, h);
        const GenericEigenvalues genEig(tempSh, nr, l);
        const double E = genEig.eigenvalue();
        const Eigenfunction eigf = tempSh.solveSchroddyByRK(x_in, x_fin, psi0(l), psiPrime0(l), E);
        elEigf.push_back(eigf);

    	m_eigenValMatrix[i][0] = nr;
    	m_eigenValMatrix[i][1] = l;
    	m_eigenValMatrix[i][2] = j;
    	m_eigenValMatrix[i][3] = E;
    }
    return elEigf;
}

ElementEigenValues Element::getLevelEigenvalue() const
{
    return m_eigenValMatrix;
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



