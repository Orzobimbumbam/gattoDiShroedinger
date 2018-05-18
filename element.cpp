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
    unsigned int nuclNum = 0;
    m_orbitalMatrixRows = 0;
    for (unsigned int i = 0; i < orbitalMatrix.size(); ++i) //loop through matrix rows
    {
        const unsigned int quantL = orbitalMatrix[i][1];
        const unsigned int degen = 2*(2*quantL + 1);
        nuclNum += degen;
        ++m_orbitalMatrixRows;
        if (nuclNum >= Parameters::NN)
        {
            const unsigned int outerShellDegen = Parameters::NN - (nuclNum - degen); // number of nucleons in outer shell
            m_levelDegen.push_back(outerShellDegen);
            break;
        }
        else
            m_levelDegen.push_back(degen);
    }

    m_eigenValMatrix.resize(m_orbitalMatrixRows, std::vector<double>(3)); //[Orzobimbumbam] hard-coded values should be avoided. If spin get added, column number will increase to 4 ?
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
        //[Orzobimbumbam] add logic for spin quantum number; remember that j = 0 means no spin-orbit interaction
        const unsigned int l = orbitalMatrix[i][1];
        const unsigned int nr = orbitalMatrix[i][0];
        ptPot -> setL(l);
        //[Orzobimbumbam] example: ptPot -> setJ(j);
        const Schroddy tempSh(*ptPot, h);
        const GenericEigenvalues genEig(tempSh, nr, l);
        const double E = genEig.eigenvalue();
        const Eigenfunction eigf = tempSh.solveSchroddyByRK(x_in, x_fin, psi0(l), psiPrime0(l), E);
        elEigf.push_back(eigf);

    	m_eigenValMatrix[i][0] = nr;
    	m_eigenValMatrix[i][1] = l;
    	m_eigenValMatrix[i][2] = E;
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



