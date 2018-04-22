#include "element.hpp"
#include "schroddy.h"
#include "parameters.h"
#include "initpot.h"
#include "eigenvalues.hpp"

#include <memory>
#include <map>
#include <fstream>

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
        if (nuclNum >= Parameters::A)
        {
            const unsigned int outerShellDegen = Parameters::A - (nuclNum - degen); // nucleons number in outer shell
            m_levelDegen.push_back(outerShellDegen);
            break;
        }
        else
            m_levelDegen.push_back(degen);
    }
    //ElementEigenValues temp(m_orbitalMatrixRows, std::vector<double>(3));
    //m_eigenValMatrix = temp;
    m_eigenValMatrix.resize(m_orbitalMatrixRows, std::vector<double>(3));
}

OrderedLevelDegeneration Element::getLevelDegeneration() const
{
    return m_levelDegen;
}

// Calculation of the eigenfunctions for each shell involved
ElementEigenfunctions Element::orbitalEigenfunction(const Schroddy& sh, const OrderedOrbitalMatrix& orbitalMatrix) const
{
    using namespace Parameters;
    ElementEigenfunctions elEigf;
    
    const std::unique_ptr<InitialPot> ptPot = sh.getInitialPotPtr() -> clone();
    const double h = sh.getH();
    
    for (unsigned int i = 0; i < m_orbitalMatrixRows; ++i)
    {
        const unsigned int l = orbitalMatrix[i][1];
        const unsigned int nr = orbitalMatrix[i][0];
        ptPot -> setL(l);
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

void writeElementEigenfunctions(const ElementEigenfunctions& elEigf, std::ostream& outStream)
{
    PsiArrayKVP kvp = elEigf[0].keyValues();
    for (int i = 0; i < elEigf[0].get().size(); ++i)
    {
        const double key = kvp[i].first;
        std::string rowEigf = std::to_string(key);
        for (const auto& it : elEigf)
        {
            rowEigf += "\t" + std::to_string(it.get()[key]);
        }
        outStream << rowEigf << std::endl;
    }
}



