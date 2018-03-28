//
//  element.cpp
//  Codice
//
//  Created by Alberto Campi on 27/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#include "element.hpp"
#include "eigenfunction.hpp"
#include "schroddy.h"
#include "parameters.h"
#include "initpot.h"
#include "eigenvalues.hpp"

#include <memory>

Element::Element(const OrderedOrbitalMatrix& orbitalMatrix)
{
    unsigned int nuclNum = 0;
    for (unsigned int i = 0; i < orbitalMatrix.size(); ++i) //loop through matrix rows
    {
        const unsigned int quantL = orbitalMatrix[i][1];
        const unsigned int degen = 2*(quantL + 1);
        nuclNum += degen;
        
        if (nuclNum > Parameters::A)
        {
            const unsigned int outerShellDegen = nuclNum - Parameters::A;
            m_levelDegen.push_back(outerShellDegen);
            break;
        }
        else
            m_levelDegen.push_back(degen);
    }
}

OrderedLevelDegeneration Element::getLevelDegeneration() const
{
    return m_levelDegen;
}

ElementEigenfunctions Element::orbitalEigenfunction(const Schroddy& sh, const OrderedOrbitalMatrix& orbitalMatrix) const
{
    using namespace Parameters;
    ElementEigenfunctions elEigf;
    
    const std::unique_ptr<InitialPot> ptPot = sh.getInitialPotPtr() -> clone();
    const double h = sh.getH();
    
    for (unsigned int i = 0; i < orbitalMatrix.size(); ++i)
    {
        const unsigned int l = orbitalMatrix[i][1];
        const unsigned int nr = orbitalMatrix[i][0];
        ptPot -> setL(l);
        const Schroddy tempSh(*ptPot, h);
        const GenericEigenvalues genEig(tempSh, nr, l);
        const double E  = genEig.eigenvalue();
        const Eigenfunction eigf = tempSh.solveSchroddyByRK(x_in, x_fin, psi0(l), psiPrime0(l), E);
        elEigf.push_back(eigf);
    }
    
    return elEigf;
}







