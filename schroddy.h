//
//  schroddy.h
//  Codice
//
//  Created by Alberto Campi on 02/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#pragma once
#include <vector>
#include <memory>
#include "initpot.h"
#include "eigenfunction.hpp"


/*=======================================================================
 * Shroedinger class + solver
 *=====================================================================*/

class Schroddy
{
public:
    Schroddy(const InitialPot& pot);
    Schroddy(const InitialPot& pot, double H);
    Schroddy(const Schroddy& sourceSh);
    Schroddy& operator=(const Schroddy& rhsSchroddy);
    
    double getH() const {return m_h;}
    void setH(double H) {m_h = H;}
    std::unique_ptr<InitialPot>& getInitialPotPtr() {return m_pot;};

    double solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, std::vector<double>& psiArray) const; //psiPrime0 is boundary condition on first derivative of eigenfunction
    
    const Eigenfunction solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E) const;

private:
    Schroddy();
    void swap(Schroddy& sourceSh);
    std::unique_ptr<InitialPot> m_pot;
    mutable double m_h;
};








