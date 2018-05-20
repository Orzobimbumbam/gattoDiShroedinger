#pragma once
#include <vector>
#include <memory>
#include "initpot.h"
#include "eigenfunction.hpp"


/*=======================================================================
 * Schroedinger class + solver
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
    const std::unique_ptr<InitialPot>& getInitialPotPtr() const {return m_pot;};

    double solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E, std::vector<double>& psiArray) const; //psiPrime0 is boundary condition on first derivative of eigenfunction
    
    const Eigenfunction solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0, double E) const;

private:
    Schroddy();
    void swap(Schroddy& sourceSh);
    std::unique_ptr<InitialPot> m_pot;
    mutable double m_h;
    
    double _spinOrbitInteraction(double x) const;
};








