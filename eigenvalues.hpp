//
//  eigenvalues.hpp
//  Codice
//
//  Created by Alberto Campi on 10/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#ifndef eigenvalues_hpp
#define eigenvalues_hpp

/*==================================================================
 * Pure virtual class for eigenvalues
 *================================================================*/
#include "schroddy.h"

class Eigenvalues
{
public:
    virtual double eigenvalue() const = 0;
    virtual ~Eigenvalues();
    
    virtual Eigenvalues* clone() const = 0;
    
private:
    
};

/*===================================================================
 * HO eigenvalues generator class
 *=================================================================*/

class HarmonicEigenvalues: public Eigenvalues
{
public:
    HarmonicEigenvalues(unsigned int n, int l);
    double eigenvalue() const override;
    //setEnergyLevel(unsigned int n);
    //setAngularMomentum(int l);
    
    Eigenvalues* clone() const override;
    
private:
    HarmonicEigenvalues();
    
    //if you want to be able to change quantum numbers within same object, implement set methods and remove const keyword
    const unsigned int m_n;
    const int m_l;
};

/*=====================================================================
 * Eigenvalues generator (shooting method) class
 *===================================================================*/

//singleton class for trial eigenvalues
class TrialEigenvalues
{
public:
    static double getEigenval1();
    static double getEigenval2();
    
    
private:
    static TrialEigenvalues* m_trialEigenvalusObj;
    const double m_eigenval1, m_eigenval2;
    TrialEigenvalues();
};


class GenericEigenvalues: public Eigenvalues
{
public:
    GenericEigenvalues (const Schroddy& sh, unsigned int nState, unsigned int lState);
    double eigenvalue() const override;
    
    Eigenvalues* clone() const override;
    
private:
    GenericEigenvalues();
    double shootingMethod(double E1, double E2, unsigned int parity) const;
    
    const Schroddy m_sh;
    unsigned int m_nState, m_lState;
    
};



#endif /* eigenvalues_hpp */
