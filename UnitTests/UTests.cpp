//
//  main.cpp
//  UnitTests
//
//  Created by Alberto Campi on 04/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#define BOOST_TEST_MODULE UTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#include "Includes.h"
#include <cmath>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;


BOOST_AUTO_TEST_SUITE(shroddy)
BOOST_AUTO_TEST_CASE(shootingHO10, *utf::tolerance(0.1))
{
    const double H = 0.001;
    //int mass_num=Parameters::A;
    
    const unsigned int k = 1;
    const unsigned int l_mom=0;
    const unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, k, l_mom);
    
    BOOST_TEST(GenEig.eigenvalue() == 3./2.*Parameters::hbar_omega);
    
}

BOOST_AUTO_TEST_CASE(shootingHO11, *utf::tolerance(0.1))
{
    const double H = 0.001;
    //int mass_num=Parameters::A;
    
    const unsigned int k = 1;
    const unsigned int l_mom = 1;
    const unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, k, l_mom);
    
    BOOST_TEST(GenEig.eigenvalue() == 5./2.*Parameters::hbar_omega);
    
}

BOOST_AUTO_TEST_CASE(shootingHO2012, *utf::tolerance(0.1))
{
    const double H = 0.001;
    //int mass_num=Parameters::A;
    
    unsigned int k = 2;
    unsigned int l_mom = 0;
    unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot10 (Parameters::mn, l_mom);
    Schroddy Sfunc10 (pot10, H);
    GenericEigenvalues GenEig10(Sfunc10, k, l_mom);
    
    BOOST_TEST(GenEig10.eigenvalue() == 7./2.*Parameters::hbar_omega);
    
    k = 1;
    l_mom = 2;
    n = 2*(k-1) + l_mom;
    
    HOPot pot02 (Parameters::mn, l_mom);
    Schroddy Sfunc02 (pot02, H);
    GenericEigenvalues GenEig02(Sfunc02, k, l_mom);
    
    BOOST_TEST(GenEig02.eigenvalue() == 7./2.*Parameters::hbar_omega);
}

BOOST_AUTO_TEST_CASE(shootingHO2113, *utf::tolerance(0.1))
{
    const double H = 0.001;
    //int mass_num=Parameters::A;
    
    unsigned int k = 2;
    unsigned int l_mom = 1;
    unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot21 (Parameters::mn, l_mom);
    Schroddy Sfunc21 (pot21, H);
    GenericEigenvalues GenEig21(Sfunc21, k, l_mom);
    
    BOOST_TEST(GenEig21.eigenvalue() == 9./2.*Parameters::hbar_omega);
    
    k = 1;
    l_mom = 3;
    n = 2*(k-1) + l_mom;
    
    HOPot pot13 (Parameters::mn, l_mom);
    Schroddy Sfunc13 (pot13, H);
    GenericEigenvalues GenEig13(Sfunc13, k, l_mom);
    
    BOOST_TEST(GenEig13.eigenvalue() == 9./2.*Parameters::hbar_omega);
}

BOOST_AUTO_TEST_CASE(shootingHORandom, *utf::tolerance(0.1))
{
    const double H = 0.001;
    //int mass_num=Parameters::A;
    
    const unsigned int k = rand()%21 + 1;
    const unsigned int l_mom = rand()%20;
    const unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, k, l_mom);
    
    BOOST_TEST(GenEig.eigenvalue() == (n + 3./2.)*Parameters::hbar_omega);
    
}






BOOST_AUTO_TEST_SUITE_END()
