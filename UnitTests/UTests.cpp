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
BOOST_AUTO_TEST_CASE(shooting, *utf::tolerance(0.1))
{
    double H = 0.001;
    //int mass_num=Parameters::A;
    
    unsigned int k = 0; //same formula as yours, just keep this even
    unsigned int l_mom=0;
    unsigned int n = 2*k + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, n, l_mom);
    
    BOOST_TEST(GenEig.eigenvalue() == 15.0);
    
}
BOOST_AUTO_TEST_SUITE_END()
