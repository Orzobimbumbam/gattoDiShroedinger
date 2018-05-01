//
//  UTests.cpp
//  UnitTests
//
//  Created by Alberto Campi on 04/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#define BOOST_TEST_MODULE UTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#include "../Codice/Includes.h"

#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <iterator>

//CALC UNIT TEST SPECIFICATIONS

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;
typedef std::vector<double> PsiArray;

//NB. set Rn = 5*Rn in parameters.h before running the test cases
BOOST_AUTO_TEST_SUITE(Calc)
/*
BOOST_AUTO_TEST_CASE(shootingWS10, *utf::tolerance(0.1))
{
    const double H = 0.1;
    //int mass_num=Parameters::A;
    
    const unsigned int k = 2;
    const unsigned int l_mom=3;
    //const unsigned int n = 2*(k-1) + l_mom;
    
    WSaxPot pot (Parameters::Rn, Parameters::a0, Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, k, l_mom);
    
    BOOST_TEST(GenEig.eigenvalue() == 4.93);
    
}*/

BOOST_AUTO_TEST_CASE(shootingHO10, *utf::tolerance(0.1))
{
    const double H = 0.1;
    //int mass_num=Parameters::A;
    
    const unsigned int k = 1;
    const unsigned int l_mom=0;
    //const unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, k, l_mom);
    
    BOOST_TEST(GenEig.eigenvalue() == 3./2.*Parameters::hbar_omega);
    
}

BOOST_AUTO_TEST_CASE(shooting_HO11, *utf::tolerance(0.1))
{
    const double H = 0.01;
    //int mass_num=Parameters::A;
    
    const unsigned int k = 1;
    const unsigned int l_mom = 1;
    //const unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, k, l_mom);
    
    BOOST_TEST(GenEig.eigenvalue() == 5./2.*Parameters::hbar_omega);
    
}

BOOST_AUTO_TEST_CASE(shooting_HO2012, *utf::tolerance(0.1))
{
    const double H = 0.001;
    //int mass_num=Parameters::A;
    
    unsigned int k = 2;
    unsigned int l_mom = 0;
    //unsigned int n = 2*(k-1) + l_mom;
    
    HOPot pot10 (Parameters::mn, l_mom);
    Schroddy Sfunc10 (pot10, H);
    GenericEigenvalues GenEig10(Sfunc10, k, l_mom);
    
    BOOST_TEST(GenEig10.eigenvalue() == 7./2.*Parameters::hbar_omega);
    
    k = 1;
    l_mom = 2;
    //n = 2*(k-1) + l_mom;
    
    HOPot pot02 (Parameters::mn, l_mom);
    Schroddy Sfunc02 (pot02, H);
    GenericEigenvalues GenEig02(Sfunc02, k, l_mom);
    
    BOOST_TEST(GenEig02.eigenvalue() == 7./2.*Parameters::hbar_omega);
}

BOOST_AUTO_TEST_CASE(shootingHO_2113, *utf::tolerance(0.1))
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

BOOST_AUTO_TEST_CASE(shooting_HORandom, *utf::tolerance(0.1))
{
    const double H = 0.1;
    //int mass_num=Parameters::A;
    
    unsigned int N = 10;
    for (unsigned int i = 0; i < N; ++i)
    {
        const unsigned int k = rand()%10 + 1;
        const unsigned int l_mom = rand()%10;
        const unsigned int n = 2*(k-1) + l_mom;
        
        HOPot pot (Parameters::mn, l_mom);
        Schroddy Sfunc (pot, H);
        GenericEigenvalues GenEig(Sfunc, k, l_mom);
        
        BOOST_TEST(GenEig.eigenvalue() == (n + 3./2.)*Parameters::hbar_omega);
    }
    
}

BOOST_AUTO_TEST_CASE(solveSchroddyByRK_HOPot1_Values)
{
    using namespace Parameters;
    const double H = 0.001;
    const unsigned int l_mom = rand()%10;
    const double E = rand()/static_cast<double>(RAND_MAX)*100;
    
    HOPot pot(Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    PsiArray psi;
    
    Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E, psi);
    const Eigenfunction eigf = Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E);
    
    PsiArrayKVP kvp = eigf.keyValues();
    PsiArray eigfValues;
    for (const auto& it : kvp)
        eigfValues.push_back(it.second);
    
    BOOST_CHECK(eigfValues == psi);
}

BOOST_AUTO_TEST_CASE(solveSchroddyByRK_HOPot1_KeyValues)
{
    using namespace Parameters;
    const double H = 0.001;
    const unsigned int l_mom = rand()%10;
    const double E = rand()/static_cast<double>(RAND_MAX)*100;
    
    HOPot pot(Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    PsiArray psi;
    
    Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E, psi);
    const Eigenfunction eigf = Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E);
    
    PsiArrayKVP kvp = eigf.keyValues();
    PsiArray eigfValues, eigfx, x;
    for (const auto& it : kvp)
    {
        eigfValues.push_back(it.second);
        eigfx.push_back(it.first);
    }
    
    x.push_back(x_in);
    for (unsigned long i = 1; i < std::abs(x_fin - x_in)/H; ++i)
        x.push_back(x[i-1] + H);
    
    BOOST_CHECK(eigfValues == psi && eigfx == x);
}

BOOST_AUTO_TEST_CASE(solveSchroddyByRK_HOPot2_Values)
{
    using namespace Parameters;
    const double H = 0.001;
    const unsigned int l_mom = rand()%10;
    const double E1 = rand()/static_cast<double>(RAND_MAX)*100;
    const double E2 = rand()/static_cast<double>(RAND_MAX)*100;
    
    HOPot pot(Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    PsiArray psi;
    
    Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E1, psi);
    const Eigenfunction eigf = Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E2);
    
    PsiArrayKVP kvp = eigf.keyValues();
    PsiArray eigfValues;
    for (const auto& it : kvp)
        eigfValues.push_back(it.second);
    
    if (E1 != E2)
        BOOST_CHECK(eigfValues != psi);
    else
        BOOST_CHECK(eigfValues == psi);
}

BOOST_AUTO_TEST_CASE(solveSchroddyByRK_HOPot3_Values)
{
    using namespace Parameters;
    const double H = 0.001;
    const unsigned int l_mom1 = rand()%10;
    const unsigned int l_mom2 = rand()%10;
    const double E = rand()/static_cast<double>(RAND_MAX)*100;
    
    HOPot pot1(Parameters::mn, l_mom1);
    HOPot pot2(Parameters::mn, l_mom2);
    Schroddy Sfunc1 (pot1, H);
    Schroddy Sfunc2 (pot2, H);
    PsiArray psi;
    
    Sfunc1.solveSchroddyByRK(x_in, x_fin, psi0(l_mom1), psiPrime0(l_mom1), E, psi);
    const Eigenfunction eigf = Sfunc2.solveSchroddyByRK(x_in, x_fin, psi0(l_mom2), psiPrime0(l_mom2), E);
    
    PsiArrayKVP kvp = eigf.keyValues();
    PsiArray eigfValues;
    for (const auto& it : kvp)
        eigfValues.push_back(it.second);
    
    if (l_mom1 != l_mom2)
        BOOST_CHECK(eigfValues != psi);
    else
        BOOST_CHECK(eigfValues == psi);
}


BOOST_AUTO_TEST_CASE(solveSchroddyByRK_WSaxPot_Values)
{
    using namespace Parameters;
    const double H = 0.001;
    const unsigned int l_mom = rand()%10;
    const double E = rand()/static_cast<double>(RAND_MAX)*100;
    
    WSaxPot pot(Rn, a0, Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    PsiArray psi;
    
    Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E, psi);
    const Eigenfunction eigf = Sfunc.solveSchroddyByRK(x_in, x_fin, psi0(l_mom), psiPrime0(l_mom), E);
    
    PsiArrayKVP kvp = eigf.keyValues();
    PsiArray eigfValues;
    for (const auto& it : kvp)
        eigfValues.push_back(it.second);
    
    BOOST_CHECK(eigfValues == psi);
}

struct KSInverseFixture : public KohnShamInverse
{
    KSInverseFixture(const KSPotential& extPotential) {setInternalMap(extPotential);};
    
    void setInternalMap(const KSPotential& extPotential) {m_KSOutPot = extPotential;};
};

struct PotOutFixture : public PotOut
{
    PotOutFixture(const KSInverseFixture& ksif): PotOut(ksif, 0, 0) {}; //implicit cast to base class object
    
    //KSInverseFixture m_ksif;
    double interpolate(double x) {return interpolatedPotential(x);};
};

BOOST_AUTO_TEST_CASE(PotOutInterpolate_straightLine)
{
    std::pair<double, double> a = std::make_pair(1, 2); //y = 2x
    std::pair<double, double> b = std::make_pair(2, 4);
    std::pair<double, double> c = std::make_pair(3, 6);
    std::pair<double, double> d = std::make_pair(4, 8);
    std::pair<double, double> e = std::make_pair(5, 10);
    KSPotential f = {a, b, c, d, e}; //fill in map with some function
    
    KSInverseFixture ksif(f);
    PotOutFixture poF(ksif);
    
    
    BOOST_CHECK(poF.interpolate(0.5) == 2.);
    BOOST_CHECK(poF.interpolate(7.) == 10.);
    BOOST_CHECK(poF.interpolate(3.5) == 7.);
}

BOOST_AUTO_TEST_CASE(PotOutInterpolate_cubic)
{
    std::pair<double, double> a = std::make_pair(-1, -1); //y = x^3
    std::pair<double, double> b = std::make_pair(-2, -8);
    std::pair<double, double> c = std::make_pair(0, 0);
    std::pair<double, double> d = std::make_pair(1, 1);
    std::pair<double, double> e = std::make_pair(2, 8);
    KSPotential f = {a, b, c, d, e}; //fill in map with some function
    
    KSInverseFixture ksif(f);
    PotOutFixture poF(ksif);
    
    const double expectedInterpolate = -9./2.;
    
    BOOST_CHECK(poF.interpolate(-1.5) == expectedInterpolate);
    BOOST_CHECK(poF.interpolate(0) == 0.);
    BOOST_CHECK(poF.interpolate(3) == 8.);
}



//more calcs test cases here...

BOOST_AUTO_TEST_SUITE_END()


//UTILS UNIT TEST SPECIFICATIONS

const std::string inputDir = "Inputs";
const std::string outputDir = "Outputs";

BOOST_AUTO_TEST_SUITE(Utils)
BOOST_AUTO_TEST_CASE(IOUtils_writeMap)
{
    
    std::pair<double, double> a = std::make_pair(1, 2); //y = 2x
    std::pair<double, double> b = std::make_pair(2, 4);
    std::pair<double, double> c = std::make_pair(3, 6);
    std::map<double, double> dMap = {a, b, c};
    
    std::ofstream out(outputDir + "/writeMapTestResult.txt");
    writeMap(dMap, out, true);
    
    std::ifstream inTestResult(outputDir + "/writeMapTestResult.txt");
    std::ifstream inExpResult(outputDir + "/writeMapExpResult.txt");
    
    std::istream_iterator<char> bTestResult(inTestResult), eTestResult;
    std::istream_iterator<char> bExpResult(inExpResult), eExpResult;
    
    BOOST_CHECK_EQUAL_COLLECTIONS(bTestResult, eTestResult, bExpResult, eExpResult);
}

BOOST_AUTO_TEST_CASE(IOUtils_readMap)
{
    
    std::pair<double, double> a = std::make_pair(1, 2); //y = 2x
    std::pair<double, double> b = std::make_pair(2, 4);
    std::pair<double, double> c = std::make_pair(3, 6);
    std::pair<double, double> d = std::make_pair(4, 8);
    std::pair<double, double> e = std::make_pair(5, 10);
    std::map<double, double> dMap = {a, b, c, d, e};
    
    std::ifstream in(inputDir + "/readMap.txt");
    std::map<double, double> rMap;
    readMap(rMap, in, true);
    
    BOOST_CHECK(rMap == dMap);
}

BOOST_AUTO_TEST_CASE(IOUtils_writeMatrix)
{
    std::vector<std::vector<double>> dMat = {{1,2,3,4}, {5,6,7,8}, {9,10,11,12}, {13,14,15,16}};
    
    std::ofstream out(outputDir + "/writeMatrixTestResult.txt");
    writeMatrix(dMat, out, false);
    
    std::ifstream inTestResult(outputDir + "/writeMatrixTestResult.txt");
    std::ifstream inExpResult(outputDir + "/writeMatrixExpResult.txt");
    
    std::istream_iterator<char> bTestResult(inTestResult), eTestResult;
    std::istream_iterator<char> bExpResult(inExpResult), eExpResult;
    
    BOOST_CHECK_EQUAL_COLLECTIONS(bTestResult, eTestResult, bExpResult, eExpResult);
}

BOOST_AUTO_TEST_CASE(IOUtils_readMatrix)
{
    std::vector<std::vector<double>> dMat = {{1,2,3,4}, {5,6,7,8}, {9,10,11,12}, {13,14,15,16}};
    
    
    std::ifstream in(inputDir + "/readMatrix.txt");
    std::vector<std::vector<double>> rMat(4);
    readMatrix(rMat, in, false);
    
    BOOST_CHECK(rMat == dMat);
}

//more Utils test cases here...

BOOST_AUTO_TEST_SUITE_END()
