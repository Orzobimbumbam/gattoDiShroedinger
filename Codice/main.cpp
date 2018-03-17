//
//  main.cpp
//  Codice
//
//  Created by Alberto Campi on 17/01/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "Includes.h"


typedef std::vector<double>::const_iterator cwalk;
typedef std::vector<double>::iterator walk;

int main(int argc, const char * argv[])
{
    
    const double H = 0.01;
    int massNumb=Parameters::A;
    
    /*==========================================================================================
    * STEP-1
    * Calculate, by Schrodinger solver, the eigenfunctions for each quantum involved state
    *========================================================================================*/
        
    //Load the matrix with quantum numbers of each state from "orbitals.txt"
    std::vector <std::vector <int>> orbitals (36);
    std::fstream in ("../orbitals.txt", std::ios::in);
    
    for (int i = 0; i < 36; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            int temp;
            in >> temp;
            orbitals[i].push_back(temp);
        }
    }
    
    /*for (int i=0; i<36; ++i)
     {
     for (int j=0; j<2; ++j)
     std::cout << orbitals[i][j];
     std::cout << std::endl;
     }*/
    
    // Solve Schrodinger equation for each involved state
    int i = 0;
    int degen = 0;
    unsigned int nuclNumb = 0;
    const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/H;
    std::vector <double> psiArray;
    std::vector<double> thdensArray(NSteps + 1);
    
    std::ofstream file;
    file.open("Outputs/eigenfunc.txt");
    
    
    while (massNumb > 0)
    {
        const unsigned int quantNr = orbitals[i][0];
        const unsigned int quantL = orbitals[i][1];
        const unsigned int quantN = 2*(quantNr - 1) + quantL;
        HOPot pot (Parameters::mn, quantL);
        Schroddy Sfunc (pot, H);
        GenericEigenvalues GenEig(Sfunc, quantNr, quantL);
        const double eigval = GenEig.eigenvalue();
        Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(quantL),
                                psiPrime0(quantL), eigval , psiArray);
        
        degen = 2*(2*quantL+1); // orbital degeneration
        massNumb -= degen;
        
        if (massNumb < 0)
            nuclNumb = massNumb + degen; // nucleons number on last orbital in not filled case
        else nuclNumb = degen; // nucleons number in filled orbital
        
        // Claculate theoretical density for involved state
        Theoreticaldensity densy;
        densy.density(psiArray,thdensArray,nuclNumb,H);
        
        cwalk wk = psiArray.begin();
        while (wk != psiArray.end())
        {
            file << quantN << "\t" << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << *wk << std::endl;
            ++wk;
        }
        
        /*evalArray.push_back(eigval);
         efuncArray.push_back(eigfunc);
         nuclnumbArray.push_back(nuclNumb);*/
        
        //dstd::cout << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << eigval << "\t" << eigfunc << std::endl;
        
        ++i;
        //std::cout << massNumb << std::endl;
    }
    file.close();
    
    file.open("Outputs/density.txt");
    double rad = Parameters::x_in;
    cwalk wk = thdensArray.begin();
    while (wk != thdensArray.end())
    {
        file << rad << "\t"<< *wk << std::endl;
        ++wk;
        rad = rad + H;
    }
    
    file.close();
    
    /*=========================================================================================
     * STEP-2
     * Calculate the theoretical density
     *=======================================================================================*/
    
    // std::vector <std::vector <double> > efunctions (63500);
    //std::vector <int> data (4);
    //std::fstream in ("eigenfunc.txt", std::ios::in);
    //in.open ("eigenfunc.txt");
    /*std::fstream in ("eigenfunc.txt", std::ios::in);
     
     for (int i = 0; i < 63500; ++i)
     {
     for (int j = 0; j < 2; ++j)
     {
 			 double temp;
 			 in >> temp;
 			 efunctions[i].push_back(temp);
     }
     //if (in.eof()) break;
     }*/
    
    /*for (int i=0; i<63500; ++i)
     {
 				for (int j=0; j<2; ++j)
     std::cout << efunctions[i][j];
     std::cout << std::endl;
     }*/
    
    
    /*std::vector<double> thDensArray;
     std::vector<double> xArray;
     std::ofstream file2("density.csv");
     const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/H;
     
     int deg = 0;
     double efunc = 0, radiusx = Parameters::x_in;
     for (unsigned long i = 0; i < NSteps + 1; ++i)
     {
     deg = efunctions[i][0];
     efunc = efunctions[i][1];
     Theoreticaldensity densy;
     double thdensity = densy.density(efunc,deg,radiusx);
     thDensArray.push_back(thdensity);
     xArray.push_back(radiusx);
     radiusx += H;
     }
     
     std::vector<double>::iterator walk2 = thDensArray.begin();
     std::vector<double>::iterator walk3 = xArray.begin();
     while (walk2 != thDensArray.end() && walk3 != xArray.end())
     {
     file2 << *walk3 << ","<< *walk2 << std::endl;
     walk2++, walk3++;
     }
     
     file2.close();*/
    
    /*=========================================================================================
     * STEP-3
     * Calculate empirical densities from MC simulations or scattering
     *=======================================================================================*/
    
    
    
    /*===============================================================================
     * TEST - Calculate first 12 levels eigenvalues for HO potential
     *==============================================================================*/
    
    /*unsigned int nr=1;
     int l_mom=2;
     unsigned int n = 2*(nr-1) + l_mom;*/
    /*std::vector <double> psiArray;
     unsigned int nmax = 5;
     
     for (unsigned int n = 0; n <= nmax; ++n)
     {
     unsigned int nr = 1;
     for (int l = n; l >= 0; l -= 2)
     {
     unsigned int l_mom = l;
     
     HOPot pot (Parameters::mn, l_mom);
     Schroddy Sfunc (pot, H);
     GenericEigenvalues GenEig(Sfunc, nr, l_mom);
     double eig = GenEig.eigenvalue();
     double eigfun= Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(l_mom),
     psiPrime0(l_mom), eig , psiArray);
     
     std::cout << n << "\t" << nr << "\t" << l << "\t" << eig << "\t" << eigfun << std::endl;
     
     ++nr;
     }
     }*/
    
    /*========================================================================================
     * TEST to find the correct eigenvalues range for shooting method
     * Set two eigenvalues E, Emax and see if there is a change of sign in eigenfunctions
     *======================================================================================*/
    /*double Emax=400, E=250, pass=0.001;
     int N=(Emax-E)/pass, l_mom=1;
     std::vector<double> arrayeval;
     std::vector<double> arrayefun;
     HOPot pot (Parameters::mn, Parameters::hbar_omega, l_mom);
     std::ofstream file("test_range.txt");
     
     for (int i=0; i<=N; ++i)
     {
    	Schroddy s(pot);
    	arrayeval.push_back(E);
    	double eigfun=s.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, E, N_step);
    	arrayefun.push_back(eigfun);
    	E+=pass;
     }
     
     std::vector <double>::iterator walk1 = arrayeval.begin();
     std::vector <double>::iterator walk2 = arrayefun.begin();
     while (walk1 != arrayeval.end() && walk2 != arrayefun.end())
     {
    	file << *walk1 << "\t\t" << *walk2 << std::endl;
    	std::cout << *walk1 << "\t\t" << *walk2 << std::endl;
    	walk1++;
    	walk2++;
     }
     
     file.close();*/
    
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
