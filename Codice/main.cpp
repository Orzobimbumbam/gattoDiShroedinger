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
/*#include <jsoncpp/json/json.h>
 #include <jsoncpp/json/reader.h>
 #include <jsoncpp/json/writer.h>
 #include <jsoncpp/json/value.h>*/
#include "Includes.h"





int main(int argc, const char * argv[]) {
    
    double H = 0.001;
    int mass_num=Parameters::A;
    
    /*==========================================================================================
     * STEP-1
     * Calculate, by Schrodinger solver, the eigenfunctions for each quantum involved state
     *========================================================================================*/
    //Load the matrix with quantum numbers of each state
    //std::vector <std::vector <int> > orbitals (10);// std::vector <int> (3,0));
    //std::vector <int> state (3);
    //std::vector <std::vector <int> > orbitals;
    /*std::fstream in ("orbitals.txt", std::ios::in);
     
     for (int i=0; i<10; ++i)
     {
     for (int j=0; j<3; ++j)
     {
     int temp;
     in >> temp;
     orbitals[i].push_back(temp);
     }
     }
     
     /*for (int i=0; i<10; ++i)
     {
     for (int j=0; j<3; ++j)
     std::cout << orbitals[i][j];
     std::cout << std::endl;
     }*/
    
    // Solve Schrodinger equation for each involved state
    /*int i=0;
     int degen = 0;
     std::vector<double> arrayeval;
     std::vector<double> arrayefun;
     
     while (mass_num > 0)
     {
    	int quantN=orbitals[i][0];
     int quantL=orbitals[i][1];
     HOPot pot (Parameters::mn, Parameters::hbar_omega, quantL);
     Schroddy Sfunc (pot, H);
     GenericEigenvalues GenEig(pot, quantN);
     double eig = GenEig.eigenvalue();
     double eigfun= Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, eig);
     
     arrayeval.push_back(eig);
     arrayefun.push_back(eigfun);
     
     degen = 2*(2*quantL+1);
     
     std::cout << "l value: "<< quantL << "\t"<< "Bisected Eigenvalue: " << eig << "\t" << "Eigenfunction: " << eigfun << std::endl;
     
     mass_num -= degen;
     i++;
     std::cout << mass_num << std::endl;
     }*/
    
    unsigned int k = 3;
    unsigned int l_mom = 0;
    unsigned int n = 2*(k-1) + l_mom; //2*(k-1) + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, k, l_mom); //we MUST look for the radial nodes only as we solve the radial part of S.
    double eig = GenEig.eigenvalue();
    std::vector<double> psiArray;
    double s = Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(l_mom),
                                       psiPrime0(l_mom), eig , psiArray);
    
    
    //std::cout << eigenfun << std::endl;
    std::cout << eig << std::endl;
    
    /*
    std::ofstream fOut;
    fOut.open("Outputs/eigenfunctions.txt", std::ios_base::out | std::ios_base::app);
    for (unsigned int i = 1; i <= 3; ++i)
    {
        for (unsigned int j = 0; j <= i; ++j)
        {
            HOPot pot_ (Parameters::mn, j);
            Schroddy Sfunc_ (pot_, H);
            GenericEigenvalues GenEig_(Sfunc_, i, j);
            double eig_ = GenEig_.eigenvalue();
            std::vector<double> psiArray_;
            double s_ = Sfunc_.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(j),
                                               psiPrime0(j), eig_ , psiArray_);
            
            for (const auto& v : psiArray_)
            {
                
                fOut << v << ";";
                
            }
            
            fOut << std::endl;
        }
    }
    
    
    
    
    /*=========================================================================================
     * STEP-2
     * Calculate the theoretical density
     *=======================================================================================*/
    
    
    

    std::cout << "Program executed successfully." << std::endl;
    return 0;
}



