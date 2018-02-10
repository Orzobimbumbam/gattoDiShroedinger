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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


int main(int argc, const char * argv[]) {
    
    //double eigenvalue1=0.93250 ;
    //double eigenvalue2=1.36256 ;
    //double error=10e-8;
    unsigned long N_step = 100;
    
    
    //const double omega = 2*Parameters::PI*Parameters::f;
    
    //for ( int i=0; i<=Parameters::angularMomentum; ++i)
    //{
    int l_mom=1;
    HOPot pot (Parameters::mn, l_mom);
    GenericEigenvalues GenEig(pot);
    double eig = GenEig.eigenvalue();
    Schroddy Sfunc (pot);
    
    //double eigfun= Sfunc.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, eig, N_step);
    
    //std::cout << /*"l value: "<< i <<*/ "Bisected Eigenvalue: " << eig << "Eigenfunction: " << eigfun << std::endl;
    //std::cout << "l value: "<< i << "Bisected Eigenvalue: " << eig << std::endl;
    
    //}
    double Emax=100, E=-100, pass=0.1;
    int N=(Emax-E)/pass;//, l_mom=1;
     //double arrayeval[N], arrayefun[N];
    std::vector<double>arrayeval(N), arrayefun(N);
     //HOPot pot (Parameters::mn, Parameters::hbar_omega, l_mom);
    
     //std::ofstream file("test_range.txt");
     //double* arraytest= new double [N];
     Schroddy s(pot);
     for (int i=0; i<N; ++i)
     {
    	
    	arrayeval[i]=E;
    	arrayefun[i]=s.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, E, N_step);
    	//file << arrayeval[i] << "      " << arrayefun[i] << std::endl;
    	//std::cout << arrayeval[i] << "      " << arrayefun[i] << std::endl;
    	E+=pass;
     }
     
     //file.close();
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
