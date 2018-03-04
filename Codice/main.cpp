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
    
    unsigned int k = 1; //same formula as yours, just keep this even
    unsigned int l_mom=0;
    unsigned int n = 2*k + l_mom;
    
    HOPot pot (Parameters::mn, l_mom);
    Schroddy Sfunc (pot, H);
    GenericEigenvalues GenEig(Sfunc, n, l_mom);
    double eig = GenEig.eigenvalue();
    //double eigenfun = Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, 14.5, vec );
    //double eigfun= Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, eig, );
    
    //std::cout << eigenfun << std::endl;
    std::cout << eig << std::endl;
    
    
    
    
    
    
    
    
    /*=========================================================================================
     * STEP-2
     * Calculate the theoretical density
     *=======================================================================================*/
    
    
    
    
    
    
    
    
    
    /*======================================================================================
     * FANCULO IL JSON!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     *====================================================================================*/
    
    /*Json::Value root;
     Json::Reader reader;
     std::ifstream ifs ("orbitals3.json");    // this read
     ifs >> root;							 // the entire file
     //Json::Value obj;
     //reader.parse(ifs, obj);
     
     const Json::Value array = root["orbitals"]["levels"];
     //const Json::Value& orbital = obj["orbital"];
     //const Json::Value& properties = obj["properties"];
     
     //std::cout << "orbital " << obj["orbital"].asString() << std::endl;
     
     
     for (int i=0; i<array.size(); ++i)
     {
    	std::cout << array[i].asInt() << std::endl;
     }
     for (int i=0; i<orbital.size(); ++i)
     {
    	std::cout << "orbitals " << orbital[i]["orbital"].asString() << std::endl;
     }
     */
    //const double omega = 2*Parameters::PI*Parameters::f;
    /*
     for ( int i=0; i<=Parameters::angularMomentum; ++i)
     {
     int l_mom=i;
     HOPot pot (Parameters::mn, Parameters::hbar_omega, l_mom);
     GenericEigenvalues GenEig(pot);
     double eig = GenEig.eigenvalue();
     Schroddy Sfunc (pot);
     double eigfun= Sfunc.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, eig, N_step);
     
     std::cout << "l value: "<< i << "\t"<< "Bisected Eigenvalue: " << eig << "\t" << "Eigenfunction: " << eigfun << std::endl;
     //std::cout << "l value: "<< i << "Bisected Eigenvalue: " << eig << std::endl;
     
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



