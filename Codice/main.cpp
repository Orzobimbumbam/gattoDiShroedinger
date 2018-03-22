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
    in.close();
    
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
     * Calculate empirical densities from MC simulations or scattering (SOG densities)
     *=======================================================================================*/
    
    //std::fstream in2 ("Ca48.txt", std::ios::in);
    in.open("../Ca48.txt");
    std::ofstream file3("Outputs/SOGdensity.txt");
    std::vector <std::vector <double>> QRparam(12);
    std::vector<double> empidensity;
    
    for (int i = 0; i < 12; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            double temp;
            in >> temp;
            QRparam[i].push_back(temp);
        }
    }
    
    in.close();
    
    /*for (int i = 0; i < 12; ++i)
     {
     for (int j=0; j<2; ++j)
     std::cout << QRparam[i][j];
     std::cout << std::endl;
     }*/
    
    SOGdensity sogdensy;
    sogdensy.sogDensity(QRparam, empidensity, H);
    
    rad = Parameters::x_in;
    for (int i = 0; i < empidensity.size(); ++i)
    {
        file3 << rad << "\t" << empidensity[i] << std::endl;
        rad += H;
    }
    
    file3.close();
    
    /*=========================================================================================
     * STEP-3
     * Calculate new potential by inverse Kohn-Sham equations
     *=======================================================================================*/
    
    std::ofstream file4("Outputs/newpotential.txt");
    std::vector<double> newpotArray;
    std::vector<double> inpotArray;
    HOPot inpotential (Parameters::mn, 0);
    
    double radius = Parameters::x_in;
    for (int i = 0; i < NSteps + 1; ++i)
    {
        double inpot = inpotential.potential(radius);
        inpotArray.push_back(inpot);
        radius += H;
    }
    
    KohnShamInverse inversion;
    inversion.KSinverse(thdensArray,empidensity,inpotArray,newpotArray);
    
    rad = Parameters::x_in;
    for (int i = 0; i < newpotArray.size(); ++i)
    {
        file4 << rad << "\t" << newpotArray[i] << std::endl;
        rad += H;
    }
    
    file4.close();
    

    
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
