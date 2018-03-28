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
#include <string>

#include "Includes.h"


typedef std::vector<double>::const_iterator cwalk;
typedef std::vector<double>::iterator walk;

int main(int argc, const char * argv[])
{
    
    const double H = 0.01;
    std::string inputPath = "Inputs/";
    std::string outputPath = "Outputs/";
    //int massNumb=Parameters::A;
    
    /*==========================================================================================
    * STEP-1
    * Calculate, by Schrodinger solver, the eigenfunctions for each quantum involved state
    *========================================================================================*/
        
    //Load the matrix with quantum numbers of each state from "orbitals.txt"
    std::vector <std::vector <unsigned int>> orbitals (36);
    std::fstream in ("../orbitals.txt", std::ios::in);
    readMatrix(orbitals, in, false);
    in.close();
    
    Element ca48(orbitals);
    
    const HOPot pot(Parameters::mn); //ground state
    const Schroddy sh(pot, H);
    ElementEigenfunctions elEigf = ca48.orbitalEigenfunction(sh, orbitals);
    Theoreticaldensity thDens;
    thDens.density(elEigf, ca48.getLevelDegeneration());
    
    std::ofstream fOut(outputPath + "refDensity.txt");
    fOut << thDens;
    
    
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
