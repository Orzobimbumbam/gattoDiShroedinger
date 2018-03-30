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


int main(int argc, const char * argv[])
{
    
    const double H = 0.01;
    std::string inputPath = "Inputs/";
    std::string outputPath = "Outputs/";
        
    //Load the matrix of quantum numbers for each state from "orbitals.txt"
    std::vector <std::vector <unsigned int>> orbitals (36);
    std::fstream in ("../orbitals.txt", std::ios::in);
    readMatrix(orbitals, in, false);
    in.close();
    
    in.open(inputPath + "ca48.txt");
    std::vector<std::vector<double>> qrParam(12);
    readMatrix(qrParam, in, false);
    
    Element ca48(orbitals);
    
    const HOPot pot(Parameters::mn); //default is ground state
    const Schroddy sh(pot, H);
    ElementEigenfunctions elEigf = ca48.orbitalEigenfunction(sh, orbitals);
    NuclearDensity NDens;
    NDens.theoreticalDensity(elEigf, ca48.getLevelDegeneration());
    NDens.sogDensity(qrParam, H);
    
    std::ofstream fOut(outputPath + "refDensity.txt");
    fOut << NDens.getTheoreticalDensity();
    
    fOut.close();
    fOut.open(outputPath + "refSogDensity.txt");
    fOut << NDens.getSOGDensity();
    
    fOut.close();
    fOut.open(outputPath + "refFirstKSPotential.txt");
    KohnShamInverse ksi(pot, H);
    ksi.KSinverse(NDens, ksi);
    fOut << ksi.getKSPot();
    
    
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
