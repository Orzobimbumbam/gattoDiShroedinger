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
#include <pthread.h>

#include "Includes.h"


int main(int argc, const char * argv[])
{
    
    const double H = 0.1;
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
    
    //Test theoretical and sog density for initial harmonic potential
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
    
    //Test Kohn-Sham inversion for initial harmonic potential
    fOut.close();
    fOut.open(outputPath + "refFirstKSPotential.txt");
    KohnShamInverse ksi(pot, H);
    KohnShamInverse tempKsi = ksi;
    ksi.KSinverse(NDens, tempKsi);
    fOut << ksi.getKSPot();
    fOut.close();
    
    unsigned long loops = 0;
    while (!NDens.hasConverged()) //simulation loop
    {
        const PotOut po(ksi, Parameters::mn);
        const Schroddy sh_(po, H);
        elEigf = ca48.orbitalEigenfunction(sh_, orbitals);
        NDens.theoreticalDensity(elEigf, ca48.getLevelDegeneration());
        KohnShamInverse tempKsi_ = ksi;
        ksi.KSinverse(NDens, tempKsi_);
        
        ++loops;
        //if (loops%10 == 0)
            std::cerr << "Convergence distance: " << NDens.distanceToConvergence()
            << " after " << loops << " iterations." << std::endl;
    }
    
    
    KSPotential finalPotential = ksi.getKSPot();
    fOut.open(outputPath + "refFinalPotential.txt");
    fOut << finalPotential;
    fOut.close();
    
    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
