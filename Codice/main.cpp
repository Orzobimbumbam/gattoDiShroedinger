//
//  main.cpp
//  Codice
//
//  Created by Alberto Campi on 17/01/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#include <iostream>
#include "Includes.h"


int main(int argc, const char * argv[]) {
    
    const unsigned int energyLevel = 1;
    const int angularMomentum = 1;
    
    const double f = 100;

    const double omega = 2*Parameters::PI*f;
    HarmonicEigenvalues eig(omega, energyLevel, angularMomentum);
    WSaxPot pot(Parameters::V0, Parameters::Rn, Parameters::a0); //this is how you call the constructors; if you need pointers, call new
    
    //Schroddy s(pot, eig); //your Schoddy object ready to go
    
    // do stuff here...
    /*
    std::fstream file;
    int i, n, l;
    
    file.open("eigenfunc.dat",std::ios::out);
    
    for (i=0; i<N; i++)
    {
        for (l=0; l<nl; l++)
        {
            HBOEigen (omega, l, hbar, n);
            rungekutta(schrodinger, r, h, 4, u);
            
            
            file << hboeigen << "\t \t"	<< u[0] << "\t" << endl;
        }
        
    }
     
     
    
    file.close();
    file.clear();
    */
    
    std::cout << "Program executed successfully. Picio!" << std::endl;
    return 0;
}
