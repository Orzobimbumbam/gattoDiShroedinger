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
    
    //double eigenvalue1=0.93250 ;
    //double eigenvalue2=1.36256 ;
    //double error=10e-8;
    //unsigned long N_step = 10;
    
    
    const double omega = 2*Parameters::PI*Parameters::f;
    HOPot pot (Parameters::mn, omega, 1) ;
    GenericEigenvalues GenEig(pot);
    double eig = GenEig.eigenvalue();
    std::cout << "Bisected Eigenvalue: " << eig << std::endl;
    
    
    
    
    std::cout << "Program executed successfully. Picio!" << std::endl;
    return 0;
}
