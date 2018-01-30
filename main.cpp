# include <iostream>
//# include "parameters.h"
# include "Includes.h"
//# include "schroddy.h"



int main(int argc, const char * argv[]) {


	double eigenvalue1=0.93250 ;
	double eigenvalue2=1.36256 ;

    //const unsigned int energyLevel = 1;
    //const int angularMomentum = 1;

    //const double f = 100;

    const double omega = 2*Parameters::PI*Parameters::f;
    HOPot pot (Parameters::mn, omega ) ;
    GenericEigenvalues eig(pot, eigenvalue1, eigenvalue2);
    //GenericEigenvalues eig(pot, Parameters::eigenvalue1, Parameters::eigenvalue2);
    //HarmonicEigenvalues eig(omega, Parameters::energyLevel, Parameters::angularMomentum);



    //WSaxPot pot(Parameters::V0, Parameters::Rn, Parameters::a0); //this is how you call the constructors; if you need pointers, call new

    Schroddy s(pot, eig); //your Schoddy object ready to go
    //double result=0;

	double result=s.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, Parameters::N_step);

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
    std::cout << result << std::endl;
    std::cout << "Program executed successfully!" << std::endl;
    return 0;
}








