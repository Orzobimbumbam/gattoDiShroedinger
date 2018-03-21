#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
#include "includes.h"


int main(int argc, const char * argv[]) {

	double H = 0.001;
	int massNumb = Parameters::A;
	clock_t start = clock(); // start time
	/*time_t t1, t2;
	t1 = time(0); //start time*/

/*==========================================================================================
 * STEP-1
 * Calculate, by Schrodinger solver, the eigenfunctions for each quantum involved state
 *========================================================================================*/

	// Load the matrix with quantum numbers of each state from "orbitals.txt"
	std::vector <std::vector <int> > orbitals (36);
	std::vector<double> thdensity;
	std::ifstream in ;
	std::string inputfile1 = Parameters::orbitalsfile;
	in.open(inputfile1);

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
	in.clear();

		/*for (int i=0; i<36; ++i)
		{
			for (int j=0; j<2; ++j)
				std::cout << orbitals[i][j];
				std::cout << std::endl;
		}*/

	// Fill orbitals and calculate energies, eigenfunctions and theoretical density for first loop
	OrbitalsFilling fill (massNumb, H);
	fill.orbFilling(orbitals, thdensity);

	std::ofstream file2("density.txt");
	double rad = Parameters::x_in;
    for (int i = 0; i < thdensity.size(); ++i)
    {
    	file2 << rad << "\t" << thdensity[i] << std::endl;
    	rad += H;
    }
    file2.close();

/*=========================================================================================
 * STEP-2
 * Calculate empirical densities from MC simulations or scattering (SOG densities)
 *=======================================================================================*/
    // Load Qi and Ri parameters for SOG
    //std::fstream in2 ("Ca48.txt", std::ios::in);
    std::string inputfile2 = Parameters::qrfile;
    in.open(inputfile2);
    std::vector <std::vector <double> > QRparam (12);
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

 	// Claculate empiric density by SOG
	SOGdensity sogdensy;
	sogdensy.sogDensity(QRparam, empidensity, H);

    std::ofstream file3("SOGdensity.txt");
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

	const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/H;
	std::vector<double> inpotArray;
	std::vector<double> psiArray;
    Theoreticaldensity densy;
    KohnShamInverse inversion;

    int count = 0; //Loop counter
    while(!densy.convergence(empidensity, thdensity))
    {
    	if (count == 0)
    	{
    		HOPot inpotential (Parameters::mn, 0);

     		double radius = Parameters::x_in;
     		for (int i = 0; i < NSteps + 1; ++i)
     		{
     			double inpot = inpotential.potential(radius);
     			inpotArray.push_back(inpot);
     			radius += H;
     		}
    	}
    	else
    	{
    		psiArray.clear();
    		thdensity.clear();
    		potOut newpot (inversion);
    		Schroddy Sfunc (newpot, H);
    		GenericEigenvalues GenEig(Sfunc, quantNr, quantL);
    		double eigval = GenEig.eigenvalue();
    		Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(quantL),
    	             psiPrime0(quantL), eigval , psiArray)
    		densy.density(psiArray,thdensity,nuclNumb,H);
    	}

 		inversion.KSinverse(thdensity,empidensity,inpotArray);

 		count++;
 		std::cout << "loop n. " << count << std::endl;
     }

    //t2 = time(0);
    clock_t end = clock(); // finish time
    //double hours = (((double)(end - start))/CLOCKS_PER_SEC)/3600;
    double hours = ((double)(end - start))/3600;

    //std::cout << "CONVERGENCE IS DONE in: " << t2 - t1 << "sec.!" << "GREAT JOB!" << std::endl;
    std::cout << "CONVERGENCE IS DONE in: " << hours << "hours!" << "GREAT JOB!" << std::endl;

    std::ofstream file4("newpotential.txt");
    rad = Parameters::x_in;
    for (int i = 0; i < newpotArray.size(); ++i)
    {
 	   	 file4 << rad << "\t" << newpotArray[i] << std::endl;
 	   	 rad += H;
    }
    file4.close();








 /*===============================================================================
  * TEST - Calculate first 12 levels eigenvalues for HO potential
  *==============================================================================*/

     /*unsigned int nr=1;
     int l_mom=2;
     unsigned int n = 2*(nr-1) + l_mom;*/
 	/*std::vector <double> psiArray;
     unsigned int nmax = 5;
     for (unsigned int n = 0; n <= nmax; ++n)
     {
     	unsigned int nr = 1;
     	for (int l = n; l >= 0; l -= 2)
     	{
     		unsigned int l_mom = l;
     		HOPot pot (Parameters::mn, l_mom);
     		Schroddy Sfunc (pot, H);
     		GenericEigenvalues GenEig(Sfunc, nr, l_mom);
     		double eig = GenEig.eigenvalue();
     		double eigfun= Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(l_mom),
                     psiPrime0(l_mom), eig , psiArray);
     		std::cout << n << "\t" << nr << "\t" << l << "\t" << eig << "\t" << eigfun << std::endl;
     		++nr;
     	}
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
