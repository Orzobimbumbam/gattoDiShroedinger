#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Includes.h"


int main(int argc, const char * argv[]) {

	double H = 0.001;
	int massNumb=Parameters::A;

/*==========================================================================================
 * STEP-1
 * Calculate, by Schrodinger solver, the eigenfunctions for each quantum involved state
 *========================================================================================*/

	//Load the matrix with quantum numbers of each state from "orbitals.txt"
	std::vector <std::vector <int> > orbitals (36);
	std::ifstream in ;
	in.open("orbitals.txt");

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

	// Solve Schrodinger equation for each involved state
	int i = 0;
	int degen = 0;
	unsigned int nuclNumb = 0;
	const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/H;
	std::vector <double> psiArray;
	std::vector<double> thdensArray(NSteps+1);
	double rad = Parameters::x_in;
    /*std::vector<double> evalArray;
    std::vector<double> efuncArray;
    std::vector<int> nuclnumbArray;*/
    std::ofstream file1("eigenfunc.txt");
    std::ofstream file2("density.txt");



     while (massNumb > 0)
	{
    	unsigned int quantNr = orbitals[i][0];
		unsigned int quantL = orbitals[i][1];
		unsigned int quantN = 2*(quantNr - 1) + quantL;
		HOPot pot (Parameters::mn, quantL);
		Schroddy Sfunc (pot, H);
		GenericEigenvalues GenEig(Sfunc, quantNr, quantL);
		double eigval = GenEig.eigenvalue();
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


	    std::vector<double>::iterator walk1 = psiArray.begin();
	    while (walk1 != psiArray.end())

	    {
	    	file1 << quantN << "\t" << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << *walk1 << std::endl;
	    	walk1++;
	    }

		//dstd::cout << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << eigval << "\t" << eigfunc << std::endl;

		i++;
	}

     for (int i = 0; i < thdensArray.size(); ++i)
     {
    	 file2 << rad << "\t" << thdensArray[i] << std::endl;
    	 rad += H;

     }

     /*std::vector<double>::iterator walk2 = thdensArray.begin();
     while (walk2 != thdensArray.end())
     {
    	 file2 << rad << "\t"<< *walk2 << std::endl;
    	 walk2++, rad += H;
     }*/

     file1.close();
     file2.close();

/*=========================================================================================
 * STEP-2
 * Calculate empirical densities from MC simulations or scattering (SOG densities)
 *=======================================================================================*/

     //std::fstream in2 ("Ca48.txt", std::ios::in);
     in.open("Ca48.txt");
     std::ofstream file3("SOGdensity.txt");
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

    std::ofstream file4("newpotential.txt");
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
