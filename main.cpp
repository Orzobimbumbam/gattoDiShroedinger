#include <iostream>
#include <fstream>
#include <vector>
#include <string>
/*#include <jsoncpp/json/json.h>
#include <jsoncpp/json/reader.h>
#include <jsoncpp/json/writer.h>
#include <jsoncpp/json/value.h>*/
#include "Includes.h"


int main(int argc, const char * argv[]) {

	double H = 0.001;
	int massNumb=Parameters::A;

/*==========================================================================================
 * STEP-1
 * Calculate, by Schrodinger solver, the eigenfunctions for each quantum involved state
 *========================================================================================*/

	//Load the matrix with quantum numbers of each state from "orbitals.txt"
	std::vector <std::vector <int> > orbitals (36);// std::vector <int> (3,0));
	//std::vector <int> state (2);
	//std::vector <std::vector <int> > orbitals;
	std::fstream in ("orbitals.txt", std::ios::in);

	for (int i = 0; i < 36; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			int temp;
			in >> temp;
			orbitals[i].push_back(temp);
		}
	}

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
	std::vector<double> thdensArray;
    /*std::vector<double> evalArray;
    std::vector<double> efuncArray;
    std::vector<int> nuclnumbArray;*/
    std::ofstream file("eigenfunc.txt");



     while (massNumb > 0)
	{
    	if(massNumb == Parameters::A)
    	{
    		for (int i = 0; i < NSteps + 1; ++i)
    		{
    			thdensArray.push_back(0);
    		}
    	}

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
	    	file << quantN << "\t" << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << *walk1 << std::endl;
	    	walk1++;
	    }

		/*evalArray.push_back(eigval);
		efuncArray.push_back(eigfunc);
		nuclnumbArray.push_back(nuclNumb);*/

		//dstd::cout << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << eigval << "\t" << eigfunc << std::endl;

		i++;
		//std::cout << massNumb << std::endl;
	}
     double rad = Parameters::x_in;
     std::ofstream file2("density.txt");
     std::vector<double>::iterator walk2 = thdensArray.begin();
     //std::vector<double>::iterator walk3 = xArray.begin();
     while (walk2 != thdensArray.end() /*&& walk3 != xArray.end()*/)
     {
         file2 << rad << "\t"<< *walk2 << std::endl;
         walk2++, rad + H;
     }

     file2.close();

     file.close();
     //in.clear();

/*=========================================================================================
 * STEP-2
 * Calculate the theoretical density
 *=======================================================================================*/

 	// std::vector <std::vector <double> > efunctions (63500);
 	 //std::vector <int> data (4);
 	 //std::fstream in ("eigenfunc.txt", std::ios::in);
 	 //in.open ("eigenfunc.txt");
 	/*std::fstream in ("eigenfunc.txt", std::ios::in);

 	 for (int i = 0; i < 63500; ++i)
 	 {
 		 for (int j = 0; j < 2; ++j)
 		 {
 			 double temp;
 			 in >> temp;
 			 efunctions[i].push_back(temp);
 		 }
 		 //if (in.eof()) break;
 	 }*/

 	/*for (int i=0; i<63500; ++i)
 			{
 				for (int j=0; j<2; ++j)
 					std::cout << efunctions[i][j];
 					std::cout << std::endl;
 			}*/


    /*std::vector<double> thDensArray;
     std::vector<double> xArray;
     std::ofstream file2("density.csv");
 	const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/H;

	int deg = 0;
	double efunc = 0, radiusx = Parameters::x_in;
	for (unsigned long i = 0; i < NSteps + 1; ++i)
	{
		deg = efunctions[i][0];
		efunc = efunctions[i][1];
		Theoreticaldensity densy;
		double thdensity = densy.density(efunc,deg,radiusx);
		thDensArray.push_back(thdensity);
		xArray.push_back(radiusx);
		radiusx += H;
	}

     std::vector<double>::iterator walk2 = thDensArray.begin();
     std::vector<double>::iterator walk3 = xArray.begin();
     while (walk2 != thDensArray.end() && walk3 != xArray.end())
     {
     	file2 << *walk3 << ","<< *walk2 << std::endl;
     	walk2++, walk3++;
     }

     file2.close();*/

/*=========================================================================================
 * STEP-3
 * Calculate empirical densities from MC simulations or scattering
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
