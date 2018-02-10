#include <iostream>
#include <fstream>
#include <vector>
/*#include <jsoncpp/json/json.h>
#include <jsoncpp/json/reader.h>
#include <jsoncpp/json/writer.h>
#include <jsoncpp/json/value.h>*/
#include "Includes.h"


int main(int argc, const char * argv[]) {

	unsigned long N_step = 1000;
	int mass_num=Parameters::A;

/*==========================================================================================
 * STEP-1
 * Calculate eigenfunctions for each quantum involved state by Schrodinger solver
 *========================================================================================*/
	//Load the matrix with quantum and degeneration numbers of each state
	//std::vector <std::vector <int> > orbitals (10, std::vector <int> (3,0));
	std::vector <int> state (3);
	std::vector <std::vector <int> > orbitals;
	std::fstream in ("orbitals.txt", std::ios::in);

	for (int i=0; i<10; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			in >> state [j];
			orbitals.push_back(state);
		}
	}

		/*for (int i=0; i<10; ++i)
		{
			for (int j=0; j<3; ++j)
				std::cout << orbitals[i][j] << std::endl;
		}*/

	// Solve Schrodinger equation for each involved state
	int i=0;
    std::vector<double> arrayeval;
    std::vector<double> arrayefun;

    while (mass_num > 0)
	{
		int l_mom=orbitals[i][1];
		HOPot pot (Parameters::mn, Parameters::hbar_omega, l_mom);
		GenericEigenvalues GenEig(pot);
		double eig = GenEig.eigenvalue();
		Schroddy Sfunc (pot);
		double eigfun= Sfunc.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, eig, N_step);

		arrayeval.push_back(eig);
		arrayefun.push_back(eigfun);

		std::cout << "l value: "<< l_mom << "\t"<< "Bisected Eigenvalue: " << eig << "\t" << "Eigenfunction: " << eigfun << std::endl;

		i++;
		mass_num -= orbitals[i][2];
		std::cout << mass_num << std::endl;
	}









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

/*========================================================================================
 * TEST to find the correct eigenvalues range for shooting method
 * Set two eigenvalues E, Emax and see if there is a change of sign in eigenfunctions
 *======================================================================================*/
    /*double Emax=400, E=250, pass=0.001;
    int N=(Emax-E)/pass, l_mom=1;
    std::vector<double> arrayeval (N);
    std::vector<double> arrayefun (N);
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
