#include <iostream>
#include <fstream>
#include <vector>
#include <jsoncpp/json/json.h>
#include "Includes.h"


int main(int argc, const char * argv[]) {
    
    //double eigenvalue1=0.93250 ;
    //double eigenvalue2=1.36256 ;
    //double error=10e-8;


	unsigned long N_step = 1000;
	std::ifstream ifs ("orbitals2.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);

    const Json::Value& orbital = obj["orbital"];
    //const Json::Value& properties = obj["properties"];

    //std::cout << "orbital " << obj["orbital"].asString() << std::endl;

    for (int i=0; i<orbital.size(); ++i)
    {
    	std::cout << "orbitals " << orbital[i]["orbital"].asString() << std::endl;
    }

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
    /*double Emax=10, E=0, pass=0.1;
    int N=(Emax-E)/pass, l_mom=1;
    std::vector<double> arrayeval (N);
    std::vector<double> arrayefun (N);
    HOPot pot (Parameters::mn, Parameters::hbar_omega, l_mom);
    std::ofstream file("test_range.txt");
    //double* arraytest= new double [N];

    for (int i=0; i<=N; ++i)
    {
    	Schroddy s(pot);
    	arrayeval.push_back(E);
    	double eigfun=s.solveShroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0, E, N_step);
    	arrayefun.push_back(eigfun);
    	//file << arrayeval << "      " << arrayefun << std::endl;
    	//std::cout << arrayeval [i] << "      " << arrayefun [i] << std::endl;
    	E+=pass;
    }

    for (int i=0; i<arrayeval.size(); ++i)
    {
    	for (int j=0; j<arrayefun.size(); ++j)
    	{
    		std::cout << arrayeval [i] << "      " << arrayefun [j] << std::endl;
    	}
    }

    file.close();*/

    std::cout << "Program executed successfully." << std::endl;
    return 0;
}
