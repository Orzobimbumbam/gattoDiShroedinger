# include <iostream>
# include <fstream>
# include <vector>
# include <cmath>
# include "includes.h"

OrbitalsFilling::OrbitalsFilling(int massNumber, const double h): m_massNumber(massNumber), m_h(h) {}

void OrbitalsFilling::orbFilling(std::vector <std::vector <int> >& orbiMatrix, std::vector<double>& thdensArray)
{
	// Solve Schrodinger equation for each involved state
	int i = 0, degen = 0;
	unsigned int nuclNumb = 0;
	const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/m_h;
	std::vector <double> psiArray;
	std::ofstream file1("eigenfunc.txt");

	while (m_massNumber > 0)
	{
	    if(m_massNumber == Parameters::A)	// initialize a vector of zero elements
	    {									// for starting sum in density
	    	for (int i = 0; i < NSteps + 1; ++i)
	    	{
	    		thdensArray.push_back(0);
	    	}
	    }

	    unsigned int quantNr = orbiMatrix[i][0];
		unsigned int quantL = orbiMatrix[i][1];
		unsigned int quantN = 2*(quantNr - 1) + quantL;

		HOPot pot (Parameters::mn, quantL);
		Schroddy Sfunc (pot, m_h);
		GenericEigenvalues GenEig(Sfunc, quantNr, quantL);
		double eigval = GenEig.eigenvalue();
		Sfunc.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, psi0(quantL),
	             psiPrime0(quantL), eigval , psiArray);

		degen = 2*(2*quantL+1); // orbital degeneration
		m_massNumber -= degen; // remaining nucleons

		if (m_massNumber < 0)
			nuclNumb = m_massNumber + degen; // nucleons number on last orbital in not filled case
		else nuclNumb = degen; // nucleons number in filled orbital

		// Claculate theoretical density for involved state
		Theoreticaldensity densy;
		densy.density(psiArray,thdensArray,nuclNumb,m_h);


		std::vector<double>::iterator walk1 = psiArray.begin();
		while (walk1 != psiArray.end())
		{
		    file1 << quantN << "\t" << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << eigval << "\t\t\t" << *walk1 << std::endl;
		    walk1++;
		}

		std::cout << quantNr << "\t" << quantL << "\t" << nuclNumb << "\t" << eigval << "\t" << std::endl;

		i++;

	}
	file1.close();
	return;

}



