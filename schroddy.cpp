//Class for Schrodinger equation solver

# include <iostream>
# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include "functions.h"
# include "initpot.h"
# include "parameters.h"
# include <cmath>
# include <fstream>

using namespace std;

class schroddy
{

double Vks=HBOPot(m, omega, rhbo);
	
//=====================================================================
// HBO eigenvalues
//=====================================================================

double HBOEigen ( double& omega, double& l, double& hbar, int& n)
{
	double hboeigen=hbar*omega*(2*n+l+(3/2));
	
	return hboeigen;
		
}	
	
//=====================================================================
// Schrodinger function for Runge-Kutta method.
//=====================================================================
	
double schrodinger (double& r, double& l, double& *Vks, double& u)
{	
	double eigenvalue=HBOEigen(omega, l, hbar, n);
	double schrod=-2*(eigenvalue-((l*(l+1))/(2*(r*r)))-*Vks)*u;
	
	return schrod;
}	

	
int main() {
	
		//const double hbar=1.0545718 * pow(10,-34);
	fstream file;
	int i, n, l;
	
	file.open("eigenfunc.dat",ios::out);
	
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
	
	
	
return 0;
}	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	}

