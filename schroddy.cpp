//Class for Schrodinger equation solver

# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
//# include <cmath>

using namespace Parameters;

HarmonicEigenvalues::HarmonicEigenvalues(double omega): m_omega(omega) {}
double HarmonicEigenvalues::eigenvalue(unsigned int n, int l) const
{
    return hbar*m_omega*(2*n+l+(3/2));
}

double Schroddy::solveShroddyByRK(const InitialPot& pot, const Eigenvalues& eigenval, double x0, double x1, double step) const
{
    //implement runge-kutta here..
    
    return 0; //change return value
}
	
//=====================================================================
// Schrodinger function for Runge-Kutta method.
//=====================================================================
/*
double schrodinger (double& r, double& l, double& *Vks, double& u)
{	
	double eigenvalue=HBOEigen(omega, l, hbar, n);
	double schrod=-2*(eigenvalue-((l*(l+1))/(2*(r*r)))-*Vks)*u;
	
	return schrod;
}	

	/*
int main() {
	
		//const double hbar=1.0545718 * pow(10,-34);
		*/
	
	
	
//return 0;
//}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	}

