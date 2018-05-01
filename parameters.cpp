#include "parameters.h"

double psi0(unsigned int l)
{
    return pow(Parameters::x_in, l + 1);
}

double psiPrime0(unsigned int l)
{
	/*if (l == 0) return 0.;
	else*/ return (l + 1)*pow(Parameters::x_in, l);
}
