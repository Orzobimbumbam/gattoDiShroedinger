# include "parameters.h"
# include "densities.h"
# include <fstream>
# include <vector>


/*============================================
 * Theoretical density
 *==========================================*/
//Theoreticaldensity::~Theoreticaldensity() {}

Theoreticaldensity::Theoreticaldensity() {}

void Theoreticaldensity::density(const std::vector<double>& psi, std::vector<double>& thDensArray, unsigned int degen, double step) const
{
	double radiusx = Parameters::x_in;
	for ( int i = 0; i < psi.size(); ++i)
	{
		double thdensity = (1/(4*Parameters::PI*(radiusx*radiusx)))*degen*(psi[i]*psi[i]);
		thDensArray[i] += thdensity;
		//thDensArray.push_back(thdensity);
		radiusx += step;
	}
	return;
}

/*===========================================
 * SOG density
 *=========================================*/

SOGdensity::SOGdensity() {}

void SOGdensity::sogDensity (const std::vector<std::vector<double>>& QRparameters, std::vector<double>& sogdensity, double h) const
{
    sogdensity.clear();
	const double alpha = sqrt(3./2)*Parameters::rp;
	const double gamma = sqrt((3./2)*(Parameters::rms*Parameters::rms));
	const double beta = sqrt((gamma*gamma)-(alpha*alpha));

	const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/h;
	for (unsigned int r = 0 ; r < NSteps +1; ++r)
	{
		const double c1 = 1/(2*pow(Parameters::PI,(3./2))*r);
		double c3 = 0;
		for (unsigned int i = 0; i < 11; ++i)
		{
			const double Ai = (Parameters::NN*QRparameters[i][1])/(1 + 2*(QRparameters[i][0]*QRparameters[i][0])/(gamma*gamma));
			const double exp1 = exp(-((r - QRparameters[i][0])/beta)*((r - QRparameters[i][0])/beta));
			const double exp2 = exp(-((r + QRparameters[i][0])/beta)*((r + QRparameters[i][0])/beta));
			const double c2 = exp1*((r + QRparameters[i][0])/(beta*beta*beta) - QRparameters[i][0]/(beta*gamma*gamma)) +
					exp2*((r - QRparameters[i][0])/(beta*beta*beta) - QRparameters[i][0]/(beta*gamma*gamma));

			c3 += Ai*c2;
		}
		const double sogdens = c1*c3;
		sogdensity.push_back(sogdens);
	}
	return;
}

