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

/*============================================
 * Convergence condition
 *==========================================*/

bool Theoreticaldensity::convergence (const std::vector<double>& empidensity, const std::vector<double>& thdensity) const
{
	double maxDiff = std::abs(thdensity[0]-empidensity[0]);
	unsigned int maxIndex = 0;
	for (int i = 0; i < empidensity.size(); ++i)
	{
		if(std::abs(thdensity[i]-empidensity[i]) > maxDiff)
		{
			maxDiff = std::abs(thdensity[i]-empidensity[i]);
			maxIndex = i;
		}
	}

	const double epsilon = empidensity[maxIndex]*0.01;
	return maxDiff < epsilon ? true : false; // convergence condition
}


/*===========================================
 * SOG density
 *=========================================*/

SOGdensity::SOGdensity() {}

void SOGdensity::sogDensity (const std::vector<std::vector<double>>& QRparameters, std::vector<double>& sogdensity, double h) const
{
    sogdensity.clear();
    std::vector<double> notNormal;
    notNormal.clear();
    double scalar = 0;
	const double alpha = sqrt(2./3.)*Parameters::rp;
	const double gamma = sqrt(2./3.)*Parameters::rms;
	const double beta = sqrt((gamma*gamma)-(alpha*alpha));

	const unsigned long NSteps = std::abs(Parameters::x_fin - Parameters::x_in)/h;
	double radiusx = Parameters::x_in;
	for (unsigned int r = 0 ; r < NSteps + 1; ++r)
	{
		const double c1 = 1/(2*pow(Parameters::PI,(3./2.))*radiusx);
		double c3 = 0;
		for (unsigned int i = 0; i < 11; ++i)
		{
			const double Ai = (Parameters::NN*Parameters::qe*QRparameters[i][1])/(1 + 2*(QRparameters[i][0]*QRparameters[i][0])/(gamma*gamma));
			const double exp1 = exp((-1)*((radiusx - QRparameters[i][0])/beta)*((radiusx - QRparameters[i][0])/beta));
			const double exp2 = exp((-1)*((radiusx + QRparameters[i][0])/beta)*((radiusx + QRparameters[i][0])/beta));
			const double c2 = exp1*(((radiusx + QRparameters[i][0])/(beta*beta*beta)) - (QRparameters[i][0]/(beta*gamma*gamma))) +
					exp2*(((radiusx - QRparameters[i][0])/(beta*beta*beta)) - (QRparameters[i][0]/(beta*gamma*gamma)));

			c3 += Ai*c2;
		}
		const double sogdens = c1*c3;
		notNormal.push_back(sogdens);
		scalar += sogdens*sogdens;
		//sogdensity.push_back(sogdens);

		radiusx += h;
	}
	// normalization
	double squared = sqrt(scalar*h);
    for (auto& it : notNormal)
    {
        it = it/squared;
        sogdensity.push_back(it);
    }
	return;
}

