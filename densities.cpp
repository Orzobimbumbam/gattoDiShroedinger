# include "parameters.h"
# include "densities.h"
# include "eigenfunction.hpp"
# include <fstream>
# include <vector>


/*============================================
 * Theoretical density
 *==========================================*/
//Theoreticaldensity::~Theoreticaldensity() {}

//Theoreticaldensity::Theoreticaldensity() {}

void Theoreticaldensity::density(const Eigenfunction& psi, unsigned int degen, double step)
{
    for (const auto& it : psi.keyValues())
	{
		const double thdensity = (1/(4*Parameters::PI*(it.first*it.first)))*degen*(it.second*it.second);
		m_thDensity[it.first] += thdensity;
	}
	return;
}

/*============================================
 * Convergence condition
 *==========================================*/

bool Theoreticaldensity::hasConverged (const std::map<double, double>& empidensity) const
{
	double maxDiff = std::abs(m_thDensity.begin() -> second - empidensity.begin() -> second);
	double xMax = 0;
    for (const auto& it : m_thDensity)
	{
		if(std::abs(it.second - empidensity.at(it.first)) > maxDiff) //access only, throw exception if not key is not found
		{
			maxDiff = std::abs(it.second - empidensity.at(it.first));
			xMax = it.first;
		}
	}

	const double epsilon = empidensity.at(xMax)*0.01;
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
		const double c1 = 1./(2*pow(Parameters::PI, (3./2.))*radiusx);

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
		//scalar += sogdens*sogdens;
		scalar += ((notNormal[r]*radiusx*radiusx)+(notNormal[r - 1]*(radiusx - h)*(radiusx - h)))*h/2.; // integral by trapezoid method for normalization
		//sogdensity.push_back(sogdens);

		radiusx += h;
	}

	// Normalization
	double norm = Parameters::NN/scalar/4/Parameters::PI;
	//double squared = sqrt(scalar*h);
    for (auto& it : notNormal)
    {
        it = it*norm;
        sogdensity.push_back(it);
    }
	return;
}

