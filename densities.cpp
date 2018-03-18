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
		const double c1 = 1./(2*pow(Parameters::PI,(3/2))*r);
		double c3 = 0;
		for (int i = 0; i < 11; ++i)
		{
			const double Ai = (Parameters::NN*QRparameters[1])/(1 + 2*(QRparameters[0]*QRparameters[0])/(gamma*gamma));
			const double exp1 = exp(-((r - QRparameters[0])/beta)*((r - QRparameters[0])/beta));
			const double exp2 = exp(-((r + QRparameters[0])/beta)*((r + QRparameters[0])/beta));
			const double c2 = exp1*((r + QRparameters[0])/(beta*beta*beta) - QRparameters[0]/(beta*gamma*gamma)) +
					exp2*((r - QRparameters[0])/(beta*beta*beta) - QRparameters[0]/(beta*gamma*gamma));

			c3 += Ai*c2;
		}

		const double sogdens = c1*c3;
		sogdensity.push_back(sogdens);
	}
	return;
}



















/*Densities::~Densities() {}
Theoreticaldensity::Theoreticaldensity(char filename): m_filename(filename), m_h(0.001) {}
double Theoreticaldensity::density(double x0, double x1, std::vector<double> thDensArray, std::vector<double> xArray) const
{
	//load eigenfunctions and nucleons number for level from file
	std::vector <std::vector <int> > efunctions;
	std::vector <int> data (4);
	std::fstream in ("m_filename", std::ios::in);
	for (int i = 0; i < in.eof(); ++i)
	{
		for (int j=0; j<4; ++j)
		{
			int temp;
			in >> temp;
			efunctions[i].push_back(temp);
		}
	}
	//calculate the theoretical density
	int deg = 0;
	double efunc = 0, thdensity = 0, radiusx = x0;
	const unsigned long NSteps = std::abs(x1 - x0)/m_h;
	for (int i = 0; i < efunctions.size(); ++i)
	{
		for (unsigned long r = 0; r < NSteps; ++r)
		{
			deg = efunctions[i][3];
			efunc = efunctions[i][4];
			thdensity = (1/(4*Parameters::PI*(radiusx*radiusx)))*deg*(efunc*efunc);
			thDensArray.push_back(thdensity);
			xArray.push_back(radiusx);
			radiusx += m_h;
		}
	}
    return thdensity;
}
std::unique_ptr<Densities> Theoreticaldensity::clone() const
{
    return std::make_unique<Theoreticaldensity> (*this); //return a derived class object through a base class pointer
}*/
