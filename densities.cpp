# include "parameters.h"
# include "densities.h"
# include <fstream>
# include <vector>


/*============================================
 * Theoretical density
 *==========================================*/
//Theoreticaldensity::~Theoreticaldensity() {}

Theoreticaldensity::Theoreticaldensity() {}

double Theoreticaldensity::density(std::vector<double> psi, std::vector<double> thDensArray, int degen, double step) const
{
	double radiusx = Parameters::x_in;
	for ( int i = 0; i < psi.size(); ++i)
	{
		double thdensity = (1/(4*Parameters::PI*(radiusx*radiusx)))*degen*(psi[i]*psi[i]);
		thDensArray[i] += thdensity;
		//thDensArray.push_back(thdensity);
		radiusx += step;
	}
	return radiusx;
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
