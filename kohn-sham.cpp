# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include "kohn-sham.h"
# include <cmath>


/*============================================
 * Theoretical density
 *==========================================*/
Densities::~Densities() {}

Theoreticaldensity::Theoreticaldensity(double eigenfunc, double x): m_eigenfunc(eigenfunc), m_x(x) {}

double Theoreticaldensity::density() const
{
	//m_eigenfunc& schroddy::normalPsi; //this line doesn't make sense..
    return (1/(4*Parameters::PI*(m_x*m_x)))*Parameters::NN*(m_eigenfunc*m_eigenfunc);// BOH!!!!!!!
}

Densities* Theoreticaldensity::clone() const
{
    return new Theoreticaldensity(*this);
}

/*=============================================
 * Kohm-Sham inverse equations
 *===========================================*/





