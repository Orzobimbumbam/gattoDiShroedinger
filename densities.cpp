# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include "kohn-sham.h"
# include <cmath>


/*============================================
 * Theoretical density
 *==========================================*/
Densities::~Densities() {}

Theoreticaldensity::Theoreticaldensity(double eigenfunc, int degen, double x): m_eigenfunc(eigenfunc), m_degen(degen), m_x(x) {}

double Theoreticaldensity::density() const
{

    return (1/(4*Parameters::PI*(m_x*m_x)))*m_degen*(m_eigenfunc*m_eigenfunc);
}

Densities* Theoreticaldensity::clone() const
{
    return new Theoreticaldensity(*this);
}



