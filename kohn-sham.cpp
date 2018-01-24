# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include "kohn-sham.h"
# include <cmath>


/*============================================
 * Theoretical density
 *==========================================*/
Theoreticaldensity::~Theoreticaldensity() {}

TheoNuclearDensity::TheoNuclearDensity(double normalPsi, double x): m_normalPsi(normalPsi) {} //this class cannot have a public constructor if virtual

double TheoNuclearDensity::theodensity(double x) const //MUST be defined in derived class
{
    return (1/(4*Parameters::PI*(x*x)))*Parameteres::NN*(m_normalPsi*m_normalPsi);// BOH!!!!!!!
}

Theoreticaldensity* TheoNuclearDensity::clone() const //same argument as before, if not virtual cloning isn't necessary
{
    return new TheoNuclearDensity(*this);
}

/*=============================================
 * Kohm-Sham inverse equations
 *===========================================*/




