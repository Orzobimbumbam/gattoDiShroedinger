# include "initpot.h"
# include "parameters.h"
# include "schroddy.h"
# include "khon-sham.h"
# include <cmath>


/*============================================
 * Theoretical density
 *==========================================*/
Theoreticaldensity::~Theoreticaldensity() {}

TheoNuclearDensity::TheoNuclearDensity(double normalPsi, double x): m_normalPsi(normalPsi), m_x(x){}

double TheoNuclearDensity::theodensity() const
{
    return (1/(4*Parameters::PI*(m_x*m_x)))*Parameteres::NN*(m_normalPsi*m_normalPsi);// BOH!!!!!!!
}

Theoreticaldensity* TheoNuclearDensity::clone() const
{
    return new TheoNuclearDensity(*this);
}

/*=============================================
 * Kohm-Sham inverse equations
 *===========================================*/




