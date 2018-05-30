#pragma once

#include<cmath>
#include <string>


double psi0(unsigned int l);
double psiPrime0(unsigned int l);

namespace Parameters
{

//Global objects - simulation dependent, to be initialised in client code
class ElementConstants
{
public:
    static void initialiseElementConstants(const std::string& fileName, const char delimiter = '_');
    
    static unsigned int NN() {return m_instancePtr -> m_NN;};
    static unsigned int NP() {return m_instancePtr -> m_NP;};
    static unsigned int A() {return m_instancePtr -> m_NP + m_instancePtr -> m_NN;};
    static double Rn() {return m_instancePtr -> m_Rn;};
    static double rms() {return m_instancePtr -> m_rms;};
    static double hBarOmega() {return m_instancePtr -> m_hBarOmega;};
    static std::string elementName() {return m_instancePtr -> m_elementName;};
    
protected: //[Orzobimbumbam] : protected only for testing purposes; fixture classes to be derived to allow pointer reset
    ElementConstants() = default;
    ElementConstants(const std::string& fileName, const char delimiter = '_');
    
    static ElementConstants* m_instancePtr;
    
    unsigned int m_NN, m_NP;
    std::string m_elementName;
    double m_Rn, m_hBarOmega, m_rms;
    
    std::string _extractFileNameKey(const std::string& fileName);
    void _initialiseParameters(const std::string& fileName, const char delimiter = '_');
};

    

class IntegrationParameters
{
public:
    static void initialiseIntegrationParameters(double x0, double x1); //[Orzobimbumbam] : h should be here too ??
    static double x0() {return m_instancePtr -> m_x0;};
    static double x1() {return m_instancePtr -> m_x1;};
    
private:
    IntegrationParameters() = default;
    IntegrationParameters(double x0, double x1);
    
    static IntegrationParameters* m_instancePtr;
    
    double m_x0, m_x1;
};
    
class NucleonType
{
public:
    static void initialiseNucleonType(bool isNeutron);
    static bool isNeutron() {return m_instancePtr -> m_isNeutron;};
    
private:
    NucleonType() = default;
    NucleonType(bool isNeutron);
    
    static NucleonType* m_instancePtr;
    
    bool m_isNeutron; //[Orzobimbumbam] : true for neutron, false for proton simulation type
};

    
//Global constants - immutable, simulation independent

//const int nucleons = 0; 							// set 0 for protons, set 1 for neutrons
//const double mp= 1.6726219e-27;                   // Proton mass [kg]
//const double mn= 1.6749273e-27;					// Neutron mass [kg]
const double mn = 939.565378;                       // Neutron mass in mnc^2 [MeV]
const double mp = 938.28;							// Proton mass in mpc^2 [MeV]
const double R0 = 1.27; 							// [fm]

// Wood-Saxon potential parameters
//const double Rn = R0*pow(A,(1./3.)); 				// Nuclear radius [fm]
const double a0 = 0.67;								// Nuclear surface thickness [fm]


// Runge-Kutta parameters
//const double x_in = 1e-6;							// Initial radius value [fm]
//const double x_fin = 3*Rn;						// Final radius value [fm]


// Other parameters
const double qe = 1.439; //1.6021e-19;				// elementary charge [MeV*fm]
//const double qe = 1.6021e-19;						// elementary charge [C]
//const double qe = 1.2;							// elementary charge [MeV fm^(1/2)]
const double PI = 4*atan(1);
const double hbar = 6.58211928e-16; 				// Reduced constant Planck [eV*s]
const double hbarc = 197.3269788;					// [MeV*fm]
//const double rms = 1.30;							// root-mean-square radius of the charge distribution (<r^2>)^1/2 [fm]
const double rp = 0.8751;							// proton rms charge radius (<rp^2>)^1/2 [fm]
const double pregamma = 0.05;						// gamma prefactor for Kohn-Sham equations

}


