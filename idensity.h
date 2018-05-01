#ifndef IDensity_h
#define IDensity_h

#include "element.hpp"
#include <map>
#include <fstream>

class Eigenfunction;
class Element;
class NuclearDensity;
typedef std::map<double, double> Density;
typedef std::pair<std::string, NuclearDensity*> NuclearDensityOutputQuery;

class NuclearDensity
{
public:
    NuclearDensity();
    
    virtual void theoreticalDensity(const ElementEigenfunctions& psi, const OrderedLevelDegeneration& degen);
    virtual Density getTheoreticalDensity() const;
    virtual bool hasConverged() const;
    
    virtual double distanceToConvergence() const;
    virtual double epsilon() const;
    
    virtual void benchmarkDensity (const std::vector<std::vector<double>>& matrix, double h = 0) = 0;
    virtual Density getBenchmarkDensity() const = 0;
    //virtual void benchmarkDensity (std::ifstream& file) = 0; //overload for experimental data from file
    
    virtual ~NuclearDensity(){};
    
private:
    Density m_thDensity;
    mutable double m_distanceToConvergenge;
    mutable double m_epsilon;
    mutable bool m_isFirstLoop;
    
protected:
    Density m_benchmarkDensity;
};


#endif /* IDensity_h */
