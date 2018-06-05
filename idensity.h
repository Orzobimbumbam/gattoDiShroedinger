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
    virtual void densityError(); //[Orzobimbumbam] : maybe should return a Density type and remove get?
    virtual Density getDensityError() const;
    virtual bool hasConverged() const;
    
    virtual double distanceToConvergence() const;
    virtual double epsilon() const;
    
    virtual void benchmarkDensity (const std::vector<std::vector<double>>& matrix, double h = 0) = 0;
    virtual Density getBenchmarkDensity() const = 0;
    
    virtual ~NuclearDensity(){};
    
private:
    Density m_thDensity;
    Density m_densError; //[Orzobimbumbam] : is this going to be used elsewhere? If not, perhaps shouldn't be stored
    mutable double m_distanceToConvergenge;
    mutable double m_epsilon;
    //mutable bool m_isFirstLoop;
    
protected:
    Density m_benchmarkDensity;
};



#endif /* IDensity_h */
