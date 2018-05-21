#ifndef energytot_hpp
#define energytot_hpp

#include "element.hpp"

typedef std::map<double, double> Laplacian;
typedef std::vector<Laplacian> Laplacians;

class EnergyTOT
{
public:
    //void energyTot (const Laplacians& LaplacVec, const ElementEigenvalues& elEigV, double h) const;
   // Laplacian get() const;
    double energyTot(const ElementEigenfunctions& elEigf, const ElementEigenvalues& elEigV, double h) const;

private:
    //void laplacian(const ElementEigenfunctions& elEigf, double h)const;
    double _laplacian(std::map<double, double>::const_iterator itt) const;
    //mutable Laplacian m_LaplacMap;
    //mutable Laplacians m_LaplacVec;

};

#endif /* energytot_hpp */


