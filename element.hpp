#ifndef element_hpp
#define element_hpp

#include <vector>
#include "IOUtils.hpp"

class Eigenfunction;
class Schroddy;

typedef std::vector<std::vector<unsigned int>> OrderedOrbitalMatrix;
typedef std::vector<Eigenfunction> ElementEigenfunctions;
typedef std::vector<unsigned int> OrderedLevelDegeneration;
typedef std::vector<std::vector<double>> ElementEigenValues;

class Element
{
public:
    Element(const OrderedOrbitalMatrix& orbitalMatrix);
    
    ElementEigenfunctions orbitalEigenfunction(const Schroddy& sh, const OrderedOrbitalMatrix& orbitalMatrix) const;
    OrderedLevelDegeneration getLevelDegeneration() const;
    ElementEigenValues getLevelEigenvalue() const;
    
private:
    OrderedLevelDegeneration m_levelDegen;
    mutable ElementEigenValues m_eigenValMatrix;
    unsigned int m_orbitalMatrixRows;
};

void writeElementEigenfunctions(const ElementEigenfunctions& elEigf, std::ostream& outStream);

#endif /* element_hpp */
