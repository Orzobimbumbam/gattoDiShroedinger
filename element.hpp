#ifndef element_hpp
#define element_hpp

#include <vector>
#include "eigenvalues.hpp"

class Eigenfunction;
class Schroddy;

typedef std::vector<std::vector<unsigned int>> OrderedOrbitalMatrix;
typedef std::vector<Eigenfunction> ElementEigenfunctions;
typedef std::vector<unsigned int> OrderedLevelDegeneration;

class Element
{
public:
    Element(const OrderedOrbitalMatrix& orbitalMatrix);
    
    ElementEigenfunctions orbitalEigenfunction(const Schroddy& sh, const OrderedOrbitalMatrix& orbitalMatrix) const;
    OrderedLevelDegeneration getLevelDegeneration() const;
    
private:
    OrderedLevelDegeneration m_levelDegen;
    unsigned int m_orbitalMatrixRows;
};

#endif /* element_hpp */
