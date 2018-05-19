#ifndef element_hpp
#define element_hpp

#include <vector>
#include "IOUtils.hpp"

class Eigenfunction;
class Schroddy;

typedef std::vector<std::vector<double>> OrderedOrbitalMatrix; //[Orzobimbumbam] : should be a double to allow half-integer j-values
typedef std::vector<Eigenfunction> ElementEigenfunctions;
typedef std::vector<unsigned long> OrderedLevelDegeneration;
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
    unsigned long m_orbitalMatrixRows;
    
    unsigned long _getLevelDegeneration(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex, unsigned long quantumAlphaSetSize) const;
    
    unsigned long _getQuantumAlphaSetSize(const OrderedOrbitalMatrix& orbitalMatrix) const;
};

void writeElementEigenfunctions(const ElementEigenfunctions& elEigf, std::ostream& outStream);

#endif /* element_hpp */
