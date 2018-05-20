#ifndef element_hpp
#define element_hpp

#include <vector>
#include <tuple>
#include "IOUtils.hpp"

class Eigenfunction;
class Schroddy;

typedef std::vector<std::vector<double>> OrderedOrbitalMatrix; //[Orzobimbumbam] : should be a double to allow half-integer j-values
typedef std::vector<Eigenfunction> ElementEigenfunctions;
typedef std::vector<unsigned long> OrderedLevelDegeneration;
typedef std::vector<std::vector<double>> ElementEigenvalues;
typedef std::tuple<unsigned int, unsigned int, double, double> LevelTuple; //[Orzobimbumbam] : must be (nr, l, j, E) in this given order!

class Element
{
public:
    Element(const OrderedOrbitalMatrix& orbitalMatrix);
    
    ElementEigenfunctions orbitalEigenfunction(const Schroddy& sh, const OrderedOrbitalMatrix& orbitalMatrix) const;
    OrderedLevelDegeneration getLevelDegeneration() const;
    ElementEigenvalues getLevelEigenvalue() const;
    
private:
    const static unsigned long maximumAlphaSetSize;
    const static unsigned long minimumAlphaSetSize;
    
    OrderedLevelDegeneration m_levelDegen;
    mutable ElementEigenvalues m_eigenvalMatrix;
    unsigned long m_orbitalMatrixRows;
    
    unsigned long _getLevelDegeneration(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex) const;
    unsigned long _getQuantumAlphaSetSize(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex) const;
    double _checkAndGetJ(const OrderedOrbitalMatrix& orbitalMatrix, unsigned long levelIndex) const;
    void _addTuple(const LevelTuple& row, unsigned long levelIndex) const;
};

void writeElementEigenfunctions(const ElementEigenfunctions& elEigf, std::ostream& outStream);

#endif /* element_hpp */
