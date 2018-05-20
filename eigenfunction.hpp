#ifndef eigenfunction_hpp
#define eigenfunction_hpp

#include <map>
#include <vector>

typedef std::map<double,double> Psi;
typedef std::vector<std::pair<double, double>> PsiArrayKVP;

class Eigenfunction
{
    
public:
    double& operator()(double key);
    //const double& operator()(double key);
    Eigenfunction& operator=(const Eigenfunction& rhsEigenfunction);
    PsiArrayKVP keyValues() const;
    Psi get() const;
    
private:
    void _swap(Eigenfunction& rhsEigenfunction);
    Psi m_psi;
};

#endif /* eigenfunction_hpp */
