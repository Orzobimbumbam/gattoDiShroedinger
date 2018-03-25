//
//  eigenfunction.hpp
//  Codice
//
//  Created by Alberto Campi on 23/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

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
    Eigenfunction& operator=(const Eigenfunction& rhsEigenfunction);
    PsiArrayKVP keyValues() const;
    Psi get() const;
    
private:
    void swap(Eigenfunction& rhsEigenfunction);
    Psi m_psi;
};

#endif /* eigenfunction_hpp */
