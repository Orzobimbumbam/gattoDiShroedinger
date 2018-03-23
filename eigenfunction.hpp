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

typedef std::map<double,double> Psi;

class Eigenfunction
{
    
public:
    double& operator()(double key);
    Eigenfunction& operator=(const Eigenfunction& rhsEigenfunction);
    
private:
    void swap(Eigenfunction& rhsEigenfunction);
    Psi m_psi;
};

#endif /* eigenfunction_hpp */
