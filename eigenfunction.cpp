//
//  eigenfunction.cpp
//  Codice
//
//  Created by Alberto Campi on 23/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#include "eigenfunction.hpp"
#include <utility>

double& Eigenfunction::operator()(double key)
{
        
    return m_psi[key];
}

void Eigenfunction::swap(Eigenfunction& rhsEigenfunction)
{
    std::swap(m_psi, rhsEigenfunction.m_psi);
}

Eigenfunction& Eigenfunction::operator=(const Eigenfunction& rhsEigenfunction)
{
    Eigenfunction tempEigenfunction(rhsEigenfunction);
    swap(tempEigenfunction);
    return *this;
}

PsiArrayKVP Eigenfunction::keyValues() const
{
    PsiArrayKVP kvp;
    for (const auto& it : m_psi)
        kvp.push_back(std::make_pair(it.first, it.second));
    
    return kvp;
}
