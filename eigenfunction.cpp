#include "eigenfunction.hpp"
//#include <utility>

double& Eigenfunction::operator()(double key)
{
        
    return m_psi[key];
}

void Eigenfunction::_swap(Eigenfunction& rhsEigenfunction)
{
    std::swap(m_psi, rhsEigenfunction.m_psi);
}

Eigenfunction& Eigenfunction::operator=(const Eigenfunction& rhsEigenfunction)
{
    Eigenfunction tempEigenfunction(rhsEigenfunction);
    _swap(tempEigenfunction);
    return *this;
}

PsiArrayKVP Eigenfunction::keyValues() const
{
    PsiArrayKVP kvp;
    for (const auto& it : m_psi)
        kvp.push_back(std::make_pair(it.first, it.second));
    
    return kvp;
}

Psi Eigenfunction::get() const
{
    return m_psi;
}


