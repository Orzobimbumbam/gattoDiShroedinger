#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>


class OrbitalsFilling
{
public:
	OrbitalsFilling(int massNumber, const double h);
    void orbFilling(std::vector <std::vector <int> >& orbiMatrix, std::vector<double>& thdensArray, std::vector<int>& quantMax);

private:
    int m_massNumber;
    const double m_h;

};




