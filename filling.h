#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>


class OrbitalsFilling
{
public:
	OrbitalsFilling(int massNumber, const double h);
    void orbFilling(std::vector <std::vector <int> >& orbiMatrix, std::vector<double>& thdensArray);

private:
    int m_massNumber;
    const double m_h;

};




