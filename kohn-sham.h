#pragma once
#include <string>
#include <vector>
#include <map>

//typedef std::pair<double, double> KVP;

class KohnShamInverse
{
public:
	KohnShamInverse();
    void KSinverse(const std::vector<double>& inTheoDensity, const std::vector<double>& empiDensity,
    	const std::vector<double>& inPot);

    void getOutPot (std::map<double, double>& outPot) const;

private:
    std::map<double, double> m_outPot;

};



