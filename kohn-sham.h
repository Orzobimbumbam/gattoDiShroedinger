#pragma once
#include <memory>
#include <string>
#include <vector>

class KohnShamInverse
{
public:
	KohnShamInverse();
    void KSinverse(const std::vector<double>& inTheoDensity, const std::vector<double>& empiDensity,
    	const std::vector<double>& inPot, std::vector<double>& outPot) const;

private:

};



