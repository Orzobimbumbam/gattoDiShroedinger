#ifndef RECAP_HPP_
#define RECAP_HPP_

class RecapFile
{
public:
	RecapFile(double h, double loops, double time, std::string fSoG, std::string fPot, std::string fDens, int pType, int edType);
	void paramRecap();

private:
	RecapFile();
	const double m_h, m_time;
	const int m_pType, m_edType;
	unsigned long m_loops;

    static std::string m_recapFilePath;
    std::string m_fSoG, m_fPot, m_fDens;
};





#endif /* RECAP_HPP_ */
