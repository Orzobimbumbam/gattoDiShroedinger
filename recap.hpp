#ifndef RECAP_HPP_
#define RECAP_HPP_

class RecapFile
{
public:
	RecapFile(double h, double loops, double time);
	void paramRecap();

private:
	RecapFile();
	const double m_h, m_time;
	unsigned long m_loops;

    static std::string m_recapFilePath;
};





#endif /* RECAP_HPP_ */
