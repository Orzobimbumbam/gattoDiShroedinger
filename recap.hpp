#ifndef RECAP_HPP_
#define RECAP_HPP_

class RecapFile
{
public:
	RecapFile(double h, double xin, double xfin);
	void paramRecap();

private:
	RecapFile();
	double m_h, m_xin, m_xfin, m_time;
	int m_loops;

};





#endif /* RECAP_HPP_ */
