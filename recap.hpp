#ifndef RECAP_HPP_
#define RECAP_HPP_

class RecapFile
{
public:
	RecapFile(double h, double loops, double time);
	void paramRecap();

private:
	RecapFile();
	double m_h, m_time;
	int m_loops;

};





#endif /* RECAP_HPP_ */
