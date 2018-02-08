#ifndef NLSolverClass_hpp
#define NLSolverClass_hpp

#include <cmath>
#include <iostream>

template <class T, double(T::*evaluate)(double) const, double(T::*fderivative)(double) const = nullptr> //using pointers to class methods allows to parametrise the methods that will be called to evaluate function and first derivative (i.e. can have different names but same signature).
//the last argument is optional, use only for NR algorithm
class NLSolver
{
public:
    NLSolver(double accuracy): m_accuracy(accuracy) {}

    void setAccuracy(double accuracy)
    {
        m_accuracy = accuracy;
    }

    double solveByBisection(const T& f, double targetValue, double a, double b) const //solver by bisection method
    {
        double y1 = (f.*evaluate)(a); //dereferencing a pointer to class method
        double y2 = (f.*evaluate)(b);
        std::cout << a << "\t" << b << std::endl;
        std::cout << y1 << "\t" << y2 << std::endl;
        if (std::abs(y1 - targetValue) <= m_accuracy)
            return y1;
        if (std::abs(y2 - targetValue) <= m_accuracy)
            return y2;
        if ((y1 - targetValue)*(y2 - targetValue) > 0)
            throw "NLSolverClass: solveByBisection : No zeros in this interval.";

        const unsigned long Nmax = 10000;
        unsigned long i = 0;
        double mid = a + 0.5*(b - a);
        double y = (f.*evaluate)(mid);

        while (std::abs(y - targetValue) > m_accuracy)
        {
            if (i == Nmax) //avoid infinite loops
            {
                std::cerr << "NLSolverClass: solveByBisection : No zeros found." << std::endl;
                break;
            }
            if ((y - targetValue)*((f.*evaluate)(a) - targetValue) < 0)
                b = mid;
            else
                a = mid;
            mid = a + 0.5*(b - a);
            y = (f.*evaluate)(mid);
            ++i;
        }
        return mid;
    }

    double solveByNR(const T& f, double x, double targetValue) const //solver by Newton-Raphson method
    {
        double y = (f.*evaluate)(x);

        const unsigned long Nmax = 1000;
        long unsigned i = 0;

        while(std::abs(y - targetValue) > m_accuracy)
        {
            double yprime = (f.*fderivative)(x);
            if (i == Nmax || yprime == 0)
            {
                std::cerr << "NLSolverClass: solveByNR : No zeros found." << std::endl;
                break;
            }
            y = (f.*evaluate)(x);
            x = x - (y - targetValue)/yprime;
            ++i;
        }
        return x;
    }

private:
    const double m_accuracy;

};

#endif /* NLSolverClass_hpp */






