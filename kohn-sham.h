#pragma once


class Theoreticaldensity
{
public:
    virtual double theodensity() const = 0;
    virtual ~Theoreticaldensity();

    virtual Theoreticaldensity* clone() const = 0;

private:

};


