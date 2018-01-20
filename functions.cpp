// Class for basic functions
/*
# include <iostream>
# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <cmath>

using namespace std;
 
class functions
{
//==============================================================
// Runge-Kutta method
//==============================================================

void rungekutta (double* (function)(double*, double), double& r, double h, int N, double *u) {
  double *k1 = (*function)(u,r);
  	double* u1=new double[N];
	for(int i=0;i<N;i++){
		u1[i]=u[i]+h*k1[i]/2.;
	}
  double *k2 = (*function)(u1,t+h/2.);
	for(int i=0;i<N;i++){
		u1[i]=u[i]+h*k2[i]/2.;
	}
  double *k3 = (*function)(u1,t+h/2.);
	for(int i=0;i<N;i++){
		u1[i]=u[i]+h*k3[i];
	}
  double *k4 = (*function)(u,r+h);
	for(int i=0;i<N;i++){
		u[i]+=h*(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
	}
  r+=h;
  delete[] x1;
  delete[] k1;
  delete[] k2;
  delete[] k3;
  delete[] k4;
}

//==============================================================
// Trapezes method
//==============================================================

double trapezes(double (*function)(double), double xmin, double xmax, int nstep) {
  double integral=0;
  double h=(xmax-xmin)/nstep;
     integral = h/2*(*function)(xmin)+ h/2*(*function)(xmax);
        for (int i=1; i<nstep; i++) {
            double x = xmin + i*h;
            integral = integral + h*(*function)(x);
        }
  return integral;
}

//===============================================================
// Derivate
//===============================================================

double f(double);
double g(double);
double derivata(double (*function)(double), double, double);

double x=0;
double h=0.001;
/*std::cout << "In quale punto vuoi calcolare la derivata? ";
std::cin >> x;
std::cout << "\n Con quale incremento?";
std::cin >> h;
std::cout << "\n La derivata nel punto scelto e' " << derivata( f, x, h);

char b;
std::cout << "Premere un tasto per terminare.";
std::cin >> b;
return 0;
}

double derivata( double (*fun)(double), double x, double h)
{
return ( fun(x+h) - fun(x) ) / h ;
}

double f(double x)
{
return x*x;
}

}*/
