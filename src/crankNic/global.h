#include <math.h>

#ifndef INC_GLOBAL
#define INC_GLOBAL

// NUMERICAL PARAMETERS

const double PI = 3.14159265358979323846; // Pi 
const double twoPI = 2*PI;      // 2*Pi

// TIMING
double tStart = 0.0;
double tEnd = 1.0;
double tWrite = 0.1;

// RESOLUTION PARAMETERS

int N         = 500;      // Size of Simulation
double rMax   = 2.0;      // outer boundary
double rMin   = 0.1;      // inner boundary
double dr     = (rMax-rMin)/(N-1.0);  // cell size
double dr2    = dr*dr;                // cell size squared

// PHYSICAL PARAMETERS
double r0   = 1.0;      // where delta-fcn starts
double q		= 0.1;			// binary mass ratio
double M		= 1.0;			// primary mass
double f		= .01;			// numerical parameter for torque density
double nu   = 0.1;      // viscosity

double max(double a, double b){return (a<b)?a:b;};
double omega_k(double r){ return sqrt(M/(r*r*r));};

// BOUNDARY CONDITIONS
const int ZERO_GRAD = 0;
const int DIRICHLET = 1;
int outer_bndry_type = 0;		
int inner_bndry_type = 0;
double outer_bndry_value = -1.0;
double inner_bndry_value = -1.0;

#endif
