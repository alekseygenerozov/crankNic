#include <math.h>

#ifndef INC_GLOBAL
#define INC_GLOBAL

// Problem Types
int problemType = 1;			// delta function, no torque

// NUMERICAL PARAMETERS

const double PI = 3.14159265358979323846; // Pi 
const double twoPI = 2*PI;      // 2*Pi

// TIMING
double tStart = 0.0;
double tEnd = 1.0;
double tWrite = 0.1;

// RESOLUTION PARAMETERS

int N         = 500;      // Size of Simulation
double lambda = 1.0;			// stretch factor of log grid
double rMax   = 2.0;      // outer boundary
double rMin   = 0.1;      // inner boundary
double dr     = (rMax-rMin)/(N-1.0);  // cell size
double dr2    = dr*dr;                // cell size squared
int STENCIL		= 0;				// For Gradient derivative term

// PHYSICAL PARAMETERS
double r0    = 1.0;     					// where delta-fcn starts
double a		 = rMax/2.0;					// Initial position of secondary
double q		 = 0.0;								// binary mass ratio
double M		 = 1.0;								// primary mass
double f		 = .01;								// numerical parameter for torque density
double n_v	 = 0.0;								// viscosity power-law index
double nu0   = -1.0;							// viscosity constant
double dhdr  = 0.1;								// r/h for disk scale height

double max(double a, double b){return (a<b)?a:b;};
double omega_k(double r){ return sqrt(M/(r*r*r));};
double nu(double r){ return (n_v==0?nu0:nu0*pow(r,n_v));};
double h(double r){ return dhdr*r;}

// BOUNDARY CONDITIONS
const int ZERO_GRAD = 0;
const int DIRICHLET = 1;
int outer_bndry_type = 0;		
int inner_bndry_type = 0;
double outer_bndry_value = -1.0;
double inner_bndry_value = -1.0;

// DEBUG PARAMS
int DEBUG_MODE = 0;							// verbose printing
double density_floor = -1.0;		// for when density goes negative

// Problem 3 Params
double p3_A = 0.0;
double p3_B = 1.0;
double p3_C = 0.0;
int	p3_CONST = 0;
double p3_courant = 1.0;

#endif
