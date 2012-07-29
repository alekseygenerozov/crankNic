#include <math.h>
#include <string>
#include <sstream>
#include <iostream>

using std::string;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;

#ifndef INC_GLOBAL
#define INC_GLOBAL

// Problem Types
const int DELTA_FCN = 1;
const int RAMPED = 2;
const int FROM_FILE = 3;
const int SQUARE_PULSE = 4;
const int RESTART = 5;
int problemType = DELTA_FCN;

// Some names
string initial_data_file = "UNSET";

// NUMERICAL PARAMETERS

const double PI = 3.14159265358979323846; // Pi 
const double twoPI = 2*PI;      // 2*Pi

// TIMING
double tStart = 0.0;
double tEnd = 1.0;
double tWrite = 0.1;
double SAFETY_NUMBER = 1.0;	// diminishes timestep

// RESOLUTION PARAMETERS

int N         = 500;      // Size of Simulation
double lambda = 1.0;			// stretch factor of log grid
double lMax   = 2.0;      // outer boundary
double lMin   = 0.1;      // inner boundary
double dl     = (lMax-lMin)/(N-1.0);  // cell size
double dl2    = dl*dl;                // cell size squared
int STENCIL		= 0;				// For Gradient derivative term

// PHYSICAL PARAMETERS
double l0    = 1.0;        // where delta-fcn starts
double l_a		 = lMax/2.0; // Initial position of secondary
double q		 = 0.0;        // binary mass ratio
double M		 = 1.0;        // primary mass
double f		 = .01;        // numerical parameter for torque density
double D0    = M*M/16.0;   // Diffsn cnst (dflt: t=tau for const visc case)
double nd    = 0.0;        // diffsn Fj power-law index
double np    = -2.0;       // diffsn l pwr-law indx (dflt = Om_k, const visc)
double dhdr  = 0.1;        // r/h for disk scale height

inline double max(const double x, const double y){return (x>y)?x:y;};
inline double min(const double x, const double y){return (x<y)?x:y;};
inline double omega_k(const double l){ return M*M/(l*l*l);};
inline double h(const double l){ return dhdr*l*l/M;};
inline double Dj(const double Fj, const double l){ 
	return D0*pow(Fj,nd)*pow(l,np);
};
double dmdl(const double Fj, const double l){return Fj/Dj(Fj,l);};

// BOUNDARY CONDITIONS
const int NEUMANN = 0;
const int DIRICHLET = 1;
const int ZERO = 0;
const int SELF_SIM = 1;
const int UNSET = -1;
int outer_bndry_type = NEUMANN;
int inner_bndry_type = NEUMANN;
int outer_bndry_laplacian = ZERO;
int inner_bndry_laplacian = ZERO;
double outer_bndry_value = UNSET;
double inner_bndry_value = UNSET;

// DEBUG PARAMS
int DEBUG_MODE = 0;							// verbose printing
double density_floor = UNSET;		// for when density goes negative

// Problem 3 Params
double p3_A = 0.0;
double p3_B = 1.0;
double p3_C = 0.0;
int	p3_CONST = 0;
double p3_courant = 1.0;

#endif
