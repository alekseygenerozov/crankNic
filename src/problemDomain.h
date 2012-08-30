#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "mr.h"
#include "unitFactors.h"

using std::string;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

#ifndef INC_DOMAIN
#define INC_DOMAIN

// Problem Types
const int RAMPED = 0;
const int FROM_FILE = 1;
const int RESTART = 2;

const int OFF = 0;
const int ON  = 1;

const double PI = 3.14159265358979323846; // Pi 
const double twoPI = 2*PI;      // 2*Pi

inline double max(const double x, const double y){return (x>y)?x:y;};
inline double min(const double x, const double y){return (x<y)?x:y;};
inline double omega_k(const double l, const double M){ return M*M/(l*l*l);};

class problemDomain {
public:
	problemDomain();
	double update_dt(double new_dt);
	bool isWriteCycle(){ return ( (t + dt) >= nextWrite ); }
	void writePlus(){ fileCount++; nextWrite += tWrite; };
	bool keepOn() const { return t < tEnd; };
	void advance() { t += dt; }

	const unitFactors units;						 // holds units of problem
	string initial_data_file;
	int problemType, debug_mode;         // verbose printing (default off)
	int fileCount;
	double nextWrite;                    // for data IO
	double tStart, tEnd, tWrite, t,dt;   // timing
	double SAFETY_NUMBER;	               // diminishes timestep
	double l0;                           // where delta-fcn starts
	double M;                            // primary mass

}; // end problemDomain

problemDomain::problemDomain()
	: tStart(0.0), tEnd(1.0), tWrite(1.0), SAFETY_NUMBER(1.0), dt(0.0), t(0.0),
	  l0(1.0), M(1.0), debug_mode(OFF), problemType(RAMPED), fileCount(0),
	  initial_data_file("UNSET"), nextWrite(0.1) {}



/*
 *	UPDATE DT
 *		Updates timestep without over-stepping the next write
 */
double problemDomain::update_dt(double new_dt)
{
	dt = new_dt; return dt;	// FIXME
	if( t + new_dt > nextWrite )
		dt = nextWrite - t;
	else
		dt = new_dt;
	return dt;
} // end update dt

#endif
