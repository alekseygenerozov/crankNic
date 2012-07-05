#include "global.h"
#include <stdio.h>
#include <math.h>

#ifndef INC_TORQUES
#define INC_TORQUES

/*
 *  TIDAL_TORQUE
 *
 *		Smoothed torque density of secondary on disk, due to
 *		to linblad resonances.
 */
double tidalTorque( double l ){

	if( q == 0.0 ){
		return 0.0;
	}// end simple cas if
	double	tmp1 = f*q*q*M*M*0.5,
					l2 = l*l,
					la2 = l_a*l_a,
					hh = h(l),
					tmp2;

	if( l2 < la2 - hh ){
		tmp2 = l2/(l2-la2);
		return -tmp1*tmp2*tmp2*tmp2*tmp2/l2;
	}
	if( l2 < la2 ){
		tmp2 = l2/hh;
		return -tmp1*tmp2*tmp2*tmp2/hh;
	}
	if( l2 < la2 + hh ){
		tmp2 = la2/hh;
		return tmp1*tmp2*tmp2*tmp2*tmp2/l2;
	}
	tmp2 = la2/(l2-la2);
	return 	tmp1*tmp2*tmp2*tmp2*tmp2/l2;
}// end tidal torque

#endif
