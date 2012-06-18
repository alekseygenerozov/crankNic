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
double tidalTorque( 	double r , 		// radial position in disk
											double a, 		// binary separation
											double hh 			// scale height of disk
									){

	if( q == 0.0 ){
		return 0.0;
	}// end simple cas if
	double	tmp1 = f*q*q*M*0.5,
					Dl;

	if( fabs(r-a) > hh )
		Dl = fabs(r-a);
	else
		Dl = hh;

	if( r < a )
		return -tmp1*r*r*r/Dl/Dl/Dl/Dl;
	
	double tmp2 = a/Dl;

	return 	tmp1*tmp2*tmp2*tmp2*tmp2/r;
}// end tidal torque

/*
 *	GAMMA
 *
 *		Dimensionless parameter defined as the ratio of the
 *		torque density's derivative to the torque-density itself,
 *		multiplied by radius.
 */
double gamma(	double r, double a, double hh ){
	if( r < a - hh )
		return 3.0 - 4.0*r/(r-a);
	if( r < a )
		return -1.0;		// !!! ASSUMES h = const * r
	if( r < a + hh )	
		return -5.0;		// !!! ASSUMES h = const * r
	return -1.0 - 4.0*r/(r-a);
}// end gamma

#endif
