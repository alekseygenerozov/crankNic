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
					lh2 = h(l)*h(l),
					tmp2;

	if( l2 < la2 - lh2 ){
		tmp2 = l2/(l2-la2);
		return -tmp1*tmp2*tmp2*tmp2*tmp2/l2;
	}
	if( l2 < la2 ){
		tmp2 = l2/lh2;
		return -tmp1*tmp2*tmp2*tmp2/lh2;
	}
	if( l2 < la2 + lh2 ){
		tmp2 = la2/lh2;
		return tmp1*tmp2*tmp2*tmp2*tmp2/l2;
	}
	tmp2 = la2/(l2-la2);
	return 	tmp1*tmp2*tmp2*tmp2*tmp2/l2;
}// end tidal torque

///*							FIXME
// *	GAMMA
// *
// *		Dimensionless parameter defined as the ratio of the
// *		torque density's derivative to the torque-density itself,
// *		multiplied by radius.
// */
//double gamma(	double r, double a, double hh ){
//	if( r < a - hh )
//		return 3.0 - 4.0*r/(r-a);
//	if( r < a )
//		return -1.0;		// !!! ASSUMES h = const * r
//	if( r < a + hh )	
//		return -5.0;		// !!! ASSUMES h = const * r
//	return -1.0 - 4.0*r/(r-a);
//}// end gamma

#endif
