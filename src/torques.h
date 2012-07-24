#include "global.h"

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
					lh2 = h(l)*M,
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

/*
 *  GAMMA
 *
 *    Dimensionless derivative of Torque, defined as:
 *
 *				gamma = Lambda_l * l / Lambda
 */
double gamma( double l ){

	if( q == 0.0 )
		return 0.0;

	double  l2=l*l,la2=l_a*l_a,lh2 = h(l)*M;
	if( l2 < la2 - lh2 )
		return 6.0 - 8.0*l2/(l2-la2);
	if( l2 < la2 )
		return -2.0;
	if( l2 < la2 + lh2 )
		return -10.0;
	return -2.0 - 8.0*l2/(l2-la2);
}// end gamma

#endif
