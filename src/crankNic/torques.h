#include "global.h"
#include <stdio.h>
#include <math.h>

/*
 *  LAMBDA
 *
 *		Smoothed torque density of secondary on disk, due to
 *		to linblad resonances.
 */
double lambda( 	double r , 		// radial position in disk
								double a, 		// binary separation
								double h 			// scale height of disk
							){
	double	tmp1 = f*q*q*M*0.5,
					Dl;

	if( fabs(r-a) > h )
		Dl = fabs(r-a);
	else
		Dl = h;

	if( r < a )
		return -tmp1*r*r*r/Dl/Dl/Dl/Dl;
	
	double tmp2 = a/Dl;

	return 	tmp1*tmp2*tmp2*tmp2*tmp2/r;
}// end lambda

/*
 *	GAMMA
 *
 *		Dimensionless parameter defined as the ratio of the
 *		torque density's derivative to the torque-density itself,
 *		multiplied by radius.
 */
double gamma(	double r, double a, double h ){
	if( r < a - h )
		return 3.0 - 4.0*r/(r-a);
	if( r < a )
		return 3.0;
	if( r < a + h )
		return -1.0;
	return -1.0 - 4.0*r/(r-a);
}// end gamma
