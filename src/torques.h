#include "global.h"

#ifndef INC_TORQUES
#define INC_TORQUES

/*
 *  TIDAL_TORQUE
 *
 *		Smoothed torque density of secondary on disk, due to
 *		to linblad resonances.
 */
double tidalTorque( const double l )
{

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

void moveSecondary( const double *l,
                    const double *Fj )
{
	if( secondary == STATIC || q == 0.0 ) return;	// do nothing

	const double LL = 5.0*sqrt(h(l_a));

	/*
	 *	We split region into four zones:
	 *			1). region to far left of secondary   ( ISCO < l < l_a - LL )
	 *			2). region just to left of secondary  ( l_a - LL < l < l_a )
	 *			3). region just to right of secondary ( l_a < l < l_a + LL )
	 *			4). region to far right of secondary  ( l_a + LL < l < l_out )
	 */

	// split region in 4

	// perform integrations

	// sum them

	// update binary position

} // end moveSecondary

#endif
