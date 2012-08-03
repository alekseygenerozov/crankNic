#include "global.h"
#include "mr.h"
#include "cubicSpline.h"

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



/*
 *	GAUSS TORQUE INT
 *
 *		Performs Gauss Quadrature of FJ*Lambda/DJ 
 *		in region a < l < b. Needs an interpolation
 *		object for FJ, from which it computes DJ.
 *
 *		GQ coefficients are stored in quadCoeffs.h
 */
double gaussTorqueInt( cubicSpline FJ_cSpline,
                       const double a,
                       const double b )
{
	return 0.0;
}// end gauss Torque Int



/*
 *	MOVE SECONDARY
 *
 *		Performs "back-reaction" of disk-secondary torque
 *		on the secondary blackhole, hardening the binary.
 *
 *		Integrates product of FJ and Torque density
 *		profile, then uses this to find dl_a/dt and updates
 *		l_a for next time step.
 */
void moveSecondary( const double *l,
                    const double *Fj,
										const double dt )
{
	if( secondary == STATIC || q == 0.0 ) return;	// do nothing

	const double LL = 5.0*sqrt(h(l_a));

	// interpolate data
	vDoub ll(N,l), FF(N,Fj);
	cubicSpline FJ_cSpline(ll,FF);

	// perform integrations in four regions, using 
	// 96-pt gaussian quadrature scheme
	double dldt = 0.0;
	dldt += gaussTorqueInt(FJ_cSpline,lMin,l_a-LL); // ( ISCO < l < l_a - LL )
	dldt += gaussTorqueInt(FJ_cSpline,l_a-LL,l_a);  // ( l_a - LL < l < l_a )
	dldt += gaussTorqueInt(FJ_cSpline,l_a,l_a+LL);  // ( l_a < l < l_a + LL )
	dldt += gaussTorqueInt(FJ_cSpline,l_a+LL,lMax); // ( l_a + LL < l < l_out )
	
	dldt *= (-1.0/(M*q));

	// update binary position
	l_a += dldt*dt;

} // end moveSecondary

#endif
