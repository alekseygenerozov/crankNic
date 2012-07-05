#include <stdio.h>
#include <math.h>
#include "global.h"
#include "torques.h"
#include "quadCoeffs.h"

/*
 *	INTEGRATOR
 *
 *		by Munier Salem
 *		July, 2012
 *
 *
 *			Implements Simpson's 1/3rd Rule on a Logarithmic Grid
 *			For Integrating a Known Integrand (must be sufficiently
 *			smooth, and evaluable at both end points -- ie closed --).
 */

double simpsonInt(double* x){

	double result = 0.0,
		evenCoeff = (2.0*lambda-1.0)/(6.0*lambda*lambda*lambda) + (2.0-lambda)/6.0,
		oddCoeff  = (lambda+1.0)*(lambda+1.0)/(6.0*lambda*lambda);

	result += (2.0-lambda)/6.0*x[0];
	for(int j=1;j<N-1;j++)
			result += ((j%2==0)?evenCoeff:oddCoeff)*x[j]*pow(lambda,j);
	result += (2.0*lambda-1.0)/6.0*pow(lambda,N-4)*x[N-1];

	return result*dr*(lambda+1.0);
}// end simpsonInt

/*
 *	BODE INTERP INT WITH TORQUE
 *
 *			by Munier Salem
 *				July, 2012
 *
 *		Integrates data points x multiplied by the tidal torque
 *		density, Lambda. For x, the accuracy is no better than
 *		the trapezoid rule, but for lambda, three extra data points
 *		are added between each cell and the integral is done using
 *		Bode's 7th-order quadrature scheme. This is designed to
 *		better handle the poorly behaved (and poorly resolved) torque
 *		profile and should greatly enhance the accuracy of the integral
 *		provided x (surface density) is smooth.
 */
double bodeInterpIntWithTorque(double *r, double* x){
	
	double result = 0.0;
	double f0,f1,f2,f3,f4,slope,dr_j,rr;

	for(int j=0; j<N-1;j++){
		// find current cell size
		dr_j = dr*pow(lambda,j);	

		// linear interpolate data for 3 pts within cell
		slope = (x[j+1] - x[j])/dr_j;
		rr = r[j  ]; f0 = x[j  ]*tidalTorque(rr,a,h(rr));
		rr = r[j+1]; f4 = x[j+1]*tidalTorque(rr,a,h(rr));
		rr = r[j] + .25*dr_j; f1 = (x[j]+.25*slope)*tidalTorque(rr,a,h(rr));
		rr = r[j] + .50*dr_j; f2 = (x[j]+.50*slope)*tidalTorque(rr,a,h(rr));
		rr = r[j] + .75*dr_j; f3 = (x[j]+.75*slope)*tidalTorque(rr,a,h(rr));

		// Sum integral, using Bode's rule (NR)
		result += (7.0*(f0+f4) + 32.0*(f1+f3) + 12.0*f2)*dr_j;
	}// end j for

	result /= 90.0;

	return result;
}// end bodeInterpIntWithTorque

/*
 *	GAUSS QUADRATURE SCHEME
 *	
 *	Coefficients included in separate file, quadCoeffs.h, and taken
 *	from Abramowitz & Stegun.
 *	
 */

double quadBit(double bit_min, double bit_max,int bit_n){
	double bMao2 = (bit_max - bit_min)/2.0,
		bPao2 = (bit_max + bit_min)/2.0;

	// figure out how many points to use
	// only options are 96, 48 and 16
	double *points;
	double *weights;
	if(bit_n == 96){
		points  = gaussQuad_x_96;
		weights = gaussQuad_w_96;
	} else if(bit_n == 48){
		points  = gaussQuad_x_48;
		weights = gaussQuad_w_48;
	} else {	// fall back on 16
		if(bit_n != 16){
			cout << "WARNING IN integrate.h quadBit()\n"
				<< "	>> bit_n improperly specified, resorting to 16" << endl;
		}// end warning if
		bit_n = 16;
		points  = gaussQuad_x_16;
		weights = gaussQuad_w_16;
	}// end n_bits if/else	

	// integrate using gauss-legendre
	double sum = 0.0,x1,x2;
	int len = bit_n/2;
	for( int i = 0 ; i < len ; i++ ){
		x1 = bPao2 + bMao2*points[i];
		x2 = bPao2 - bMao2*points[i];
		sum += weights[i]*(tidalTorque(x1,a,h(x1))+tidalTorque(x2,a,h(x2)));
	}// end i for
	return sum*bMao2;
}
	
double gaussIntWithTorque(int n_pts){

	// we break our integral into regions ...
	double gauss_scale_mult = 5.0,
		steep_rMin = max(a-gauss_scale_mult*h(a),rMin),
		steep_rMax = min(a+gauss_scale_mult*h(a),rMax),
		sum = 0.0;

	if( rMin < steep_rMin)
		sum += quadBit(rMin,steep_rMin,n_pts);
	if( steep_rMax < rMax )
		sum += quadBit(steep_rMax,rMax,n_pts);
	sum += quadBit(steep_rMin,a,n_pts);
	sum += quadBit(a,steep_rMax,n_pts);

	return sum;
}// end gaussIntWithTorque
