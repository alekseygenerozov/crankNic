#include <math.h>

#include "global.h"

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
		evenCoeff = (2.0*lambda-1.0)/(6.0*lambda*lambda) + (2.0-lambda)/6.0,
		oddCoeff  = (lambda+1.0)*(lambda+1.0)/(6.0*lambda*lambda);

	result += (2.0-lambda)/6.0*x[0];
	for(int j=1;j<N-1;j++)
			result += ((j%2==0)?evenCoeff:oddCoeff)*x[j]*pow(lambda,j);
	result += (2.0*lambda-1.0)/(6.0*lambda)*pow(lambda,N-4)*x[N-1];

	return result*dr*(lambda+1.0);
}// end simpsonInt 
