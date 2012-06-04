#include <math.h>

#include "global.h"
#include "torques.h" 

#ifndef INC_TIMESTEP
#define INC_TIMESTEP

/*
 *  CALCULATE TIME STEP
 *
 *		Finds optimal time interval to deliver accurate
 *		evolution
 */
double calculateTimeStep( 	
					double *r,				// radial position in disk
					double *sigma,		// surface density profile (currently not needed)
					double a, 				// binary separation
					double dr					// resolution of smallest grid cell
){

	double r2,nu_j,Lambda_j,dt=0.0,omega,dtMin=1E5;

	for( int j = 0 ; j < N ; j++ ){

		r2 = r[j]*r[j];
		nu_j = nu(r[j]);
		Lambda_j = tidalTorque(r[j],a,h(r[j]));
		omega = omega_k(r[j]);

		// viscous diffusion timescale
		dt = .5*dr/nu_j*pow(lambda,2.0*j);
		dtMin = min(dtMin,dt);

		// torque & viscous advection timescale
		dt =  fabs(r[j]*dr/(9.0/4.0*nu_j - Lambda_j/omega))*pow(lambda,1.0*j);
		dtMin = min(dtMin,dt);

		// torque exponential timescale
		dt = fabs(r2*omega/(Lambda_j*(1.5+gamma(r[j],a,h(r[j])))));
		dtMin = min(dtMin,dt);	
	} // end j for

	return SAFETY_NUMBER*dtMin;
}// end calculate time step

#endif
