#include <math.h>

#include "global.h"
//#include "torques.h" FIXME

#ifndef INC_TIMESTEP
#define INC_TIMESTEP

/*
 *  CALCULATE TIME STEP
 *
 *		Finds optimal time interval to deliver accurate
 *		evolution
 */
double calculateTimeStep( 	
					double *l,		// specific angular momentum
					double *Fj,		// angular momentum flux
					double l_a, 		// binary separation
					double dl			// resolution of smallest grid cell
){

	// manual for problem-type 3
	if( problemType == 3 )
		return p3_courant*dl;

	double l2,nu_j,Lambda_j,dt=0.0,omega,dtMin=1E5;

	for( int j = 0 ; j < N ; j++ ){

		l2 = l[j]*l[j];
		nu_j = nu(l[j]);
		//Lambda_j = tidalTorque(r[j],a,h(r[j]));
		//omega = omega_k(r[j]);

		// viscous diffusion timescale
		dt = .5*dl/Dj(l[j])*pow(lambda,j); //FIXME
		dtMin = min(dtMin,dt);

//		// torque & viscous advection timescale
//		dt =  fabs(r[j]*dr/(9.0/4.0*nu_j - Lambda_j/omega))*pow(lambda,1.0*j);
//		dtMin = min(dtMin,dt);
//
//		// torque exponential timescale
//		dt = fabs(r2*omega/(Lambda_j*(1.5+gamma(r[j],a,h(r[j])))))*.001;	// FIXME
//		dtMin = min(dtMin,dt);	
	} // end j for

	return SAFETY_NUMBER*dtMin;
}// end calculate time step

#endif
