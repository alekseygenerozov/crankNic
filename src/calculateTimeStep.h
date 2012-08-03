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
					vDoub_i &l,    // specific angular momentum
					vDoub_i &Fj,   // angular momentum flux
					double l_a,    // binary separation
					double dl      // resolution of smallest grid cell
){

	// manual for problem-type 3
	if( problemType == SQUARE_PULSE )
		return p3_courant*dl;

	double dt=0.0,dtMin=1E5;

	for( int j = 0 ; j < N ; j++ ){

		// viscous diffusion timescale
		dt = .5*dl/Dj(Fj[j],l[j])*pow(lambda,j); //FIXME
		dtMin = min(dtMin,dt);

	} // end j for

	return SAFETY_NUMBER*dtMin;
}// end calculate time step

#endif
