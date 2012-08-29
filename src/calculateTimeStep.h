#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"

#ifndef INC_TIMESTEP
#define INC_TIMESTEP

/*
 *  CALCULATE TIME STEP
 *
 *		Finds optimal time interval to deliver accurate
 *		evolution
 */
double calculateTimeStep( problemDomain &domain,
                        gasDisk &disk,
                        secondaryBH &secondary )
{
	double dt=0.0, dtMin=1E5;

	for( size_t j = 0 ; j < disk.N ; j++ ){
		dt = .5*disk.dl*disk.dl/disk.Dj(j)*pow(disk.lambda,(int)j); //FIXME
		dtMin = min(dtMin,dt);
	} // end j for

	domain.update_dt( domain.SAFETY_NUMBER*dtMin );

	return domain.dt;
}// end calculate time step

#endif
