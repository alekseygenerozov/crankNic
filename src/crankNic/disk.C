#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "cnSolver.h"
#include "torques.h"
#include "readWrite.h"
#include "initialize.h"
#include "calculateTimeStep.h"
//#include "integrator.h" FIXME

int main(){

	int status = EXIT_SUCCESS,
	    fileCount = 0;

	// read in from parameter file, params.in
	if( EXIT_SUCCESS != (status = readParams()))
		return status;

	// Create arrays for data
	double l[N];			// specific angular momentum
	double Fj[N];			// angular momentum flux

	// Intialize r and sigma
	if(EXIT_SUCCESS != (status = initialize(l,Fj)))
		return status;

	// Initialize timing parameters
	double t         = tStart,
	       dt        = calculateTimeStep(l,Fj,l_a,dl),
	       nextWrite = tStart + tWrite;
	fprintf(stderr,"\t>> tStart = %g\n\t>> initial dt = %g\n\t>> tEnd = %g\n", tStart,dt,tEnd);
	fprintf(stdout,"tStart = %g\ndt = %g\ntEnd = %g\n", tStart,dt,tEnd);

	// print ICs & Parameters we'll use
	if(EXIT_SUCCESS != (status = writeParams()))
		return status;
	if(EXIT_SUCCESS != (status = writeStandard(fileCount++,l,Fj,l_a,t)))
		return status;

	// intialize our Crank-Nicolson solver	
	cnSolver solver;

	while(t<tEnd){		// main loop (in time)

		// take a time step	
		dt = calculateTimeStep(l,Fj,l_a,dl);
		if(EXIT_SUCCESS != (status = solver.step(l,Fj,t,dt,l_a,(t+dt)>=nextWrite))){
			writeStandard(-1,l,Fj,l_a,t);
			return status;
		}
		t += dt;

		// check if we write out
		if( t >= nextWrite ){
			nextWrite += tWrite;
			if(EXIT_SUCCESS != (status = writeStandard(fileCount++,l,Fj,l_a,t)))
				return status;
			fprintf(stderr,"		> dt = %g\n", dt);
		} // end write if
	}// end time-step loop

	// print last file:
	if(EXIT_SUCCESS != (status = writeStandard(fileCount,l,Fj,l_a,t)))
		return status;
	
	return status;
} // end main
