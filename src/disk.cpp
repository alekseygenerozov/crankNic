#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "cnSolver.h"
#include "torques.h"
#include "readWrite.h"
#include "initialize.h"
#include "calculateTimeStep.h"
#include "mass.h"

//#include "integrator.h" FIXME

int main(int argc, char *argv[]){

	int status = EXIT_SUCCESS,
	    fileCount = 0;

	// read in from parameter file, params.in
	if( EXIT_SUCCESS != (status = readParams()))
		return status;

	// Create arrays for data
	vDoub l(N);                  // specific angular momentum
	vDoub Fj(N);                 // angular momentum flux
	double t=0,dt=0,nextWrite=0;  // timing

	// Intialize r and sigma
	if(EXIT_SUCCESS != (status = initialize(argc,argv,l,Fj,fileCount,t)))
		return status;
	nextWrite = t + tWrite;

	// intialize our Crank-Nicolson solver	
	cnSolver solver;

	// print ICs & Parameters we'll use
	if(EXIT_SUCCESS != (status = writeParams()))
		return status;
	if( problemType != RESTART ){
		if(EXIT_SUCCESS != (status = writeStandard(fileCount++,l,Fj,solver,t)))
			return status;
	} else
		fileCount++;
	
	writeMass(l,Fj,t,solver);

	while(t<tEnd){		// main loop (in time)

		// take a time step	
		dt = calculateTimeStep(l,Fj,l_a,dl);
		if( t + dt >= nextWrite )
			cout << "		>> dt = " << dt << endl;
		if(EXIT_SUCCESS != (status = solver.step(l,Fj,t,dt,l_a,(t+dt)>=nextWrite))){
			writeStandard(-1,l,Fj,solver,t);
			return status;
		}
		t += dt;

		// check if we write out
		if( t >= nextWrite ){
			nextWrite += tWrite;
			if(EXIT_SUCCESS != (status = writeStandard(fileCount++,l,Fj,solver,t)))
				return status;
			writeMass(l,Fj,t,solver);
		} // end write if
	}// end time-step loop

	return status;
} // end main
