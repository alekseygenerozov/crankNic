#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "cnSolver.h"
#include "readWrite.h"
#include "initialize.h"
#include "calculateTimeStep.h"
#include "mass.h"

int main(int argc, char *argv[]){

	int status = EXIT_SUCCESS;

	// intialize domain, secondary & disk
	problemDomain domain;
	gasDisk disk;
	secondaryBH secondary;

	// read in from parameter file, params.in
	if( EXIT_SUCCESS != (status = readParams(domain,disk,secondary)))
		return status;

	// Intialize r and sigma
	if(EXIT_SUCCESS != (status = initialize(argc,argv,domain,disk,secondary)))
		return status;

	// intialize our Crank-Nicolson solver	
	cnSolver solver(disk);

	// print ICs & Parameters we'll use
	if(EXIT_SUCCESS != (status = writeParams(domain,disk,secondary)))
		return status;
	if( domain.problemType != RESTART ){
		if(EXIT_SUCCESS != (status = writeStandard(domain,disk,secondary,solver)))
			return status;
	} else
		domain.writePlus();
	
	writeMass(domain,disk,secondary,solver);

	while(domain.keepOn()){		// main loop (in time)

		// calculate new timestep	
		calculateTimeStep(domain,disk,secondary);
		
		// take a step
		if(EXIT_SUCCESS != (status = solver.step(domain,disk,secondary))){
			writeStandard(domain,disk,secondary,solver,true);
			return status;
		}
		
		// update secondary's position
		secondary.moveSecondary(disk,domain.dt,domain.M);

		domain.advance();

		// check if we write out data
		if( domain.isWriteCycle() ){
			cout << "		>> dt = " << domain.dt << endl;
			if(EXIT_SUCCESS != (status = writeStandard(domain,disk,secondary,solver)))
				return status;
			writeMass(domain,disk,secondary,solver);
		} // end write if

	}// end time-step loop

	return status;
} // end main
