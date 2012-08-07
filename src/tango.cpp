#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "cnSolver.h"
#include "readWrite.h"
#include "initialize.h"
#include "calculateTimeStep.h"
#include "mass.h"

int main(int argc, char *argv[])
{
	int status = EXIT_SUCCESS;

	problemDomain domain;
	gasDisk disk;
	secondaryBH secondary;

	if(EXIT_SUCCESS!=(status=readParams(domain,disk,secondary))) return status;
	if(EXIT_SUCCESS!=(status=initialize(argc,argv,domain,disk,secondary))) return status;

	// intialize our Crank-Nicolson solver	
	cnSolver solver(disk);

	// print ICs & Parameters we'll use
	if(EXIT_SUCCESS != (status = writeParams(domain,disk,secondary))) return status;
	if( domain.problemType != RESTART ){
		if(EXIT_SUCCESS != (status = writeStandard(domain,disk,secondary,solver)))
			return status;
	} else
		domain.writePlus();
	
	writeMass(domain,disk,secondary,solver);

	while(domain.keepOn()){		// main loop (in time)

		calculateTimeStep(domain,disk,secondary);   // calculate new timestep	
		
		if(EXIT_SUCCESS != (status = solver.step(domain,disk,secondary))){	// take a step
			writeStandard(domain,disk,secondary,solver,true);
			return status;
		} // end solver error if
		
		secondary.moveSecondary(disk,domain.dt,domain.M);   // update secondary's position
		domain.advance();                                   // update current time

		if( domain.isWriteCycle() ){		// write out data
			cout << "		>> dt = " << domain.dt << endl;
			if(EXIT_SUCCESS != (status = writeStandard(domain,disk,secondary,solver))) return status;
			writeMass(domain,disk,secondary,solver);
			secondary.writeOut(domain);
		} // end write if

	}// end time-step loop

	return status;
} // end main
