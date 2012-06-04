#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "cnSolver.h"
#include "torques.h"
#include "readWrite.h"
#include "initialize.h"
#include "calculateTimeStep.h"

int main(){

	// read in from parameter file, params.in
	int status = readParams();
	if( status != EXIT_SUCCESS ){
		fprintf(stderr,"ERROR READING INPUT FILE\n");
		return EXIT_FAILURE;
	}

	// Update global grid params
	if( lambda == 1.0 ){
	  dr = (rMax-rMin)/(N-1.0);
	} else {
		dr = (rMax-rMin)*(lambda-1.0)/(pow(lambda,N-1)-1.0);
	}
	fprintf(stderr,"dr = %g\n",dr);
  dr2  = dr*dr;

	int fileCount = 0;

	// Create arrays for data
	double r[N];				// radial position
	double sigma[N];		// surface density (azimuthally averaged)
	double mDot[N];			// mass flow (for printing only) FIXME
	double tVisc[N];		// timescale of viscosity (printing only)
	double tTork[N];		// timescale of torque (printing only)
	VecDoub sNew(N);		// sigma of current time step

	// Intialize r and sigma
	status = initialize(r,sigma);
	if( EXIT_SUCCESS != status ) return status;

	// print ICs & Parameters we'll use
	double tork[N];
	for(int j=0;j<N;j++){
		tork[j] = tidalTorque(r[j],a,h(r[j]));
		mDot[j] = 3.0*PI*nu(r[j])*sigma[j];
		tVisc[j] = 2.0/3.0*r[j]*r[j]/nu(r[j]);
		tTork[j] = omega_k(r[j])*r[j]*r[j]/fabs(tork[j]);
	}
	writeParams();
	writeStandard(fileCount++,N,r,sigma,tork,mDot,tVisc,tTork);

	// Initialize timing parameters
	double t=tStart,
		dt= calculateTimeStep(r,sigma,a,dr),
		nextWrite = tStart + tWrite;
	if( problemType == 3 )
		dt = p3_courant*dr;
	fprintf(stderr,"\t>> dt = %g\n",dt);
	bool keepOn = true;


	// intialize our Crank-Nicolson solver	
	cnSolver solver;

	// step through time ...
	while( (t<tEnd) && keepOn ){

		dt = calculateTimeStep(r,sigma,a,dr);	
		t += dt;
		solver.step(r,sigma,sNew,t,dt,a,t>=nextWrite);
				
		// check for negatives
		for( int j = 0 ; j < N ; j++ ){
	    if( sNew[j] < 0.0 ){
				if( density_floor < 0.0 ){
		      fprintf(stderr,"ERROR: Density negative @ j = %d\n",j);
		      fprintf(stderr,"\t>> t = %g , tStart = %g, dt= %g \n",t,tStart,dt);
		      keepOn = false;
					for(int j=0;j<N;j++){
						sigma[j]=sNew[j];
					}
		      writeOut("ERROR_OUT.dat",N,r,sigma,tork);
					return EXIT_FAILURE;
		    } else {
					sNew[j] = density_floor;	// if floor enabled
				} // end floor if/else
			}// end negative density if
		} // end j for

		// update sigma
		double r2,nu_j;
	  for( int j = 0 ; j < N ; j++ ){
			r2 = r[j]*r[j];
			nu_j = nu(r[j]);
	    sigma[j] = sNew[j];
			mDot[j] = 3.0*PI*nu_j*sigma[j];
			tVisc[j] = 2.0/3.0*r2/nu_j;
			tTork[j] = omega_k(r[j])*r2/fabs(tork[j]);
	  }

		if( t >= nextWrite ){
			nextWrite += tWrite;
			if(EXIT_SUCCESS != writeStandard(fileCount,N,r,sigma,tork,mDot,tVisc,tTork)){
				fprintf(stderr,"ERROR IN SOLVER -- Failed to open output file #%d\n",fileCount);
				return EXIT_FAILURE;
			}// end error if
			fprintf(stderr,"	>> Wrote Output #%d at T = %e\n",fileCount++,t);
			fprintf(stderr,"		> dt = %g\n", dt);
		} // end write if

	}// end time-step loop

	// print last file:
	if(EXIT_SUCCESS != writeStandard(fileCount,N,r,sigma,tork,mDot,tVisc,tTork)){
		fprintf(stderr,"ERROR IN SOLVER -- Failed to open output file #%d\n",fileCount);
		return EXIT_FAILURE;
	}// end error if
	fprintf(stderr,"	>> Wrote Output #%d at T = %e\n",fileCount++,t);

	
	return status;
} // end main
