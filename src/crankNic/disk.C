#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "cnSolver.h"
#include "torques.h"
#include "readWrite.h"
#include "initialize.h"

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
	double mDot[N];			// mass flow (for printing only)
	VecDoub sNew(N);		// sigma of current time step

	// Intialize r and sigma
	status = initialize(r,sigma);
	if( EXIT_SUCCESS != status ) return status;

	// print ICs & Parameters we'll use
	double tork[N];
	for(int j=0;j<N;j++){
		tork[j] = tidalTorque(r[j],a,h(r[j]));
		mDot[j] = 3.0*PI*nu(r[j])*sigma[j];
	}
	writeParams();
	writeStandard(fileCount++,N,r,sigma,tork,mDot);

	// Initialize timing parameters
	double t,
		dt= 0.01*dr/nu(r[0]),
		nextWrite = tStart + tWrite;
	if( problemType == 3 )
		dt = p3_courant*dr;
	int Nt = 1 + (int)((tEnd-tStart)/dt);
	fprintf(stderr,"\t>> dt = %g\n\t>> Time Steps: %d\n", dt,Nt);
	bool keepOn = true;


	// intialize our Crank-Nicolson solver	
	cnSolver solver;

	// step through time ...
	for( int i = 0 ; i < Nt && keepOn; i++ ){

		t = i*dt + tStart;	
		solver.step(r,sigma,sNew,t,dt,a,t>=nextWrite);
				
		// check for negatives
		for( int j = 0 ; j < N ; j++ ){
	    if( sNew[j] < 0.0 ){
				if( density_floor < 0.0 ){
		      fprintf(stderr,"ERROR: Density negative @ i = %d\n",j);
		      fprintf(stderr,"\t>> t = %g , Nt = %d\n",t,i);
		      keepOn = false;
					for(int j=0;j<N;j++){
						sigma[j]=sNew[j];
						mDot[j] = 3.0*PI*nu(r[j])*sigma[j];
					}
		      writeOut("ERROR_OUT.dat",N,r,sigma,tork,mDot);
					return EXIT_FAILURE;
		    } else {
					sNew[j] = density_floor;	// if floor enabled
				} // end floor if/else
			}// end negative density if
		} // end j for

		// update sigma
	  for( int j = 0 ; j < N ; j++ ){
	    sigma[j] = sNew[j];
			mDot[j] = 3.0*PI*nu(r[j])*sigma[j];
	  }

		if( t >= nextWrite ){
			nextWrite += tWrite;
			if(EXIT_SUCCESS != writeStandard(fileCount,N,r,sigma,tork,mDot)){
				fprintf(stderr,"ERROR IN SOLVER -- Failed to open output file #%d\n",fileCount);
				return EXIT_FAILURE;
			}// end error if
			fprintf(stderr,"	>> Wrote Output #%d at T = %e\n",fileCount++,t);
		} // end write if

	}// end time-step loop

	// print last file:
	if(EXIT_SUCCESS != writeStandard(fileCount,N,r,sigma,tork,mDot)){
		fprintf(stderr,"ERROR IN SOLVER -- Failed to open output file #%d\n",fileCount);
		return EXIT_FAILURE;
	}// end error if
	fprintf(stderr,"  >> Wrote Output #%d at T = %e\n",fileCount++,t);

	
	return status;
} // end main
