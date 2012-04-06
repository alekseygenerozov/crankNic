#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "cnSolver.h"
#include "diagSolvers.h"
#include "torques.h"
#include "readWrite.h"

int main(){
	
	// read in from parameter file, params.in
	int status = readParams();
	if( status != EXIT_SUCCESS ){
		fprintf(stderr,"ERROR READING INPUT FILE\n");
		return EXIT_FAILURE;
	}
  dr   = (rMax-rMin)/(N-1.0);		// update globals 
  dr2  = dr*dr;

	int fileCount = 0;

	// Create arrays for data
	double r[N];				// radial position
	double sigma[N];		// surface density (azimuthally averaged)
	VecDoub sNew(N);		// sigma of current time step

	double a = 1.0;	// binary separation
	double h = .03;	// disk scale height

	// We intialize from file ...
	FILE* fp = fopen("analytic_T0.dat","r");
	if(!fp){
		fprintf(stderr,"ERROR IN DISK.C --- Failed to Open IC file:\n");
		return EXIT_FAILURE;
	} // end error if
	double tmp1,tmp2;
	for( int i = 0 ; i < N ; i++ ){
		fscanf(fp,"%lg",&tmp1);
		fscanf(fp,"%lg",&tmp2);
		r[i] = tmp1;
		sigma[i] = tmp2;
	}// end i for	
	fclose(fp);

	// If Dirichlet BVs and no value set, use IC file:
	if(-1.0==outer_bndry_value)
		outer_bndry_value = sigma[N-1];
	if(-1.0==inner_bndry_value)
		inner_bndry_value = sigma[0];
	
	double tork[N];
	for(int i=0;i<N;i++)
		tork[i] = tidalTorque(r[i],a,h);

	// print ICs & Parameters we'll use
	writeParams();
	writeStandard(fileCount++,N,r,sigma,tork);

	double t,
		dt= 0.01*dr/nu,
		nextWrite = tStart + tWrite;
	int Nt = 1 + (int)((tEnd-tStart)/dt);
	fprintf(stderr,"\t>> Time Steps: %d\n", Nt);
	bool keepOn = true;

	// intialize our Crank-Nicolson solver	
	cnSolver solver;

	// step through time ...
	for( int i = 0 ; i < Nt && keepOn; i++ ){

		t = i*dt + tStart;	
		solver.step(r,sigma,sNew,t,dt,a,h);
				
		// check for negatives
		for( int j = 0 ; j < N ; j++ )
	    if( sNew[j] < 0.0 ){
	      fprintf(stderr,"ERROR: Density negative @ i = %d\n",j);
	      fprintf(stderr,"\t>> t = %g , Nt = %d\n",t,i);
	      keepOn = false;
				for(int j=0;j<N;j++)
					sigma[j]=sNew[j];
	      writeOut("ERROR_OUT.dat",N,r,sigma,tork);
				return EXIT_FAILURE;
	    }// end error if

		// update sigma
	  for( int j = 0 ; j < N ; j++ ){
	    sigma[j] = sNew[j];
	  }

		if( t >= nextWrite ){
			nextWrite += tWrite;
			if(EXIT_SUCCESS != writeStandard(fileCount,N,r,sigma,tork)){
				fprintf(stderr,"ERROR IN SOLVER -- Failed to open output file #%d\n",fileCount);
				return EXIT_FAILURE;
			}// end error if
			fprintf(stderr,"	>> Wrote Output #%d at T = %e\n",fileCount++,t);
		} // end write if

	}// end time-step loop

	// print last file:
	if(EXIT_SUCCESS != writeStandard(fileCount,N,r,sigma,tork)){
		fprintf(stderr,"ERROR IN SOLVER -- Failed to open output file #%d\n",fileCount);
		return EXIT_FAILURE;
	}// end error if
	fprintf(stderr,"  >> Wrote Output #%d at T = %e\n",fileCount++,t);

	
	return status;
} // end main
