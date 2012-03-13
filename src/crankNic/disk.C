#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "diagSolvers.h"
#include "readWrite.h"
#include "torques.h"

int main(){

	int status = readParams();
	if( status != EXIT_SUCCESS ){
		fprintf(stderr,"ERROR READING INPUT FILE\n");
		return EXIT_FAILURE;
	}

	int fileCount = 0;

	double r[N];
	double sigma[N];
	double sNew[N];		// sigma of t + 1

	double a = 1.0;	// binary separation FIXME
	double h = .03;	// disk scale height
	
	FILE* fp = fopen("analytic_T0.dat","r");
	double tmp1,tmp2;

	/*
	 *   We intialize from file ...
	 */

	for( int i = 0 ; i < N ; i++ ){
		fscanf(fp,"%lg",&tmp1);
		fscanf(fp,"%lg",&tmp2);
		r[i] = tmp1;
		sigma[i] = tmp2;
	}	
	fclose(fp);

	double tork[N];
	for(int i=0;i<N;i++)
		tork[i] = lambda(r[i],a,h);

	// print ICs
	writeStandard(fileCount++,N,r,sigma,tork);

	double t,
		dt= 0.01*dr/nu,
		nextWrite = tStart + tWrite;
	int Nt = 1 + (int)((tEnd-tStart)/dt);
	fprintf(stderr,"\t>> Time Steps: %d\n", Nt);
	bool keepOn = true;

	double alpha = 3.0*nu*dt/(2.0*dr2), delR, beta;

	fprintf(stderr,"alpha = %e \ndelMin = %e\n",alpha,dr/.1);
	
	// Vectors for Crank-Nicolson solver
	double	d[N],	  // RHS of matrix eq
					L[N], 	// left (lower) diagonal of matrix
					R[N],		// right (upper) " "
					C[N];		// central ""

	for( int j = 0 ; j < N ; j++ ){
		delR = dr/r[j];
		beta = lambda(r[j],a,h)*dt/(omega_k(r[j])*dr2);
		L[j] = -alpha*(1.0-0.75*delR);
		C[j] = 1.0 + 2.0*alpha;
		R[j] = -alpha*(1.0+0.75*delR);
	} // end j for

	// boundary conditions (zero-gradient)
	R[0] = 1.0;	C[0] = -1.0;	
	L[N-1] = -1.0; C[N-1] = 1.0;


	for( int i = 0 ; i < Nt && keepOn; i++ ){

		t = i*dt + tStart;

		// ----- SOLVER HERE
		for( int j = 1 ; j < N-1 ; j++ ){
			delR = dr/r[j];
			beta = lambda(r[j],a,h)*dt/(omega_k(r[j])*dr2);

			L[j] = -(alpha*(1.0-0.75*delR)+0.5*beta*delR);
			C[j] = 1.0 + 2.0*alpha + delR*delR*beta*(1.5+gamma(r[j],a,h));
			R[j] = -(alpha*(1.0+0.75*delR)-0.5*beta*delR);

			d[j] = (alpha*(1.0+0.75*delR)-0.5*beta*delR)*sigma[j+1] 
							+ (1.0-2.0*alpha-beta*delR*delR*(1.5+gamma(r[j],a,h)))*sigma[j]
							+ (alpha*(1.0-0.75*delR)+0.5*beta*delR)*sigma[j-1];
		} // end j for

		// boundary conditions
		C[0] 		= -1.0;
		C[N-1] 	= 1.0;
		d[0] = 0;
		d[N-1] = 0;
	
		solveMatrix(N,L,C,R,d,sNew);
		
		// ----- END SOLVER

		// check for negatives
		for( int j = 0 ; j < N ; j++ )
	    if( sNew[j] < 0.0 ){
	      fprintf(stderr,"ERROR: Density negative @ i = %d\n",j);
	      fprintf(stderr,"\t>> t = %g , Nt = %d\n",t,i);
	      keepOn = false;
	      writeOut("ERROR_OUT.dat",N,r,sNew,tork);
	    }// end error if

		// update sigma
	  for( int j = 0 ; j < N ; j++ ){
	    sigma[j] = sNew[j];
	  }


		if( t >= nextWrite ){
			nextWrite += tWrite;
			if(EXIT_SUCCESS != writeStandard(fileCount++,N,r,sigma,tork)){
				fprintf(stderr,"ERROR IN SOLVER -- Failed to open output file #%d\n",fileCount-1);
				return EXIT_FAILURE;
			}// end error if
		} // end write if

	}// end time-step loop
	
	return status;
}
