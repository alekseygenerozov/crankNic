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

	double r[N];
	double sigma[N];
	double sNew[N];		// sigma of t + 1
	
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

	double Lambda[N], Gamma[N];
	for( int i = 0 ; i < N ; i++ ){
		Lambda[i] = lambda( r[i] , 1.0 , .03 );
		Gamma[i]	= gamma(	r[i] , 1.0 , .03 );
	}
	writeOut("TorqueProfile.dat",r,Lambda,N);
	writeOut("GammaProfile.dat", r,Gamma ,N);

	return EXIT_SUCCESS;

	writeOut("T008.dat",r,sigma,N);	// FIXME

 double  tMin = .008/(12.0*nu)*r0*r0,
    tMax = .512/(12.0*nu)*r0*r0,
		t,
		dt= 0.01*dr/nu;
	int Nt = 1 + (int)((tMax-tMin)/dt);
	fprintf(stderr,"\t>> Time Steps: %d\n", Nt);
	bool keepOn = true;

	double alpha = 3.0*nu*dt/(2.0*dr2), delR;

	fprintf(stderr,"alpha = %e \ndelMin = %e\n",alpha,dr/.1);
	
	// Vectors for Crank-Nicolson solver
	double	d[N],	  // RHS of matrix eq
					L[N], 	// left (lower) diagonal of matrix
					R[N],		// right (upper) " "
					C[N];		// central ""

	for( int j = 0 ; j < N ; j++ ){
		delR = dr/r[j];
		L[j] = -alpha*(1.0-0.75*delR);
		C[j] = 1.0 + 2.0*alpha;
		R[j] = -alpha*(1.0+0.75*delR);
	} // end j for

	// boundary conditions (zero-gradient)
	R[0] = 1.0;	C[0] = -1.0;	
	L[N-1] = -1.0; C[N-1] = 1.0;


	for( int i = 0 ; i < Nt && keepOn; i++ ){

		t = i*dt + tMin;

		// ----- SOLVER HERE
		for( int j = 1 ; j < N-1 ; j++ ){
			delR = dr/r[j];
			C[j] = 1.0 + 2.0*alpha;
			d[j] = alpha*(1.0+0.75*delR)*sigma[j+1] + (1.0-2.0*alpha)*sigma[j]
						+ alpha*(1.0-0.75*delR)*sigma[j-1];
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
	      writeOut("ERROR_OUT.dat",r,sNew,N);
	    }// end error if

		// update sigma
	  for( int j = 0 ; j < N ; j++ ){
	    sigma[j] = sNew[j];
	  }


		if( i == (int)(Nt*(32.0-8.0)/512.0)){writeOut("T032.dat",r,sigma,N);}		
		if( i == (int)(Nt*(128.0-8.0)/512.0)){writeOut("T128.dat",r,sigma,N);}

	}// end time-step loop
	
	writeOut("T512.dat",r,sigma,N);

	return status;
}
