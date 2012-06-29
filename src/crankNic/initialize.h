#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"

int initialize( double *l, double *Fj ){

	// Check N is odd
//	if(N%2==0){
//		N++;
//		fprintf(stderr,"WARNING in initialize:\n");
//		fprintf(stderr,"\t\t>> Integrator requires odd # of grid cells,");
//		fprintf(stderr,"\t\t   setting N = %d\n",N);
//	}

  // calcualte innermost grid cell size
  if( lambda == 1.0 ){
    dl = (lMax-lMin)/(N-1.0);
  } else {
    dl = (lMax-lMin)*(lambda-1.0)/(pow(lambda,N-1)-1.0);
  }
  dl2  = dl*dl;
  fprintf(stderr,"dl = %g\n",dl);
  fprintf(stdout,"dl = %g\n",dl);

	// setup grid
	for( int j = 0 ; j < N ; j++ ){
      if( lambda == 1.0 ){
        l[j] = lMin + j*dl;
      } else {
        if( j == 0 )
          l[j] = lMin;
        else
          l[j] = lMin + dl*(pow(lambda,j)-1.0)/(lambda-1.0);
      }// end lambda if/else
	}// end j for

	/*
	 * PROBLEM 1
	 *
	 *		Delta function response for S&S Disk,
	 *		constant viscosity, zero torque
	 *
	 *				>> Can be compared to analytic
	 *					 expression
	 */
	if( problemType == DELTA_FCN ){
		// We intialize from file ...
		FILE* fp = fopen("analytic_T0.dat","r");
		if(!fp){
			fprintf(stderr,"ERROR IN INITIALIZE.H --- Failed to Open IC file:\n");
			return EXIT_FAILURE;
		} // end error if
		double tmp1,tmp2;
		for( int i = 0 ; i < N ; i++ ){
			fscanf(fp,"%lg",&tmp1);
			fscanf(fp,"%lg",&tmp2);
			l[i] = tmp1;
			Fj[i] = tmp2;
		}// end i for	
		fclose(fp);
	} // end delta-function initialize

	/*
	 * PROBLEM 2
	 *
	 *		Steady-state ramp of F_J
	 *			F = l			(assumes M-dot infty = 1.0)
	 */
	else if( problemType == RAMPED ){
		for( int j = 0; j < N ; ++j )
			Fj[j] = l[j];
	} // end ramp test problem

	else if( problemType == FROM_FILE ){
		;	//FIXME
	}

	/*
	 *	PROBLEM 4 -- Square Pulse
	 *		To test the equation
	 *			u_t = c/r * u_x	
	 *					(i.e. the advective term)
	 */
//	else if( problemType == SQUARE_PULSE ){
//		for( int i = 0 ; i < N ; i++ ){
//			if( r[i] > (rMax-rMin)*0.7 && r[i] < (rMax-rMin)*(0.8) )
//				sigma[i] = 1.0;
//			else
//				sigma[i] = 0.1;
//		}// end i for
//	}// end square test problem

	// Check outer boundary type/value
	if(-1.0==outer_bndry_value)
		if( outer_bndry_type == DIRICHLET )
			outer_bndry_value = Fj[N-1];
		else if( outer_bndry_type == NEUMANN )
		{
			outer_bndry_value = 0.0;
			fprintf(stderr,"WARNING -- Outer bndry value not set.\n");
			fprintf(stderr,"\t>> Setting value to 0.0 (zero mass flow\n");
		}
	if( DIRICHLET < outer_bndry_type )
	{
		fprintf(stderr,"ERROR IN INITIALIZE.H\n\t>> Outer boundry type improperly set\n");
		return EXIT_FAILURE;
	}
	
	// Check inner boundary type/value
	if(-1.0==inner_bndry_value)
		if( inner_bndry_type == DIRICHLET )
			inner_bndry_value = Fj[0];
		else if( inner_bndry_type == NEUMANN )
		{
			inner_bndry_value = 0.0;
			fprintf(stderr,"WARNING -- Inner bndry value not set.\n");
			fprintf(stderr,"\t>> Setting value to 0.0 (zero mass flow\n");
		}
	if( DIRICHLET < inner_bndry_type )
	{
		fprintf(stderr,"ERROR IN INITIALIZE.H\n\t>> Inner boundry type improperly set\n");
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
} // end initialize
