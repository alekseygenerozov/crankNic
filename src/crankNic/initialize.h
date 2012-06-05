#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"

int initialize( double *r, double *sigma ){

  // calcualte innermost grid cell size
  if( lambda == 1.0 ){
    dr = (rMax-rMin)/(N-1.0);
  } else {
    dr = (rMax-rMin)*(lambda-1.0)/(pow(lambda,N-1)-1.0);
  }
  dr2  = dr*dr;
  fprintf(stderr,"dr = %g\n",dr);

	// setup grid
	for( int j = 0 ; j < N ; j++ ){
      if( lambda == 1.0 ){
        r[j] = rMin + j*dr;
      } else {
        if( j == 0 )
          r[j] = rMin;
        else
          r[j] = rMin + dr*(pow(lambda,j)-1.0)/(lambda-1.0);
      }// end lambda if/else
	}// end j for

	// if not explicitly set, normalize viscosity at outer boundary
	if(-1.0==nu0){
		fprintf(stderr,"WARNING --- nu0 not explicitly set, normalizing nu at rMax\n");
		nu0 = (n_v==0?1.0:pow(1.0/r[N-3],n_v));
	}// end nu if

	/*
	 * PROBLEM 1
	 *
	 *		Delta function response for S&S Disk,
	 *		constant viscosity, zero torque
	 *
	 *				>> Can be compared to analytic
	 *					 expression
	 */
	if( problemType == 1 ){
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
			r[i] = tmp1;
			sigma[i] = tmp2;
		}// end i for	
		fclose(fp);
	} // end delta-function initialize

	/*
	 * PROBLEM 2
	 *
	 *		Static Secondary, constant viscosity
	 *			Following Liu & Schapiro
	 */
	else if( problemType == 2 ){
		for( int i = 0 ; i < N ; i++ )
			sigma[i] = 1.0/nu(r[i]);
		sigma[0] = 1E-8;	
	} // end linTorq test problem

	/*
	 *	PROBLEM 3 -- Square Pulse
	 *		To test the equation
	 *			u_t = c/r * u_x	
	 *					(i.e. the advective term)
	 */
	else if( problemType == 3 ){
		for( int i = 0 ; i < N ; i++ ){
			if( r[i] > (rMax-rMin)*0.7 && r[i] < (rMax-rMin)*(0.8) )
				sigma[i] = 1.0;
			else
				sigma[i] = 0.1;
		}// end i for
	}// end square test problem

	// If Dirichlet BVs and no value set, use IC file:
	if(-1.0==outer_bndry_value)
		outer_bndry_value = sigma[N-1];
	if(-1.0==inner_bndry_value)
		inner_bndry_value = sigma[0];

	return EXIT_SUCCESS;

} // end main
