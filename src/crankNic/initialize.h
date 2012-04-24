#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"

int initialize( double *r, double *sigma ){

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
		for( int i = 0 ; i < N ; i++ ){
			if( lambda == 1.0 ){
				r[i] = rMin + i*dr;
			} else {
		    if( i == 0 )
					r[i] = rMin;
		    else
					r[i] = rMin + dr*( pow(lambda,i) - 1.0 )/(lambda-1.0);
			}// end lambda if/else
			if( r[i] > 0.5 && r[i] < 1.5 )
				sigma[i] = r[i];
			else
				sigma[i] = 0.0;
		}// end i for	
	} // end linTorq test problem

	// If Dirichlet BVs and no value set, use IC file:
	if(-1.0==outer_bndry_value)
		outer_bndry_value = sigma[N-1];
	if(-1.0==inner_bndry_value)
		inner_bndry_value = sigma[0];

	return EXIT_SUCCESS;

} // end main
