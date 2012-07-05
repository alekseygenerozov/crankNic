#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "calculateTimeStep.h"

#ifndef INC_INITIALIZE
#define INC_INITIALIZE 

int initialize( double *l, double *Fj , double &t ){

  // calcualte innermost grid cell size
  if( lambda == 1.0 ){
    dl = (lMax-lMin)/(N-1.0);
  } else {
    dl = (lMax-lMin)*(lambda-1.0)/(pow(lambda,N-1)-1.0);
  }
  dl2  = dl*dl;
	cout << "dl = " << dl << endl;
	cerr << "dl = " << dl << endl;

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
			cerr << "ERROR IN INITIALIZE.H --- Failed to Open IC file" << endl;
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

	/*
	 *	FROM_FILE
	 *
	 *		Initializes the surface density from
	 *		a datafile
	 */
	else if( problemType == FROM_FILE ){
		
		FILE *fp = fopen(initial_data_file.c_str(),"r");
		if(!fp){
			cerr << "ERROR IN INITIALIZE.H --- Failed to Open IC file" << initial_data_file << endl;
			return EXIT_FAILURE;
		} // end error if
		
		double tmp1,tmp2;
		for( int j = 0 ; j < N ; ++j ){
			fscanf(fp,"%lg",&tmp1);
			fscanf(fp,"%lg",&tmp2);
			Fj[j] = tmp2;
		}// end i for 
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
			cout << "WARNING -- Outer bndry value not set\n"
				<< "	>> Setting value to 0.0 (zero mass flow" << endl;
		}
	if( DIRICHLET < outer_bndry_type )
	{	
		cerr << "ERROR IN INITIALIZE.H\n"
			<< "	>> Outer boundry type improperly set" << endl;
		return EXIT_FAILURE;
	}
	
	// Check inner boundary type/value
	if(-1.0==inner_bndry_value)
		if( inner_bndry_type == DIRICHLET )
			inner_bndry_value = Fj[0];
		else if( inner_bndry_type == NEUMANN )
		{
			inner_bndry_value = 0.0;
			cout << "WARNING -- Inner bndry value not set.\n" 
				<< "	>> Setting value to 0.0 (zero mass flow" << endl;
		}
	if( DIRICHLET < inner_bndry_type )
	{
		cerr << "ERROR IN INITIALIZE.H\n"
			<< "	>> Inner boundry type improperly set" << endl;
		return EXIT_FAILURE;
	}


	/*
	 *  ------- Initialize Timing
	 */
	if( tStart >= tEnd ){
		cout << "ERROR IN INITIALIZE: tStart cannot exceed tEnd" << endl;
		return EXIT_FAILURE;
	}// end time error if
	
	t = tStart;
	cout << "tStart = " << tStart << endl 
		<< "tEnd = " << tEnd << endl
		<< "tWrite = " << tWrite << endl
		<< "Initial dt = " << calculateTimeStep(l,Fj,l_a,dl) << endl;

	return EXIT_SUCCESS;
} // end initialize

#endif
