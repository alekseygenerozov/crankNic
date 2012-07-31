#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "disk.h"
#include "calculateTimeStep.h"

#ifndef INC_INITIALIZE
#define INC_INITIALIZE 

int initialize( int argc, char **argv, double *l, double *Fj , int &fileCount, double &t ){

	fileCount = 0;
	t = 0.0;

	/*
 	 *	----------- Process command-line arguments
	 */
	const string restart_str  = "-r";
	const string fileInit_str = "-i";
	const string std_err_mssg = "ERROR IN INITIALIZE -- Improper command line inputs";
	const string res_and_init_err = "ERROR IN INITIALIZE -- Cannot specify restart and IC file";

	if( argc > 2 ){

		if( (argc-1) % 2 != 0 ){
			cerr << std_err_mssg << endl;
			return EXIT_FAILURE;
		} // end error if

		for( size_t i = 1 ; i < argc ; i += 2 ){
			if( argv[i][0] == '-' ){
				if( argv[i] == restart_str ){
					if( ! (initial_data_file == "UNSET" )){
						cerr << res_and_init_err << endl;
						return EXIT_FAILURE;
					} // end initial data error
					problemType = RESTART;
					initial_data_file = argv[i+1];
				} else if( argv[i] == fileInit_str ){
					if( problemType == RESTART ){
						cerr << res_and_init_err << endl;
						return EXIT_FAILURE;
					} // end restart error
					initial_data_file = argv[i+1];
					problemType = FROM_FILE;
				} else {
					cerr << std_err_mssg << endl;
					return EXIT_FAILURE;
				}// end case if/else
			} else {
				cerr << std_err_mssg << endl;
				return EXIT_FAILURE;
			} // end argv if/else
		}// end i for

	} else if( argc == 2 ){
		initial_data_file = argv[1];
		problemType = FROM_FILE;
	} // end argc if/else


	// calcualte innermost grid cell size
	if( lambda == 1.0 ){
		dl = (lMax-lMin)/(N-1.0);
	} else {
		dl = (lMax-lMin)*(lambda-1.0)/(pow(lambda,N-1)-1.0);
	}
	dl2  = dl*dl;
	cerr << "dl = " << dl << endl;
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
	 * DELTA FUNCTION PROBLEM
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
	 * RAMPED INITIAL CONDITIONS
	 *
	 *		Steady-state ramp of F_J
	 */
	else if( problemType == RAMPED ){
		double mdot = 3.0*PI/( 1.0 - lMin/lMax );
		cerr << "Mdot = " << mdot << endl;
		for( int j = 0; j < N ; ++j )
			Fj[j] = mdot*( l[j] - lMin );
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
		
		fclose(fp);
	}

	/*
	 *	RESTART
	 *
	 *		Initialize from a previous data dump and continue
	 */
	else if( problemType == RESTART ) {

		int MAX_STRING_LENGTH = 200;
		char line[MAX_STRING_LENGTH];

		FILE *fp = fopen(initial_data_file.c_str(),"r");

		if(!fp){
			cerr << "ERROR IN INITIALIZE.H --- Failed to Open Data File: " << initial_data_file << endl;
			return EXIT_FAILURE;
		}// end error if

		double tmp1, tmp2;
		
		fgets(line,MAX_STRING_LENGTH,fp); // read in fileNum
		sscanf(line,"# N = %d",&fileCount);

		fgets(line,MAX_STRING_LENGTH,fp);	// read in current time
		sscanf(line,"# t = %lg",&t);

		fgets(line,MAX_STRING_LENGTH,fp);	// swallow header line

		for( int j = 0 ; j < N ; ++j ){
			fgets(line,MAX_STRING_LENGTH,fp);
			sscanf(line,"%lg\t%lg",&tmp1,&tmp2);
			Fj[j] = tmp2;
		}// end j for

		fclose(fp);
		cerr << "Restarting at t = " << t << ", fileNum = " << fileCount << endl;
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
			cerr << "WARNING -- Outer bndry value not set\n"
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
			cerr << "WARNING -- Inner bndry value not set.\n" 
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
		cerr << "ERROR IN INITIALIZE: tStart cannot exceed tEnd" << endl;
		return EXIT_FAILURE;
	}// end time error if

	if( problemType != RESTART){	
		t = tStart;
		cerr << "tStart = " << tStart << endl;
	}// end non-restart if 

	cerr << "tEnd = " << tEnd << endl
		<< "tWrite = " << tWrite << endl
		<< "Initial dt = " << calculateTimeStep(l,Fj,l_a,dl) << endl;

	return EXIT_SUCCESS;
} // end initialize

#endif
