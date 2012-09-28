#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "calculateTimeStep.h"

#ifndef INC_INITIALIZE
#define INC_INITIALIZE 

int initialize( int argc, 
                char **argv, 
                problemDomain &domain,
                gasDisk &disk,
                secondaryBH &secondary )
{
	
	int status = EXIT_SUCCESS;
	
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
					if( ! ( domain.initial_data_file == "UNSET" )){
						cerr << res_and_init_err << endl;
						return EXIT_FAILURE;
					} // end initial data error
					domain.problemType = RESTART;
					domain.initial_data_file = argv[i+1];
				} else if( argv[i] == fileInit_str ){
					if( domain.problemType == RESTART ){
						cerr << res_and_init_err << endl;
						return EXIT_FAILURE;
					} // end restart error
					domain.initial_data_file = argv[i+1];
					domain.problemType = FROM_FILE;
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
		domain.initial_data_file = argv[1];
		domain.problemType = FROM_FILE;
	} // end argc if/else

	// build grid
	disk.buildGrid();

	/*
	 * RAMPED INITIAL CONDITIONS
	 *
	 *		Steady-state ramp of F_J
	 */
	if( domain.problemType == RAMPED ){
		double mdot = 3.0*PI/( 1.0 - disk.lMin/disk.lMax );
		cerr << "Mdot = " << mdot << endl;
		for( int j = 0; j < disk.N ; ++j )
			disk.Fj[j] = mdot*( disk.l[j] - disk.lMin );
	} // end ramp test problem

	/*
	 *	FROM_FILE
	 *
	 *		Initializes the surface density from
	 *		a datafile
	 */
	else if( domain.problemType == FROM_FILE ){
		
		FILE *fp = fopen(domain.initial_data_file.c_str(),"r");
		if(!fp){
			cerr << "ERROR IN INITIALIZE.H --- Failed to Open IC file" << domain.initial_data_file << endl;
			return EXIT_FAILURE;
		} // end error if
		
		double tmp1,tmp2;
		for( size_t j = 0 ; j < disk.N ; ++j ){
			fscanf(fp,"%lg",&tmp1);
			fscanf(fp,"%lg",&tmp2);
			disk.Fj[j] = tmp2;
		}// end i for 
		
		fclose(fp);
	}

	/*
	 *	RESTART
	 *
	 *		Initialize from a previous data dump and continue
	 */
	else if( domain.problemType == RESTART ) {

		int MAX_STRING_LENGTH = 200;
		char line[MAX_STRING_LENGTH];

		FILE *fp = fopen(domain.initial_data_file.c_str(),"r");

		if(!fp){
			cerr << "ERROR IN INITIALIZE.H --- Failed to Open Data File: " << domain.initial_data_file << endl;
			return EXIT_FAILURE;
		}// end error if

		double tmp1, tmp2;
		
		fgets(line,MAX_STRING_LENGTH,fp); // read in fileNum
		sscanf(line,"# N = %d",&domain.fileCount);

		fgets(line,MAX_STRING_LENGTH,fp);	// read in current time
		sscanf(line,"# t = %lg",&domain.t);

		fgets(line,MAX_STRING_LENGTH,fp);	// swallow header line

		for( int j = 0 ; j < disk.N ; ++j ){
			fgets(line,MAX_STRING_LENGTH,fp);
			sscanf(line,"%lg\t%lg",&tmp1,&tmp2);
			disk.Fj[j] = tmp2;
		}// end j for

		fclose(fp);
		cerr << "Restarting at t = " << domain.t << ", fileNum = " << domain.fileCount << endl;
	}

	// Check outer boundary type/value
	if(UNSET==disk.outer_bndry_value)
		if( disk.outer_bndry_type == DIRICHLET )
			disk.outer_bndry_value = disk.Fj[disk.N-1];
		else if( disk.outer_bndry_type == NEUMANN )
		{
			disk.outer_bndry_value = 0.0;
			cerr << "WARNING -- Outer bndry value not set\n"
				<< "	>> Setting value to 0.0 (zero mass flow" << endl;
		}
	if( DIRICHLET < disk.outer_bndry_type )
	{	
		cerr << "ERROR IN INITIALIZE.H\n"
			<< "	>> Outer boundry type improperly set" << endl;
		return EXIT_FAILURE;
	}
	
	// Check inner boundary type/value
	if(UNSET==disk.inner_bndry_value)
		if( disk.inner_bndry_type == DIRICHLET )
			disk.inner_bndry_value = disk.Fj[0];
		else if( disk.inner_bndry_type == NEUMANN )
		{
			disk.inner_bndry_value = 0.0;
			cerr << "WARNING -- Inner bndry value not set.\n" 
				<< "	>> Setting value to 0.0 (zero mass flow" << endl;
		}
	if( DIRICHLET < disk.inner_bndry_type )
	{
		cerr << "ERROR IN INITIALIZE.H\n"
			<< "	>> Inner boundry type improperly set" << endl;
		return EXIT_FAILURE;
	}

	// Initialize DJ, H and T
	for( int j = 0 ; j < disk.N ; ++j ){
		if( disk.visc_model == PWR_LAW ){
			disk.DJ[j] = disk.D0*pow(disk.l[j],disk.np)*pow(disk.Fj[j],disk.nd);
			disk.H[j]  = disk.dhdr*disk.l[j]*disk.l[j]; // ASSUMES M = 1
		}
		else if( disk.visc_model == BETA_DISK ){  // FIXME !!!
			disk.DJ[j] = 1.0;
			disk.T[j]  = 1.0;
			disk.H[j]  = 1.0;
		} // end visc_model if/else
	}// end j for


	/*
	 *  ------- Initialize Timing
	 */
	if( domain.tStart >= domain.tEnd ){
		cerr << "ERROR IN INITIALIZE: tStart cannot exceed tEnd" << endl;
		return EXIT_FAILURE;
	}// end time error if

	if( domain.problemType != RESTART){	
		domain.t = domain.tStart;
		cerr << "tStart = " << domain.tStart << endl;
	}// end non-restart if 

	domain.nextWrite = domain.tStart;

	cerr << "tEnd = " << domain.tEnd << endl
		<< "tWrite = " << domain.tWrite << endl
		<< "Initial dt = " << calculateTimeStep(domain,disk,secondary) << endl;

	return status;
} // end initialize

#endif
