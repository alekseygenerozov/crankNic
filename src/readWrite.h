#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "udSolver.h"

#ifndef INC_READ_PARAMS
#define INC_READ_PARAMS

int readParams( problemDomain &domain, 
                gasDisk &disk, 
                secondaryBH &secondary )
{
	int MAX_STRING_LENGTH = 200;
	char line[MAX_STRING_LENGTH];

	FILE *fp = fopen("params.in","r");
	if( ! fp ){
		cerr << "ERROR IN READ PARAMS -- File params.in cannot be opened" << endl;
		return EXIT_FAILURE;
	}

	int nV = 0, tmp;
	while( fgets(line, MAX_STRING_LENGTH, fp) ){

		nV += sscanf(line, "problemType = %d",&domain.problemType);

		// dimensions
		nV += sscanf(line, "N = %d",&tmp);
			disk.N = tmp;
		nV += sscanf(line, "lambda = %lg",&disk.lambda);			
		nV += sscanf(line, "lMax = %lg",&disk.lMax);			
		nV += sscanf(line, "lMin = %lg",&disk.lMin);
		nV += sscanf(line, "STENCIL = %d",&disk.STENCIL);			

		// physical params
		nV += sscanf(line, "l0 = %lg", &domain.l0);
		nV += sscanf(line, "M = %lg",  &domain.M);

		nV += sscanf(line, "visc_model = %d",  &disk.visc_model);
		nV += sscanf(line, "D0 = %lg",         &disk.D0);
		nV += sscanf(line, "nd = %lg",         &disk.nd);
		nV += sscanf(line, "np = %lg",         &disk.np);
		nV += sscanf(line, "dhdr = %lg",       &disk.dhdr);

		nV += sscanf(line, "l_a = %lg",     &secondary.l_a);
		nV += sscanf(line, "q = %lg",       &secondary.q);
		nV += sscanf(line, "f = %lg",       &secondary.f);
		nV += sscanf(line, "position = %d", &secondary.position);
		nV += sscanf(line, "GW_loss = %d",  &secondary.GW_loss);

		// timing
		nV += sscanf(line, "tStart = %lg",        &domain.tStart);
		nV += sscanf(line, "tEnd = %lg",          &domain.tEnd);
		nV += sscanf(line, "tWrite = %lg",        &domain.tWrite);
		nV += sscanf(line, "SAFETY_NUMBER = %lg", &domain.SAFETY_NUMBER);

		// Boundary Conditions
		nV += sscanf(line, "outer_bndry_type = %d",      &disk.outer_bndry_type);
		nV += sscanf(line, "inner_bndry_type = %d",      &disk.inner_bndry_type);
		nV += sscanf(line, "outer_bndry_laplacian = %d", &disk.outer_bndry_laplacian);
		nV += sscanf(line, "inner_bndry_laplacian = %d", &disk.inner_bndry_laplacian);
		nV += sscanf(line, "outer_bndry_value = %lg",    &disk.outer_bndry_value);
		nV += sscanf(line, "inner_bndry_value = %lg",    &disk.inner_bndry_value);

		// Debug params
		nV += sscanf(line, "DEBUG_MODE = %d", &domain.debug_mode);
		nV += sscanf(line, "density_floor = %lg", &disk.density_floor);

	} // end read while	
	
	fclose(fp);

	return EXIT_SUCCESS;
}// end readParams

/*
 *	WRITE_PARAMS
 *
 *		Writes out parameters we'll use throughout the simulation 
 *		to file params.out
 */
int writeParams( problemDomain &domain, 
                 gasDisk &disk, 
                 secondaryBH &secondary )
{

	FILE *fp;
	if(!(fp=fopen("params.out","w"))){
		cerr << "ERROR IN WRITE PARAMS" << endl 
			<< "	>> Failed to open params.out" << endl;
		return EXIT_FAILURE;
	} // end error if
	
	// Initialization
	if( domain.problemType == FROM_FILE )
		cout << "// Initialized from file " << domain.initial_data_file << endl;

	// Problem Type
	fprintf(fp,"problemType	= %d\n", domain.problemType);
	fprintf(fp,"\n");

	// Dimensions
	fprintf(fp,"N       = %d\n", (int)disk.N);
	fprintf(fp,"lambda  = %g\n", disk.lambda);
	fprintf(fp,"lMax    = %g\n", disk.lMax);
	fprintf(fp,"lMin    = %g\n", disk.lMin);
	fprintf(fp,"STENCIL = %d\n", disk.STENCIL);
	fprintf(fp,"\n");

	// Physical Params
	fprintf(fp,"l0   = %g\n", domain.l0);
	fprintf(fp,"M    = %g\n", domain.M);

	fprintf(fp,"visc_model = %d\n", disk.visc_model);
	fprintf(fp,"D0         = %g\n", disk.D0);
	fprintf(fp,"nd         = %g\n", disk.nd);
	fprintf(fp,"np         = %g\n", disk.np);
	fprintf(fp,"dhdr       = %g\n", disk.dhdr);

	fprintf(fp,"l_a      = %g\n", secondary.l_a);
	fprintf(fp,"q        = %g\n", secondary.q);
	fprintf(fp,"f        = %g\n", secondary.f);
	fprintf(fp,"position = %d\n", secondary.position);
	fprintf(fp,"GW_loss  = %d\n", secondary.GW_loss);
	fprintf(fp,"\n");

	// Timing
	fprintf(fp,"tStart        = %g\n", domain.tStart);
	fprintf(fp,"tEnd          = %g\n", domain.tEnd);
	fprintf(fp,"tWrite        = %g\n", domain.tWrite);
	fprintf(fp,"SAFETY_NUMBER = %g\n", domain.SAFETY_NUMBER);
	fprintf(fp,"\n");

	// Boundary Conditions
	fprintf(fp,"outer_bndry_type      = %d\n", disk.outer_bndry_type);
	fprintf(fp,"inner_bndry_type      = %d\n", disk.inner_bndry_type);
	fprintf(fp,"outer_bndry_laplacian = %d\n", disk.outer_bndry_laplacian);
	fprintf(fp,"inner_bndry_laplacian = %d\n", disk.inner_bndry_laplacian);
	fprintf(fp,"outer_bndry_value     = %g\n", disk.outer_bndry_value);
	fprintf(fp,"inner_bndry_value     = %g\n", disk.inner_bndry_value);
	fprintf(fp,"\n");

	// Debug Params
	fprintf(fp,"DEBUG_MODE    = %d\n", domain.debug_mode);
	fprintf(fp,"density_floor = %g\n", disk.density_floor);
	fprintf(fp,"\n");	

	fclose(fp);

	return EXIT_SUCCESS;
} // end writeParams

// small helper function ...
string intToStr(int i){
	stringstream ss;
	if( i < 10 )
		ss << "00" << i;
	else if( i < 100 )
		ss << "0" << i;
	else 
		ss << i;
	return ss.str();
} // end intToStr

/*
 *	WRITE STANDARD
 *
 *		>> Writes data to file, given # of data dump
 *		   also calculates and prints torque and viscous 
 *		   timescales and mass flow rate
 *
 *		>> Negative datafile # implies an error print
 */
int writeStandard(	problemDomain &domain,
                    gasDisk &disk,
                    secondaryBH &secondary,
										udSolver &solver,
										bool errorPrint = false
){

	int status = EXIT_SUCCESS;
	char str[3];
	string fileName;

	if( errorPrint )
		fileName = "ERROR.dat";
	else
		fileName = "outputFiles/T" + intToStr(domain.fileCount) + ".dat";

  FILE* fp = fopen(fileName.c_str(),"w");
	if(NULL == fp){
		cerr << "ERROR IN WRITE STANDARD" << endl
			<< "	>> File " << fileName << " failed to open." << endl;
		cerr << "\t\t!!! Make sure directory 'outputFiles' exists !!!" << endl; 
		return EXIT_FAILURE;
	}

	// Print current data file #, current time and column headers
	fprintf(fp,"# N = %d\n",domain.fileCount);
	fprintf(fp,"# t = %g\n",domain.t);
	fprintf(fp,"#	l		Fj		Lambda	Mdot\n");

	// print data
	for( size_t j = 0 ; j < disk.N ; j++ ){
		fprintf(fp,"%e\t%e\t%e\t%e\n",disk.l[j],disk.Fj[j],
			secondary.torque(disk,disk.l[j],domain.M),
			solver.Mdot(domain,disk,secondary,j));
	}// end j for

	if(EXIT_SUCCESS == (status = fclose(fp)))
		cerr << "	>> Wrote Output #" << domain.fileCount << " at T = " << domain.t << endl;
	else {
		cerr << "ERROR IN STD WRITE --- file did not close properly" << endl;
		return status;
	}// end error if/else

	domain.writePlus();

	return status;
} // end write standard

#endif
