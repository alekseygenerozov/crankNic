#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "torques.h"
#include "cnSolver.h"

#ifndef INC_READ_PARAMS
#define INC_READ_PARAMS

int readParams(){

	int MAX_STRING_LENGTH = 200;
	char line[MAX_STRING_LENGTH];

	FILE *fp = fopen("params.in","r");
	if( ! fp ){
		cerr << "ERROR IN READ PARAMS -- File params.in cannot be opened" << endl;
		return EXIT_FAILURE;
	}

	int nV = 0;
	while( fgets(line, MAX_STRING_LENGTH, fp) ){

		nV += sscanf(line, "problemType = %d",&problemType);

		// dimensions
		nV += sscanf(line, "N = %d",&N);			
		nV += sscanf(line, "lambda = %lg",&lambda);			
		nV += sscanf(line, "lMax = %lg",&lMax);			
		nV += sscanf(line, "lMin = %lg",&lMin);
		nV += sscanf(line, "STENCIL = %d",&STENCIL);			

		// physical params
		nV += sscanf(line, "l0 = %lg",&l0);
		nV += sscanf(line, "D0 = %lg",&D0);
		nV += sscanf(line, "nd = %lg",&nd);
		nV += sscanf(line, "np = %lg",&np);
		nV += sscanf(line, "dhdr = %lg",&dhdr);
		nV += sscanf(line, "l_a = %lg",&l_a);
		nV += sscanf(line, "q = %lg",&q);
		nV += sscanf(line, "M = %lg",&M);
		nV += sscanf(line, "f = %lg",&f);

		// timing
		nV += sscanf(line, "tStart = %lg",&tStart);
		nV += sscanf(line, "tEnd = %lg",&tEnd);
		nV += sscanf(line, "tWrite = %lg",&tWrite);
		nV += sscanf(line, "SAFETY_NUMBER = %lg",&SAFETY_NUMBER);

		// Boundary Conditions
		nV += sscanf(line, "outer_bndry_type = %d",&outer_bndry_type);
		nV += sscanf(line, "inner_bndry_type = %d",&inner_bndry_type);
		nV += sscanf(line, "outer_bndry_laplacian = %d",&outer_bndry_laplacian);
		nV += sscanf(line, "inner_bndry_laplacian = %d",&inner_bndry_laplacian);
		nV += sscanf(line, "outer_bndry_value = %lg",&outer_bndry_value);
		nV += sscanf(line, "inner_bndry_value = %lg",&inner_bndry_value);

		// Debug params
		nV += sscanf(line, "DEBUG_MODE = %d",&DEBUG_MODE);
		nV += sscanf(line, "density_floor = %lg",&density_floor);

		// Problem 3 Params
		nV += sscanf(line, "p3_A = %lg",&p3_A);
		nV += sscanf(line, "p3_B = %lg",&p3_B);
		nV += sscanf(line, "p3_C = %lg",&p3_C);
		nV += sscanf(line, "p3_courant = %lg",&p3_courant);
		nV += sscanf(line, "p3_CONST = %d",&p3_CONST);

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
int writeParams(){

	FILE *fp;
	if(!(fp=fopen("params.out","w"))){
		cerr << "ERROR IN WRITE PARAMS" << endl 
			<< "	>> Failed to open params.out" << endl;
		return EXIT_FAILURE;
	} // end error if
	
	// Initialization
	if( problemType == FROM_FILE )
		cout << "// Initialized from file " << initial_data_file << endl;

	// Problem Type
	fprintf(fp,"problemType	= %d\n",problemType);
	fprintf(fp,"\n");

	// Dimensions
	fprintf(fp,"N       = %d\n",N);
	fprintf(fp,"lambda  = %g\n",lambda);
	fprintf(fp,"lMax    = %g\n",lMax);
	fprintf(fp,"lMin    = %g\n",lMin);
	fprintf(fp,"STENCIL = %d\n",STENCIL);
	fprintf(fp,"\n");

	// Physical Params
	fprintf(fp,"l0   = %g\n",l0);
	fprintf(fp,"D0   = %g\n",D0);
	fprintf(fp,"nd   = %g\n",nd);
	fprintf(fp,"np   = %g\n",np);
	fprintf(fp,"dhdr = %g\n",dhdr);
	fprintf(fp,"l_a  = %g\n",l_a);
	fprintf(fp,"q    = %g\n",q);
	fprintf(fp,"M    = %g\n",M);
	fprintf(fp,"f    = %g\n",f);
	fprintf(fp,"\n");

	// Timing
	fprintf(fp,"tStart        = %g\n",tStart);
	fprintf(fp,"tEnd          = %g\n",tEnd);
	fprintf(fp,"tWrite        = %g\n",tWrite);
	fprintf(fp,"SAFETY_NUMBER = %g\n",SAFETY_NUMBER);
	fprintf(fp,"\n");

	// Boundary Conditions
	fprintf(fp,"outer_bndry_type      = %d\n",outer_bndry_type);
	fprintf(fp,"inner_bndry_type      = %d\n",inner_bndry_type);
	fprintf(fp,"outer_bndry_laplacian = %d\n",outer_bndry_laplacian);
	fprintf(fp,"inner_bndry_laplacian = %d\n",inner_bndry_laplacian);
	fprintf(fp,"outer_bndry_value     = %g\n",outer_bndry_value);
	fprintf(fp,"inner_bndry_value     = %g\n",inner_bndry_value);
	fprintf(fp,"\n");

	// Debug Params
	fprintf(fp,"DEBUG_MODE    = %d\n",DEBUG_MODE);
	fprintf(fp,"density_floor = %g\n",density_floor);
	fprintf(fp,"\n");	

	// Problem 3 Params	
	fprintf(fp,"p3_A       = %g\n",p3_A);
	fprintf(fp,"p3_B       = %g\n",p3_B);
	fprintf(fp,"p3_C       = %g\n",p3_C);
	fprintf(fp,"p3_courant = %g\n",p3_courant);
	fprintf(fp,"p3_CONST   = %d\n",p3_CONST);	

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
int writeStandard(	int fileNum,     	// datafile #
										double *l,       	// radius
										double *Fj,   		// surface density
										cnSolver &solver,	// matrix solver struct
										double t         	// time
){

	int status = EXIT_SUCCESS;
	char str[3];
	string fileName;

	if( fileNum < 0 )
		fileName = "ERROR.dat";
	else
		fileName = "outputFiles/T" + intToStr(fileNum) + ".dat";

  FILE* fp = fopen(fileName.c_str(),"w");
	if(NULL == fp){
		cerr << "ERROR IN WRITE STANDARD" << endl
			<< "	>> File " << fileName << " failed to open." << endl; 
		return EXIT_FAILURE;
	}

	for( int j = 0 ; j < N ; j++ ){
		fprintf(fp,"%e\t%e\t%e\t%e\n",l[j],Fj[j],tidalTorque(l[j]),solver.Mdot(Fj,l,j));
	}// end j for

	if(EXIT_SUCCESS == (status = fclose(fp)))
		cout << "	>> Wrote Output #" << fileNum << " at T = " << t << endl;
	else {
		cerr << "ERROR IN STD WRITE --- file did not close properly" << endl;
		return status;
	}// end error if/else

	return status;
} // end write standard

#endif