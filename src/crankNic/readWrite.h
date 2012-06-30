#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "torques.h"

#ifndef INC_READ_PARAMS
#define INC_READ_PARAMS

int readParams(){

	int MAX_STRING_LENGTH = 200;
	char line[MAX_STRING_LENGTH];

	FILE *fp = fopen("params.in","r");
	if( ! fp )
		return EXIT_FAILURE;

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

	fprintf(stderr,"%d variables read from file\n",nV);

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
		return EXIT_FAILURE;
		fprintf(stderr,"ERROR IN WRITE PARAMS\n	>> Failed to open params.out\n");
	} // end error if
	
	// Initialization
	if( problemType == FROM_FILE )
		fprintf(fp,"// Initialized from file %s\n",initial_data_file.c_str());

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
int intToStr(int i, char *str){
	if( i < 10 )
		sprintf(str,"00%d",i);
	else if( i < 100 )
		sprintf(str,"0%d",i);
	else 
		sprintf(str,"%d",i);
	return EXIT_SUCCESS;
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
int writeStandard(	int fileNum,     // datafile #
										double* l,       // radius
										double* Fj,   // surface density
										double l_a,        // binary separation
										double t         // time
){

	int status = EXIT_SUCCESS;
	char str[3];
	char fileName[100];

	if( fileNum < 0 ){	// error print
		sprintf(fileName,"ERROR.dat");
	} else {	// normal print
		if(EXIT_SUCCESS != intToStr(fileNum,str)){
			return EXIT_FAILURE;
		}// end error if
		sprintf(fileName,"outputFiles/T%s.dat",str);
	}// end i if

  FILE* fp = fopen(fileName,"w");
	if(NULL == fp){
		fprintf(stderr,"ERROR IN WRITE STANDARD\n\t>> File %s failed to open.\n",fileName); 
		return EXIT_FAILURE;
	}

	double tork;
	for( int j = 0 ; j < N ; j++ ){
		tork = tidalTorque(l[j],l_a,h(l[j]));
		fprintf(fp,"%e\t%e\t%e\n",l[j],Fj[j],tork);
	}// end j for

	status = fclose(fp);
	if( EXIT_SUCCESS == status ){
		fprintf(stderr,"	>> Wrote Output #%d at T = %e\n",fileNum,t);
	}// end status print

	return status;
}

#endif
