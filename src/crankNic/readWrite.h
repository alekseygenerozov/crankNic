#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"

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
		nV += sscanf(line, "rMax = %lg",&rMax);			
		nV += sscanf(line, "rMin = %lg",&rMin);

		// physical params
		nV += sscanf(line, "r0 = %lg",&r0);
		nV += sscanf(line, "nu0 = %lg",&nu0);
		nV += sscanf(line, "n_v = %lg",&n_v);
		nV += sscanf(line, "dhdr = %lg",&dhdr);
		nV += sscanf(line, "a = %lg",&a);
		nV += sscanf(line, "q = %lg",&q);
		nV += sscanf(line, "M = %lg",&M);
		nV += sscanf(line, "f = %lg",&f);

		// timing
		nV += sscanf(line, "tStart = %lg",&tStart);
		nV += sscanf(line, "tEnd = %lg",&tEnd);
		nV += sscanf(line, "tWrite = %lg",&tWrite);

		// Boundary Conditions
		nV += sscanf(line, "outer_bndry_type = %d",&outer_bndry_type);
		nV += sscanf(line, "inner_bndry_type = %d",&inner_bndry_type);
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
	
	// Problem Type
	fprintf(fp,"problemType	= %d\n",problemType);
	fprintf(fp,"\n");

	// Dimensions
	fprintf(fp,"N      = %d\n",N);
	fprintf(fp,"lambda = %g\n",lambda);
	fprintf(fp,"rMax   = %g\n",rMax);
	fprintf(fp,"rMin   = %g\n",rMin);
	fprintf(fp,"\n");

	// Physical Params
	fprintf(fp,"r0   = %g\n",r0);
	fprintf(fp,"nu0  = %g\n",nu0);
	fprintf(fp,"n_v  = %g\n",n_v);
	fprintf(fp,"dhdr = %g\n",dhdr);
	fprintf(fp,"a    = %g\n",a);
	fprintf(fp,"q    = %g\n",q);
	fprintf(fp,"M    = %g\n",M);
	fprintf(fp,"f    = %g\n",f);
	fprintf(fp,"\n");

	// Timing
	fprintf(fp,"tStart = %g\n",tStart);
	fprintf(fp,"tEnd   = %g\n",tEnd);
	fprintf(fp,"tWrite = %g\n",tWrite);
	fprintf(fp,"\n");

	// Boundary Conditions
	fprintf(fp,"outer_bndry_type  = %d\n",outer_bndry_type);
	fprintf(fp,"inner_bndry_type  = %d\n",inner_bndry_type);
	fprintf(fp,"outer_bndry_value = %g\n",outer_bndry_value);
	fprintf(fp,"inner_bndry_value = %g\n",inner_bndry_value);
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

int writeOut(	char* fileName,
							int n, 
							double* r, 
							double* f1,
							double* f2 = NULL,	// Optional extra fields 
							double* f3 = NULL, 
							double* f4 = NULL	
						){
  FILE* fp = fopen(fileName,"w");
	if(f4)
		for( int i = 0; i < n ; i++ )
			fprintf(fp,"%e\t%e\t%e\t%e\t%e\n",r[i],f1[i],f2[i],f3[i],f4[i]);
	else if(f3)
		for( int i = 0; i < n ; i++ )
			fprintf(fp,"%e\t%e\t%e\t%e\n",r[i],f1[i],f2[i],f3[i]);
	else if(f2)
		for( int i = 0; i < n ; i++ )
			fprintf(fp,"%e\t%e\t%e\n",r[i],f1[i],f2[i]);
	else
		for( int i = 0; i < n ; i++ )
			fprintf(fp,"%e\t%e\n",r[i],f1[i]);
			
  fclose(fp);
  return EXIT_SUCCESS;
}// end writeOut

int intToStr(int i, char *str){
	if( i < 10 )
		sprintf(str,"00%d",i);
	else if( i < 100 )
		sprintf(str,"0%d",i);
	else 
		sprintf(str,"%d",i);
	return EXIT_SUCCESS;
} // end intToStr

int writeStandard(	int i,
										int n,
										double* r,
										double* f1, 
										double* f2 = NULL,
										double* f3 = NULL,
										double* f4 = NULL ){
	char str[3];
	char fileName[100];
	if(EXIT_SUCCESS != intToStr(i,str)){
		return EXIT_FAILURE;
	}
	sprintf(fileName,"outputFiles/T%s.dat",str);
	return writeOut(fileName,n,r,f1,f2,f3,f4);
}

#endif
