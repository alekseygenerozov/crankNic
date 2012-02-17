#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"

int readParams(){

	int MAX_STRING_LENGTH = 200;
	char line[MAX_STRING_LENGTH];

	FILE *fp = fopen("params.dat","r");
	if( ! fp )
		return EXIT_FAILURE;

	int nV = 0;
	while( fgets(line, MAX_STRING_LENGTH, fp) ){

		nV += sscanf(line, "N = %d",&N);			
		nV += sscanf(line, "rMax = %lg",&rMax);			
		nV += sscanf(line, "rMin = %lg",&rMin);
		nV += sscanf(line, "r0 = %lg",&r0);
		nV += sscanf(line, "nu = %lg",&nu);

	} // end read while	
	
	fclose(fp);

	// update global variables dependent on these ...
	dr   = (rMax-rMin)/(N-1.0);
	dr2  = dr*dr;



	fprintf(stderr,"%d variables read from file\n-----------------------\n",nV);
	fprintf(stderr,"N			= %d\n",N);
	fprintf(stderr,"rMax	= %g\n",rMax);
	fprintf(stderr,"rMin	= %g\n",rMin);
	fprintf(stderr,"r0		= %g\n",r0);
	fprintf(stderr,"nu		= %g\n",nu);

	return EXIT_SUCCESS;
} // end readParams

int writeOut(char* fileName, double* r, double* f, int n){
  FILE* fp = fopen(fileName,"w");
  for( int i = 0; i < n ; i++ )
    fprintf(fp,"%g\t%g\n",r[i], f[i]);
  fclose(fp);
  return 0;
}// end writeOut

