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
		nV += sscanf(line, "q = %lg",&q);
		nV += sscanf(line, "M = %lg",&M);
		nV += sscanf(line, "f = %lg",&f);

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
	fprintf(stderr,"q			= %g\n",q);
	fprintf(stderr,"M			= %g\n",M);
	fprintf(stderr,"f			= %g\n",f);

	return EXIT_SUCCESS;
} // end readParams

int writeOut(	char* fileName,
							int n, 
							double* r, 
							double* f1, 
							double* f2 = NULL	// optional second field
						){
  FILE* fp = fopen(fileName,"w");
  for( int i = 0; i < n ; i++ )
		if(f2)
	    fprintf(fp,"%e\t%e\t%e\n",r[i],f1[i],f2[i]);
		else
	    fprintf(fp,"%e\t%e\n",r[i],f1[i]);
			
  fclose(fp);
  return 0;
}// end writeOut

