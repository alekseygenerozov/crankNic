#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// NUMERICAL PARAMETERS

const double PI = 3.14159265358979323846;	// Pi 
const double twoPI = 2*PI;			// 2*Pi

// RESOLUTION PARAMETERS

int N 		= 500;				// Size of Simulation
double rMax 	= 2.0;		// outer boundary
double rMin	= 0.1;			// inner boundary
double dr		= (rMax-rMin)/(N-1.0);	// cell size
double dr2	= dr*dr; 		// cell size squared

// PHYSICAL PARAMETERS
double r0		= 1.0;			// where delta-fcn starts
double nu		= 0.1;			// viscosity

// FUNCTIONS

void solveMatrix(int,double*,double*,double*,double*,double*);
int writeOut(char*,double*,double*);
double max(double a, double b){return (a<b)?a:b;};
int readParams();

int main(){

	int status = readParams();
	if( status != EXIT_SUCCESS ){
		fprintf(stderr,"ERROR READING INPUT FILE\n");
		return EXIT_FAILURE;
	}

	double r[N];
	double sigma[N];
	double sNew[N];		// sigma of t + 1
	
	FILE* fp = fopen("analytic_T0.dat","r");
	double tmp1,tmp2;

	/*
	 *   We intialize from file ...
	 */

	for( int i = 0 ; i < N ; i++ ){
		fscanf(fp,"%lg",&tmp1);
		fscanf(fp,"%lg",&tmp2);
		r[i] = tmp1;
		sigma[i] = tmp2;
	}	
	fclose(fp);

	writeOut("T008.dat",r,sigma);	// FIXME

 double  tMin = .008/(12.0*nu)*r0*r0,
    tMax = .512/(12.0*nu)*r0*r0,
		t,
		dt= 0.01*dr/nu;
	int Nt = 1 + (int)((tMax-tMin)/dt);
	fprintf(stderr,"\t>> Time Steps: %d\n", Nt);
	bool keepOn = true;

	double alpha = 3.0*nu*dt/(2.0*dr2), delR;

	fprintf(stderr,"alpha = %e \ndelMin = %e\n",alpha,dr/.1);
	
	// Vectors for Crank-Nicolson solver
	double	d[N],	  // RHS of matrix eq
					L[N], 	// left (lower) diagonal of matrix
					R[N],		// right (upper) " "
					C[N];		// central ""

	for( int j = 0 ; j < N ; j++ ){
		delR = dr/r[j];
		L[j] = -alpha*(1.0-0.75*delR);
		C[j] = 1.0 + 2.0*alpha;
		R[j] = -alpha*(1.0+0.75*delR);
	} // end j for

	// boundary conditions (zero-gradient)
	R[0] = 1.0;	C[0] = -1.0;	
	L[N-1] = -1.0; C[N-1] = 1.0;


	for( int i = 0 ; i < Nt && keepOn; i++ ){

		t = i*dt + tMin;

		// ----- SOLVER HERE
		for( int j = 1 ; j < N-1 ; j++ ){
			delR = dr/r[j];
			C[j] = 1.0 + 2.0*alpha;
			d[j] = alpha*(1.0+0.75*delR)*sigma[j+1] + (1.0-2.0*alpha)*sigma[j]
						+ alpha*(1.0-0.75*delR)*sigma[j-1];
		} // end j for

		// boundary conditions
		C[0] 		= -1.0;
		C[N-1] 	= 1.0;
		d[0] = 0;
		d[N-1] = 0;
	
		solveMatrix(N,L,C,R,d,sNew);
		
		// ----- END SOLVER

		// check for negatives
		for( int j = 0 ; j < N ; j++ )
	    if( sNew[j] < 0.0 ){
	      fprintf(stderr,"ERROR: Density negative @ i = %d\n",j);
	      fprintf(stderr,"\t>> t = %g , Nt = %d\n",t,i);
	      keepOn = false;
	      writeOut("ERROR_OUT.dat",r,sNew);
	    }// end error if

		// update sigma
	  for( int j = 0 ; j < N ; j++ ){
	    sigma[j] = sNew[j];
	  }


		if( i == (int)(Nt*(32.0-8.0)/512.0)){writeOut("T032.dat",r,sigma);}		
		if( i == (int)(Nt*(128.0-8.0)/512.0)){writeOut("T128.dat",r,sigma);}
/*		if( i == 2){writeOut("T032.dat",r,sigma);}		
		if( i == 4){writeOut("T064.dat",r,sigma);}		
		if( i == 8){writeOut("T128.dat",r,sigma);}		*/

	}// end time-step loop
	
	writeOut("T512.dat",r,sigma);

	return status;
}

int writeOut(char* fileName, double* r, double* f){
	FILE* fp = fopen(fileName,"w");
	for( int i = 0; i < N ; i++ )
	  fprintf(fp,"%g\t%g\n",r[i], f[i]);
	fclose(fp);
	return 0;
}// end writeOut



void solveMatrix (int n, double *a, double *b, double *c, double *v, double *x)
{

	/*
	 * n - number of equations
	 * a - sub-diagonal -- indexed from 1..n-1
	 * b - the main diagonal
	 * c - sup-diagonal -- indexed from 0..n-2
	 * v - right part
	 * x - the answer
	 */

	for (int i = 1; i < n; i++){
		double m = a[i]/b[i-1];
		b[i] = b[i] - m*c[i-1];
		v[i] = v[i] - m*v[i-1];
	}

	x[n-1] = v[n-1]/b[n-1];

	for (int i = n - 2; i >= 0; i--)
		x[i]=(v[i]-c[i]*x[i+1])/b[i];
	
} // end solve matrix

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
