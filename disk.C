#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// NUMERICAL PARAMETERS

const double PI = 3.14159265358979323846;	// Pi 
const double twoPI = 2*PI;			// 2*Pi

// RESOLUTION PARAMETERS

const int N 		= 500;				// Size of Simulation
const double rMax 	= 2.0;		// outer boundary
const double rMin	= 0.1;			// inner boundary
const double dr		= (rMax-rMin)/(N-1.0);	// cell size
const double dr2	= dr*dr; 		// cell size squared

// PHYSICAL PARAMETERS
const double r0		= 1.0;			// where delta-fcn starts
const double nu		= 0.1;			// viscosity

// FUNCTIONS

void solveMatrix(int,double*,double*,double*,double*,double*);
int writeOut(char*,double*,double*);
double max(double a, double b){return (a<b)?a:b;};


int main(){


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
		sigma[i] = tmp2 + .1;
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
	R[0] = -1.0;	C[0] = 1.0;	
	L[N-1] = -1.0; C[N-1] = 1.0;


	for( int i = 0 ; i < Nt && keepOn; i++ ){

		t = i*dt + tMin;

		// ----- SOLVER HERE
		for( int j = 1 ; j < N-2 ; j++ ){
			delR = dr/r[j];
			d[j] = alpha*(1.0+0.75*delR)*sigma[j+1] + (1.0-2.0*alpha)*sigma[j]
						+ alpha*(1.0-0.75*delR)*sigma[j-1];
		} // end j for

		// boundary conditions
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
		if( i == (int)(Nt*(64.0-8.0)/512.0)){writeOut("T064.dat",r,sigma);}
		if( i == (int)(Nt*(128.0-8.0)/512.0)){writeOut("T128.dat",r,sigma);}
	}// end time-step loop
	
	writeOut("T512.dat",r,sigma);

	return 0;
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

	// make some temp vectors ...
	static double bb[N];
	static double vv[N];
	for(int i=0;i<n;i++){
		bb[i] = b[i];
		vv[i] = v[i];
	}
        /**
         * n - number of equations
         * a - sub-diagonal -- indexed from 1..n-1
         * b - the main diagonal
         * c - sup-diagonal -- indexed from 0..n-2
         * v - right part
         * x - the answer
         */

        for (int i = 1; i < n; i++){
                double m = a[i]/b[i-1];
                bb[i] = b[i] - m*c[i-1];
                vv[i] = v[i] - m*v[i-1];
        }
        x[n-1] = vv[n-1]/bb[n-1];
        for (int i = n - 2; i >= 0; i--)
                x[i]=(vv[i]-c[i]*x[i+1])/bb[i];
} // end solve matrix
