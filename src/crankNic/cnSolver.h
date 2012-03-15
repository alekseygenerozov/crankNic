#include <stdio.h>
#include <math.h>

#include "global.h"
#include "diagSolvers.h"
#include "readWrite.h"
#include "torques.h"

#ifndef CN_SOLVER
#define CN_SOLVER

struct cnSolver{
	double 	*d,		// RHS of matrix eq
					*L,		// left (lower) diagonal of matrix
					*R,		// right (upper) ""
					*C;		// central ""
	cnSolver();
	int step(double *r,double *sigma,double *sNew,double t,double dt,double &a,double &h);
};

// Constructor
cnSolver::cnSolver() 
	: d(new double[N]),L(new double[N]),R(new double[N]),C(new double[N]) {} 

/*
 *	STEP
 *
 *		Solves our equations for the next time step,
 *		using a Crank-Nicolson scheme and NR3 matrix solver
 *
 */
int cnSolver::step( 
									double *r, 			// radius
									double *sigma, 	// current surface density
									double *sNew,		// updated surface density
									double t, 			// time
									double dt, 			// width of time step
									double &a,				// binary separation
									double &h				// disk scale height
){

	double delR, beta, alpha = 3.0*nu*dt/(2.0*dr2);

	// Build vectors for matrix solver
	for( int j = 1 ; j < N-1 ; j++ ){
		delR = dr/r[j];
		beta = lambda(r[j],a,h)*dt/(omega_k(r[j])*dr2);

		L[j] = -(alpha*(1.0-0.75*delR)+0.5*beta*delR);
		C[j] = 1.0 + 2.0*alpha + delR*delR*beta*(1.5+gamma(r[j],a,h));
		R[j] = -(alpha*(1.0+0.75*delR)-0.5*beta*delR);

		d[j] = (alpha*(1.0+0.75*delR)-0.5*beta*delR)*sigma[j+1] 
						+ (1.0-2.0*alpha-beta*delR*delR*(1.5+gamma(r[j],a,h)))*sigma[j]
						+ (alpha*(1.0-0.75*delR)+0.5*beta*delR)*sigma[j-1];
	} // end j for

	// update boundary conditions
	if( ZERO_GRAD == inner_bndry_type ){			// Inner Boundary
		R[0] = 1.0;
		C[0] = -1.0;
		d[0] = 0;
	} else if( DIRICHLET == inner_bndry_type ){
		R[0] = 0.0;
		C[0] = 1.0;
		d[0] = inner_bndry_value;
	} else {
		fprintf(stderr,"ERROR --- Inner Bndry Type Improperly Specified as %d \n",
			inner_bndry_type);
			return EXIT_FAILURE;
	} // end outer BC if/else

	if( ZERO_GRAD == outer_bndry_type ){			// Outer Boundary
		L[N-1] = -1.0;
		C[N-1] = 1.0;
		d[N-1] = 0;
	} else if( DIRICHLET == outer_bndry_type){
		L[N-1] = 0.0;
		C[N-1] = 1.0;
		d[N-1] = outer_bndry_value;
	} else {
		fprintf(stderr,"ERROR --- Outer Bndry Type Improperly Specified as %d \n",
			outer_bndry_type);
		return EXIT_FAILURE;
	} // end outer BC if/else

	solveMatrix(N,L,C,R,d,sNew);
		
	return EXIT_SUCCESS;
} // end solve

#endif
