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
					*L2,		// twice leftward (lower) diagonal of matrix
					*L1,		// once leftward (lower) diagonal of matrix
					*C,			// central ""
					*R1,		// once rightward (upper) ""
					*R2;		// twice rightward (upper) ""
	double coeffs[8];		// finite difference coefficients
	cnSolver();
	int step(double *r,double *sigma,double *sNew,double t,double dt,double &a,double &h);
};

// Constructor
cnSolver::cnSolver() 
	: d(new double[N]),L2(new double[N]),L1(new double[N]),
			C(new double[N]),R1(new double[N]),R2(new double[N])
{
	
	// prep some CN constants based on log-grid stretch factor
	double 	l2 = lambda*lambda,
					lp1 = lambda + 1.0,
					tmp1 = l2+lambda+1.0,
					tmp2 = tmp1*(l2+1.0)*lp1*lp1;

	// for laplacian term ...
	coeffs[0] = -2.0/lambda*(2.0*l2-1.0)/tmp2;									// j+2
	coeffs[1] = 2.0/lambda*lp1*(2.0*lambda-1.0)/tmp1;						// j+1
	coeffs[3] = -2.0*pow(lambda,3)*(lambda-2.0)*lp1/tmp1;				// j-1
	coeffs[4] = 2.0*pow(lambda,7)*(l2-2.0)/tmp2;								// j-2
	coeffs[2] = -1.0*(coeffs[0]+coeffs[1]+coeffs[3]+coeffs[4]);	// j

	// for gradient term ...
	tmp1 = (l2+1.0)/lambda;
	coeffs[5] = 1.0/tmp1;
	coeffs[6] = (l2-1.0)/tmp1;
	coeffs[7] = -l2/tmp1;
} 

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

	double delR, beta, alpha = 3.0*nu*dt/(2.0*dr2),tmp0,tmp1,tmp2;

	// Build vectors for matrix solver
	for( int j = 1 ; j < N-1 ; j++ ){
		
		delR = dr/r[j];
		beta = lambda(r[j],a,h)*dt/(omega_k(r[j])*dr2);
		
		tmp0 = pow(lambda,-2.0*j)*alpha;
		tmp1 = pow(lambda,-1.0*j)*delR/4.0*(3.0*alpha-2.0*beta);
		tmp2 = -1.0*beta*delR*delR*(3.0/2.0+gamma(r[j],a,h));	

		L2[j] = -tmp0*coeffs[4]; 
		L1[j] = -tmp0*coeffs[3]-tmp1*coeffs[7];
		 C[j] = -tmp0*coeffs[2]-tmp1*coeffs[6]-tmp2+1.0;
		R1[j] = -tmp0*coeffs[1]-tmp1*coeffs[5];
		R2[j] = -tmp0*coeffs[0];

		d[j] = tmp0*coeffs[0]*sigma[j+2] + (tmp0*coeffs[1]+tmp1*coeffs[5])*sigma[j+1]
						+ (tmp0*coeffs[2]+tmp1*coeffs[6]+tmp2+1.0)*sigma[j]
						+ (tmp0*coeffs[3]+tmp1*coeffs[7])*sigma[j-1] + tmp0*coeffs[4]*sigma[j-2];

	} // end j for

	// update boundary conditions   FIXME
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
