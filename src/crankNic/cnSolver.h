#include <stdio.h>
#include <math.h>

#include "global.h"
#include "nr3.h"
#include "banded.h"
#include "readWrite.h"
#include "torques.h"

#ifndef CN_SOLVER
#define CN_SOLVER

struct cnSolver{
	double 	*d;					// RHS of matrix eq
	double coeffs[9];		// finite difference coefficients
	MatDoub M;					// Matrix of Crank-Nicolson Scheme
	cnSolver();
	int step(double *r,double *sigma,double *sNew,double t,double dt,double &a,double &h);
};

// Constructor
cnSolver::cnSolver() 
	: d(new double[N]),M(N,5)
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
	coeffs[5] = 1.0/tmp1;															// j+1
	coeffs[7] = -l2;																	// j-1
	coeffs[8] = pow(lambda,5)/(lp1*tmp1);							// j-2
	coeffs[6] = -1.0*(coeffs[5]+coeffs[7]+coeffs[8]);	// j
}// end constructor 

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
	static const double L2=0,L1=1,C=2,R1=3,R2=4;

	// Build vectors for matrix solver
	for( int j = 2 ; j < N-3 ; j++ ){
		
		delR = dr/r[j];
		beta = tidal_torque(r[j],a,h)*dt/(omega_k(r[j])*dr2);
		
		tmp0 = pow(lambda,-2.0*j)*alpha;
		tmp1 = pow(lambda,-1.0*j)*delR/4.0*(3.0*alpha-2.0*beta);
		tmp2 = -1.0*beta*delR*delR*(3.0/2.0+gamma(r[j],a,h));

		M[j][L2] = -tmp0*coeffs[4]-tmp1*coeffs[8];						// Second sub-diagonal
		M[j][L1] = -tmp0*coeffs[3]-tmp1*coeffs[7];						// First sub-diagonal
		M[j][C]  = -tmp0*coeffs[2]-tmp1*coeffs[6]-tmp2+1.0;		// central band
		M[j][R1] = -tmp0*coeffs[1]-tmp1*coeffs[5];						// first super-diagonal
		M[j][R2] = -tmp0*coeffs[0];														// second super-diagonal

		// RHS vector
		d[j] = tmp0*coeffs[0]*sigma[j+2] + (tmp0*coeffs[1]+tmp1*coeffs[5])*sigma[j+1]
						+ (tmp0*coeffs[2]+tmp1*coeffs[6]+tmp2+1.0)*sigma[j]
						+ (tmp0*coeffs[3]+tmp1*coeffs[7])*sigma[j-1] + (tmp0*coeffs[4]+tmp1*coeffs[8])*sigma[j-2];
	} // end j for

	// update boundary conditions
	if( ZERO_GRAD == inner_bndry_type ){			// Inner Boundary

		double lp1 = lambda+1.0,l2=lambda*lambda,
			tmp1 = lambda*lambda+lambda+1.0;

		M[0][R1] = lp1;
		M[0][R2] = -1.0/lp1;
		M[0][C]  = -1.0*(M[0][R1]+M[0][R2]);

		M[1][L1] = pow(lambda,3)*(lambda+2.0)/tmp1;
		M[1][R1] = (l2+lambda-1.0)/lambda;
		M[1][R2] = (lambda-1.0)/lambda/tmp1;
		M[1][C]  = -1.0*(M[1][L1]+M[1][R1]+M[1][R2]);

	} else if( DIRICHLET == inner_bndry_type ){
		fprintf(stderr,"ERROR --- Dirichlet Bndry Type not implemented yet.\n");	// FIXME
	} else {
		fprintf(stderr,"ERROR --- Inner Bndry Type Improperly Specified as %d \n",
			inner_bndry_type);
			return EXIT_FAILURE;
	} // end outer BC if/else

	if( ZERO_GRAD == outer_bndry_type ){			// Outer Boundary

		M[N-2][L2] = l2*l2*(lambda-1.0)/tmp1;
		M[N-2][L1] = -lambda*(l2-lambda-1.0);
		M[N-2][R1] = (2.0*lambda+1.0)/tmp1;
		M[N-2][C]  = -1.0*(M[N-2][L2]+M[N-2][L1]+M[N-2][R1]);
 
		M[N-1][L2] = -lambda*lp1;
		M[N-1][L1] = l2*lambda/lp1;
		M[N-1][C]  = -(M[N-1][L2]+M[N-1][L1]);
 
	} else if( DIRICHLET == outer_bndry_type){
		fprintf(stderr,"ERROR --- Dirichlet Bndry Type not implemented yet.\n");	// FIXME
	} else {
		fprintf(stderr,"ERROR --- Outer Bndry Type Improperly Specified as %d \n",
			outer_bndry_type);
		return EXIT_FAILURE;
	} // end outer BC if/else

	solveMatrix(N,L,C,R,d,sNew);
		
	return EXIT_SUCCESS;
} // end solve

#endif
