#include <stdio.h>
#include <math.h>
#include <stdio.h>

#include "global.h"
#include "nr3.h"
#include "banded.h"
#include "readWrite.h"
#include "torques.h"

#ifndef CN_SOLVER
#define CN_SOLVER

struct cnSolver{
	VecDoub	d;						// RHS of matrix eq
	VecDoub FjNew;				// angular momentum flux of new timestep
	double coeffs[10];		// finite difference coefficients
	MatDoub M;						// Matrix of Crank-Nicolson Scheme
	cnSolver();
	int step(double *l,double *Fj,double t,double dt,double &l_a,bool dWrite);
};

// Constructor
cnSolver::cnSolver() 
	: d(N),M(N,5),FjNew(N)
{
	
	// prep some CN constants based on log-grid stretch factor
	double 	la2 = lambda*lambda,
					l3 = la2*lambda,l4=l3*lambda,l5=l4*lambda,l6=l3*l3,
					lp1 = lambda+1.0,
					tmp1 = la2+lambda+1.0,
					tmp2 = 1.0+lambda-l3+l5+l6;

	// setup second derivative stencil
	if( 1 == STENCIL ){	// 3rd order, metastable
		double tmp3 = tmp1*(la2+1.0)*lp1*lp1;

		coeffs[0] = -2.0/lambda*(2.0*la2-1.0)/tmp3;                  // j+2
		coeffs[1] = 2.0/lambda*lp1*(2.0*lambda-1.0)/tmp1;           // j+1
		coeffs[3] = -2.0*pow(lambda,3)*(lambda-2.0)*lp1/tmp1;       // j-1
		coeffs[4] = 2.0*pow(lambda,7)*(la2-2.0)/tmp3;                // j-2
		coeffs[2] = -1.0*(coeffs[0]+coeffs[1]+coeffs[3]+coeffs[4]); // j

		if( lambda > 1.3 ){
			cout << "WARNING in cnSolver.h constructor:" << endl
				<< "	>> Stretch factor too large for stability of 3rd order stencil" << endl
				<< "		( try setting lambda < 1.3 or change STENCIL to 0 )" << endl;
		} // end lambda error warning
	} else {	// 2nd order, robustly stable
		coeffs[0] = -2.0*l5*(lambda-1.0)/(lp1*(1.0+l6)*tmp1);									// j+2
		coeffs[1] = 2.0*lambda*(1+2.0*lambda-la2-2.0*l3+2.0*l5+l6)/(lp1*tmp2);	// j+1
		coeffs[3] = 2.0*la2*(1+2.0*lambda-2.0*l3-l4+2.0*l5+l6)/(lp1*tmp2);			// j-2
		coeffs[4] = -coeffs[0];																								// j
		coeffs[2] = -(coeffs[1]+coeffs[3]);																		// j-1
		
		if( 0 != STENCIL ){
			cout << "WARNING in cnSolver.h constructor:" << endl
				<< "	>> Stencil improperly specified, resorting to 2nd order" << endl
				<< "		(i.e. STENCIL = 0)" << endl;
		}// end error if
	} // end STENCIL if for gradient term	

	// for gradient term ...
	coeffs[5] = -1.0/(lambda*lp1*(1.0+la2)*tmp1);								// j+2
	coeffs[6] = lp1/(lambda*tmp1);															// j+1
	coeffs[8] = -pow(lambda,3)*lp1/tmp1;												// j-1
	coeffs[9] = pow(lambda,7)/(lp1*(1.0+la2)*tmp1);							// j-2
	coeffs[7] = -1.0*(coeffs[5]+coeffs[6]+coeffs[8]+coeffs[9]);	// j
	
}// end constructor 

/*
 *	STEP
 *
 *		Solves our equations for the next time step,
 *		using a Crank-Nicolson scheme and NR3 matrix solver
 *
 */
int cnSolver::step( 
									double *l, 			// specific angular momentum
									double *Fj, 		// current a.m. flux
									double t, 			// time
									double dt, 			// width of time step
									double &l_a,			// binary separation
									bool dWrite			// debug write step
){

	int status = EXIT_SUCCESS;
	double tmp0,tmp1,tmp2;
	static const int L2=0,L1=1,C=2,R1=3,R2=4;

	if(DEBUG_MODE && dWrite){
		cout << endl << endl << "# ---------------------------------------------------------" << endl
			<< "#l		D_J		tmp0		tmp1" << endl;
	} // end debug if

	// Build vectors for matrix solver
	for( int j = 2 ; j < N-2 ; j++ ){
	
		tmp0 = pow(lambda,-2.0*j)*.5*dt/dl2*Dj(Fj[j],l[j]);
		tmp1 = 0.0;//pow(lambda,-1.0*j)*delR*(alpha*(2.0*n_v+1.5)-beta); FIXME

		// Coefficiencts manually set for problem-type 3
		if( problemType == SQUARE_PULSE ){
			tmp0 = pow(lambda,-2.0*j)*p3_A;
			if( ! p3_CONST ){
				tmp1 = pow(lambda,-1.0*j)*p3_B/l[j];
			} else {
				tmp1 = pow(lambda,-1.0*j)*p3_B;
			}// end const if
		} // end problem 3 if

		if(DEBUG_MODE && dWrite ){	
			cout << l[j] << "	" << Dj(Fj[j],l[j]) << "	" << tmp0 << "	" << tmp1 << endl;
		}// end debug if


		M[j][L2] = -tmp0*coeffs[4]-tmp1*coeffs[9];     // Second sub-diagonal
		M[j][L1] = -tmp0*coeffs[3]-tmp1*coeffs[8];     // First sub-diagonal
		M[j][C]  = -tmp0*coeffs[2]-tmp1*coeffs[7]+1.0; // central band
		M[j][R1] = -tmp0*coeffs[1]-tmp1*coeffs[6];     // first super-diagonal
		M[j][R2] = -tmp0*coeffs[0]-tmp1*coeffs[5];     // second super-diagonal

		// RHS vector
		d[j] =  (tmp0*coeffs[0]+tmp1*coeffs[5]    )*Fj[j+2] 
		      + (tmp0*coeffs[1]+tmp1*coeffs[6]    )*Fj[j+1]
		      + (tmp0*coeffs[2]+tmp1*coeffs[7]+1.0)*Fj[j  ]
		      + (tmp0*coeffs[3]+tmp1*coeffs[8]    )*Fj[j-1]
		      + (tmp0*coeffs[4]+tmp1*coeffs[9]    )*Fj[j-2];
	} // end j for

	/*
	 * --------- Update boundary conditions ...
	 */
	double lp1 = lambda+1.0,la2=lambda*lambda,
		grad_const = 0.0, laplace_const = 0.0, laplace_val = 0.0;
	tmp1 = lambda*lambda+lambda+1.0, tmp2 = lambda*lambda-lambda-1.0;

	if( NEUMANN == inner_bndry_type ){			// ####### INNER BOUNDS

		// Fixed gradient
		grad_const = 1.0/(dl*lambda);

		M[0][R1] =  grad_const*lp1;
		M[0][R2] = -grad_const/lp1;
		M[0][C]  = -(M[0][R1]+M[0][R2]);
		d[0] = inner_bndry_value;

		// Zero laplacian
		laplace_const = 2.0/(la2*dl2*lp1);
		laplace_val = 0.0;
		if( inner_bndry_laplacian == SELF_SIM && t > 0.0){
			double x = dl/sqrt(4.0*D0*t);
			laplace_val = (1.0-inner_bndry_value)*sqrt(1.0/(PI*D0*t))*exp(-x*x);
		}

		M[1][L1] =  laplace_const*pow(lambda,3)*(lambda+2.0)/tmp1;
		M[1][R1] =  laplace_const*(la2+lambda-1.0)/lambda;
		M[1][R2] = -laplace_const*(lambda-1.0)/lambda/tmp1;
		M[1][C]  = -(M[1][L1]+M[1][R1]+M[1][R2]);
		d[1] = laplace_val;

	} else if( DIRICHLET == inner_bndry_type ){

		// Constant Value
		M[0][R1] = 0.0;
		M[0][R2] = 0.0;
		M[0][C]  = 1.0;
		d[0] = inner_bndry_value;

		// Zero laplacian
		M[1][L1] = pow(lambda,3)*(lambda+2.0)/tmp1;
		M[1][R1] = (la2+lambda-1.0)/lambda;
		M[1][R2] = (lambda-1.0)/lambda/tmp1;
		M[1][C]  = -1.0*(M[1][L1]+M[1][R1]+M[1][R2]);
		d[1] = 0.0;

	} else {
		cerr << "ERROR --- Inner Bndry Type Improperly Specified as " 
			<< inner_bndry_type << endl;
		return EXIT_FAILURE;
	} // end outer BC if/else

	if( NEUMANN == outer_bndry_type ){			// ########### OUTER BOUNDS

		// Laplacian
		laplace_const = 2.0*lambda*pow(lambda,-2.0*(N-2.0))/lp1/dl2;
		laplace_val = 0.0;

		M[N-2][L2] =  laplace_const*la2*la2*(lambda-1.0)/tmp1;
		M[N-2][L1] = -laplace_const*lambda*tmp2;
		M[N-2][R1] =  laplace_const*(2.0*lambda+1.0)/tmp1;
		M[N-2][C]  = -1.0*(M[N-2][L2]+M[N-2][L1]+M[N-2][R1]);
		d[N-2] = laplace_val;
 
		// Fixed gradient
		grad_const = pow(lambda,-(N-1.0))/tmp1/dl;

		M[N-1][L2] = -grad_const*la2*la2/lp1;
		M[N-1][L1] =  grad_const*lp1;
		M[N-1][C]  = -(M[N-1][L2]+M[N-1][L1]);
		d[N-1] = outer_bndry_value;
	
	} else if( DIRICHLET == outer_bndry_type){

		// Zero Laplacian
		M[N-2][L2] = la2*la2*(lambda-1.0)/tmp1;
		M[N-2][L1] = -lambda*tmp2;
		M[N-2][R1] = (2.0*lambda+1.0)/tmp1;
		M[N-2][C]  = -1.0*(M[N-2][L2]+M[N-2][L1]+M[N-2][R1]);
		d[N-2] = 0.0;

		// Constant Value
		M[N-1][L2] = 0.0;
		M[N-1][L1] = 0.0;
		M[N-1][C]  = 1.0;
		d[N-1] = outer_bndry_value;

	} else {
		cerr << "ERROR --- Outer Bndry Type Improperly Specified as " 
			<< outer_bndry_type << endl;
		return EXIT_FAILURE;
	} // end outer BC if/else

	// Solve Matrix
	Bandec banded(M,2,2);
	banded.solve(d,FjNew);

	// Check for negative
	for( int j = 0 ; j < N ; j++ ){
		if( FjNew[j] < 0.0 ){
			if( density_floor < 0.0 ){
				cerr << "ERROR IN CN SOLVER: Density negative @ j = " << j << endl
					<< "	>> t = " << t << ", tStart = " << tStart << ", dt= " << dt << endl;
				status = EXIT_FAILURE;
			} else {
				FjNew[j] = density_floor;  // if floor enabled
				cout << "WARNING IN CN SOLVER: Density negative @ j = " << j << endl 
					<< "	>> t = " << t << ", tStart = " << tStart << ", dt = " << dt << endl
					<< "		Density Floor of " << density_floor << " activated" << endl;
			} // end floor if/else
		}// end negative density if
	} // end j for
	
	// Copy new density into sigma
	for( int j = 0 ; j < N ; j++ )
		Fj[j] = FjNew[j];
	
	return status;
} // end solve

#endif
