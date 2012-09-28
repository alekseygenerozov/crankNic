#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "mr.h"
#include "nr3.h"
#include "banded.h"

#ifndef CN_SOLVER
#define CN_SOLVER

class cnSolver{
public:
	cnSolver(const gasDisk& disk );
	int step(problemDomain &domain, gasDisk &disk, secondaryBH &secondary);
	double Mdot(const problemDomain &domain, const gasDisk &disk, 
		const secondaryBH &secondary, const int j ) const;
private:
	static const unsigned int STENCIL_SIZE = 5;              // # of cells in stencil (per time step)
	static const unsigned int CNTR = 2;                      // center of stencil (jth grid cell)
	static const unsigned int JP2=4,JP1=3,J=2,JM1=1,JM2=0,   // indexing for stencil
    R2=JP2,R1=JP1,C=J,L1=JM1,L2=JM2;                // indexing for Crank Nicolson matrix

	VecDoub	d;                           // RHS of matrix eq
	VecDoub FjNew;                       // angular momentum flux of new timestep
	double grad_coeffs[STENCIL_SIZE];    // finite diff coefficients for 1st deriv
	double laplace_coeffs[STENCIL_SIZE]; // ""                           2nd deriv
	double const_coeffs[STENCIL_SIZE];   // ""                           0th deriv
	MatDoub M;                           // Matrix of Crank-Nicolson Scheme
};

// Default Constructor
cnSolver::cnSolver(const gasDisk& disk) 
	: d(disk.N),M(disk.N,STENCIL_SIZE),FjNew(disk.N)
{
	const int N = disk.N;
	const double lambda = disk.lambda;
	const int STENCIL = disk.STENCIL;
	
	// prep some CN constants based on log-grid stretch factor
	double 	la2 = lambda*lambda,
					l3 = la2*lambda,l4=l3*lambda,l5=l4*lambda,l6=l3*l3,
					lp1 = lambda+1.0,
					tmp1 = la2+lambda+1.0,
					tmp2 = 1.0+lambda-l3+l5+l6;

	// setup second derivative stencil
	if( 1 == STENCIL ){	// 3rd order, metastable
		double tmp3 = tmp1*(la2+1.0)*lp1*lp1;

		laplace_coeffs[JP2] = -2.0/lambda*(2.0*la2-1.0)/tmp3;
		laplace_coeffs[JP1] =  2.0/lambda*lp1*(2.0*lambda-1.0)/tmp1;
		laplace_coeffs[JM1] = -2.0*pow(lambda,3)*(lambda-2.0)*lp1/tmp1;
		laplace_coeffs[JM2] =  2.0*pow(lambda,7)*(la2-2.0)/tmp3;
		laplace_coeffs[J] = -( laplace_coeffs[JM2]+laplace_coeffs[JM1]
		                      +laplace_coeffs[JP2]+laplace_coeffs[JP1] );

		if( lambda > 1.3 ){
			cout << "WARNING in cnSolver.h constructor:" << endl
				<< "	>> Stretch factor too large for stability of 3rd order stencil" << endl
				<< "		( try setting lambda < 1.3 or change STENCIL to 0 )" << endl;
		} // end lambda error warning
	} else {	// 2nd order, robustly stable
		laplace_coeffs[JP2] = -2.0*l5*(lambda-1.0)/(lp1*(1.0+l6)*tmp1);
		laplace_coeffs[JP1] = 2.0*lambda*(1+2.0*lambda-la2-2.0*l3+2.0*l5+l6)/(lp1*tmp2);
		laplace_coeffs[JM1] = 2.0*la2*(1+2.0*lambda-2.0*l3-l4+2.0*l5+l6)/(lp1*tmp2);
		laplace_coeffs[JM2] = -laplace_coeffs[JP2];
		laplace_coeffs[J  ] = -(laplace_coeffs[JP1] + laplace_coeffs[JM1]);
		
		if( 0 != STENCIL ){
			cout << "WARNING in cnSolver.h constructor:" << endl
				<< "	>> Stencil improperly specified, resorting to 2nd order" << endl
				<< "		(i.e. STENCIL = 0)" << endl;
		}// end error if
	} // end STENCIL if for gradient term	

	// for gradient term ...
	grad_coeffs[JP2] = -1.0/(lambda*lp1*(1.0+la2)*tmp1);
	grad_coeffs[JP1] = lp1/(lambda*tmp1);
	grad_coeffs[JM1] = -pow(lambda,3)*lp1/tmp1;
	grad_coeffs[JM2] = pow(lambda,7)/(lp1*(1.0+la2)*tmp1);
	grad_coeffs[J  ] = -( grad_coeffs[JM2] + grad_coeffs[JM1]
	                     +grad_coeffs[JP1] + grad_coeffs[JP2] );

	// for constant term ...
	const_coeffs[JP2] = 0.0;
	const_coeffs[JP1] = 0.0;
	const_coeffs[J  ] = 1.0;
	const_coeffs[JM1] = 0.0;
	const_coeffs[JM2] = 0.0;
	
}// end constructor 

/*
 *	STEP
 *
 *		Solves our equations for the next time step,
 *		using a Crank-Nicolson scheme and NR3 matrix solver
 *
 */
int cnSolver::step( problemDomain &domain,
                    gasDisk &disk,
                    secondaryBH &secondary ) 
{
	int status = EXIT_SUCCESS;
	size_t N = disk.N;
	double tmp0,tmp1,tmp2,lambda = disk.lambda;
	int offset;

	if( domain.debug_mode && domain.isWriteCycle() ){
		cout << endl << endl << "# -----------------------" 
			<< "----------------------------------" << endl
			<< "#l		D_J		tmp0		tmp1" << endl;
	} // end debug if

	/*
	 *	Build matrix M and solution vector d elements, based on finite 
	 *	diff coeffs, current specific angular momentum flux, diffusion
	 *	coeffs and torque density profile
	 */
	for( int j = 2 ; j < disk.N-2 ; j++ ){
	
		tmp0 = pow(lambda,-2.0*j)*.5*domain.dt/disk.dl2*disk.DJ[j]/(1.0-disk.nd);
		tmp1 = -pow(lambda,-1.0*j)*.5*domain.dt/disk.dl*disk.DJ[j]/(1.0-disk.nd);

		if( domain.debug_mode && domain.isWriteCycle() ){	
			cout << disk.l[j] << "	" << disk.DJ[j] << "	" << tmp0 << "	" << tmp1 << endl;
		}// end debug if

		d[j] = 0.0;
		for( int k = 0 ; k < STENCIL_SIZE ; ++k ){
			offset = j - CNTR + k;
			tmp2 = tmp1*secondary.torque(disk,offset)/disk.DJ[offset];

			M[j][k] = -tmp0*laplace_coeffs[k] - tmp2*grad_coeffs[k] + const_coeffs[k];
			d[j] +=  ( tmp0*laplace_coeffs[k] + tmp2*grad_coeffs[k] + const_coeffs[k] )*disk.Fj[offset];
		}// end k for

	} // end j for

	/*
	 * --------- Update boundary conditions ...
	 */
	double lp1 = lambda+1.0,la2=lambda*lambda,
		grad_const = 0.0, laplace_const = 0.0, laplace_val = 0.0;
	tmp1 = lambda*lambda+lambda+1.0, tmp2 = lambda*lambda-lambda-1.0;

	if( NEUMANN == disk.inner_bndry_type ){			// ####### INNER BOUNDS

		// Fixed gradient
		grad_const = 1.0/(disk.dl*lambda);

		M[0][R1] =  grad_const*lp1;
		M[0][R2] = -grad_const/lp1;
		M[0][C]  = -(M[0][R1]+M[0][R2]);
		d[0] = disk.inner_bndry_value;

		// Zero laplacian
		laplace_const = 2.0/(la2*disk.dl2*lp1);
		laplace_val = 0.0;
		if( disk.inner_bndry_laplacian == SELF_SIM && domain.t > 0.0){
			double x = disk.dl/sqrt(4.0*disk.D0*domain.t);
			laplace_val = (1.0-disk.inner_bndry_value)*sqrt(1.0/(PI*disk.DJ[1]*domain.t))*exp(-x*x);
		}

		M[1][L1] =  laplace_const*pow(lambda,3)*(lambda+2.0)/tmp1;
		M[1][R1] =  laplace_const*(la2+lambda-1.0)/lambda;
		M[1][R2] = -laplace_const*(lambda-1.0)/lambda/tmp1;
		M[1][C]  = -(M[1][L1]+M[1][R1]+M[1][R2]);
		d[1] = laplace_val;

	} else if( DIRICHLET == disk.inner_bndry_type ){

		// Constant Value
		M[0][R1] = 0.0;
		M[0][R2] = 0.0;
		M[0][C]  = 1.0;
		d[0] = disk.inner_bndry_value;

		// Zero laplacian
		M[1][L1] = pow(lambda,3)*(lambda+2.0)/tmp1;
		M[1][R1] = (la2+lambda-1.0)/lambda;
		M[1][R2] = -(lambda-1.0)/lambda/tmp1;
		M[1][C]  = -(M[1][L1]+M[1][R1]+M[1][R2]);
		d[1] = 0.0;

	} else {
		cerr << "ERROR --- Inner Bndry Type Improperly Specified as " 
			<< disk.inner_bndry_type << endl;
		return EXIT_FAILURE;
	} // end outer BC if/else

	if( NEUMANN == disk.outer_bndry_type ){			// ########### OUTER BOUNDS
	
		// Laplacian
		laplace_const = 2.0*lambda*pow(lambda,-2.0*(N-2.0))/lp1/disk.dl2;
		laplace_val = 0.0;

		M[N-2][L2] =  laplace_const*la2*la2*(lambda-1.0)/tmp1;
		M[N-2][L1] = -laplace_const*lambda*tmp2;
		M[N-2][R1] =  laplace_const*(2.0*lambda+1.0)/tmp1;
		M[N-2][C]  = -1.0*(M[N-2][L2]+M[N-2][L1]+M[N-2][R1]);
		d[N-2] = laplace_val;
 
		// Fixed gradient
		grad_const = pow(lambda,-(N-1.0))/tmp1/disk.dl;

		M[N-1][L2] = -grad_const*la2*la2/lp1;
		M[N-1][L1] =  grad_const*lp1;
		M[N-1][C]  = -(M[N-1][L2]+M[N-1][L1]);
		d[N-1] = disk.outer_bndry_value;
	
	} else if( DIRICHLET == disk.outer_bndry_type){

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
		d[N-1] = disk.outer_bndry_value;

	} else {
		cerr << "ERROR --- Outer Bndry Type Improperly Specified as " 
			<< disk.outer_bndry_type << endl;
		return EXIT_FAILURE;
	} // end outer BC if/else

	// Solve Matrix
	Bandec banded(M,2,2);
	banded.solve(d,FjNew);

	// Check for negative
	for( int j = 0 ; j < N ; j++ ){
		if( FjNew[j] < 0.0 ){
			if( disk.density_floor < 0.0 ){
				cerr << "ERROR IN CN SOLVER: Density negative @ j = " << j << endl
					<< "	>> t = " << domain.t << ", tStart = " << domain.tStart << ", dt= " << domain.dt << endl;
				status = EXIT_FAILURE;
			} else {
				FjNew[j] = disk.density_floor;  // if floor enabled
				cout << "WARNING IN CN SOLVER: Density negative @ j = " << j << endl 
					<< "	>> t = " << domain.t << ", tStart = " << domain.tStart << ", dt = " << domain.dt << endl
					<< "		Density Floor of " << disk.density_floor << " activated" << endl;
			} // end floor if/else
		}// end negative density if
	} // end j for
	
	// Copy new density into sigma
	for( int j = 0 ; j < N ; j++ )
		disk.Fj[j] = FjNew[j];
	
	return status;
} // end solve

/*
 *	MDOT
 *	
 *		Calculates the mass flux into a given radius:
 *
 *			M_dot = dF_J/dl - F * Lambda / D_J
 *
 */
double cnSolver::Mdot( const problemDomain &domain,
                       const gasDisk &disk, 
                       const secondaryBH &secondary, 
                       const int j ) const
{

	if( j < 2 || j > disk.N - 3 ){
//		cerr << "WARNING -- Cannot compute mass flux with 2 cells of bounds" << endl; FIXME
		return 0.0;
	} // end j if
	
	double grad_const = pow(disk.lambda,-1.0*j)/disk.dl, tmp=0.0;
	for( int k = 0 ; k < STENCIL_SIZE ; ++k )
		tmp += grad_const*grad_coeffs[k]*disk.Fj[ j - CNTR + k ];

	tmp -= disk.Fj[j]*secondary.torque(disk,j)/disk.DJ[j];
	
	return tmp;
} // end mDOt

#endif
