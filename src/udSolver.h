#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "mr.h"

#ifndef UD_SOLVER
#define UD_SOLVER

class udSolver{
public:
	udSolver(const gasDisk& disk );
	int step(problemDomain &domain, gasDisk &disk, secondaryBH &secondary);
	double Mdot(const problemDomain &domain, const gasDisk &disk, 
		const secondaryBH &secondary, const int j ) const;
private:
	static const unsigned int STENCIL_SIZE = 5;              // # of cells in stencil (per time step)
	static const unsigned int CNTR = 2;                      // center of stencil (jth grid cell)
	static const unsigned int JP2=4,JP1=3,J=2,JM1=1,JM2=0;   // indexing for stencil
	static const unsigned int A=0,B=1,C=2;                   // for bndry_laplace indexing

	double laplace_coeffs[STENCIL_SIZE];  // finite diff coefficients for 2nd deriv
	double bndry_laplace[3];							// ""      near bndry       "" ""
	double grad_coeffs[STENCIL_SIZE];     // ""                           1st deriv
};

// Default Constructor
udSolver::udSolver(const gasDisk& disk)
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

	// for boundary laplace:
	bndry_laplace[A] = 2.0*la2/lp1;
	bndry_laplace[C] = 2.0*lambda/lp1;
	bndry_laplace[B] = -(bndry_laplace[A]+bndry_laplace[C]);

	// for gradient term ...
	grad_coeffs[JP2] = -1.0/(lambda*lp1*(1.0+la2)*tmp1);
	grad_coeffs[JP1] = lp1/(lambda*tmp1);
	grad_coeffs[JM1] = -pow(lambda,3)*lp1/tmp1;
	grad_coeffs[JM2] = pow(lambda,7)/(lp1*(1.0+la2)*tmp1);
	grad_coeffs[J  ] = -( grad_coeffs[JM2] + grad_coeffs[JM1]
	                     +grad_coeffs[JP1] + grad_coeffs[JP2] );

}// end constructor 

/*
 *	STEP
 *
 *		Solves our equations for the next time step,
 *		using an explicit upwind differencing scheme
 */
int udSolver::step( problemDomain &domain,
                    gasDisk &disk,
                    secondaryBH &secondary ) 
{
	int status = EXIT_SUCCESS;
	size_t N = disk.N;
	double tmp0,lambda=disk.lambda,dFoDdt,FoD,l2=lambda*lambda;

	if( domain.debug_mode && domain.isWriteCycle() ){
//		cout << endl << endl << "# -----------------------" 
//			<< "----------------------------------" << endl
//			<< "#l		D_J		tmp0		tmp1" << endl;
	} // end debug if

	for( int j = 2 ; j < disk.N-2 ; j++ ){

		dFoDdt = 0.0;
		tmp0 = pow(lambda,-2.0*j)/disk.dl2;

		if( domain.debug_mode && domain.isWriteCycle() ){	
//			cout << disk.l[j] << "	" << disk.Dj(j) << "	" << tmp0 << "	" << tmp1 << endl;
		}// end debug if

		// diffusion term
		for( int k = 0 ; k < STENCIL_SIZE ; ++k )
			dFoDdt += tmp0*laplace_coeffs[k]*disk.Fj[j-CNTR+k];

		// upwind differencing for the advective term
		double old = dFoDdt;
		if( disk.l[j] < secondary.l_a ){ 
			dFoDdt -= (   disk.Fj[j+1]*secondary.torque(disk,disk.l[j+1],domain.M)/disk.Dj(j+1)
			            - disk.Fj[j  ]*secondary.torque(disk,disk.l[j  ],domain.M)/disk.Dj(j  )
			          )/( disk.l[j+1]-disk.l[j]);
			cout << "change @ " << j << ": " << dFoDdt - old << endl;
		} else {
			dFoDdt -= (   disk.Fj[j  ]*secondary.torque(disk,disk.l[j  ],domain.M)/disk.Dj(j  )
			            - disk.Fj[j-1]*secondary.torque(disk,disk.l[j-1],domain.M)/disk.Dj(j-1)
			          )/( disk.l[j]-disk.l[j-1]);
			cout << "change @ " << j << ": " << dFoDdt - old << endl;
		} // end upwind if/else

		// find new Fj/Dj
		FoD = disk.Fj[j]/disk.Dj(j) + domain.dt*dFoDdt;

		// Find Fj via Dj
		if( disk.nd == 0 )
			disk.Fj[j] = FoD*disk.D0*pow(disk.l[j],disk.np);
		else 
			disk.Fj[j] = pow(FoD*disk.D0*pow(disk.l[j],disk.np),1.0/(1.0-disk.nd));
	} // end j for

	/*
	 * --------- Update boundary conditions:
	 */
	double laplace_val,tmp;

	if( NEUMANN == disk.inner_bndry_type ){		// INNER BNDRY SWITCH
		
		// fixed laplacian
		laplace_val = 0.0;
		if( disk.inner_bndry_laplacian == SELF_SIM && domain.t > 0.0){	// for Rafikov, 2012 solt'n
			double x = disk.dl/sqrt(4.0*disk.D0*domain.t);
			laplace_val = (1.0-disk.inner_bndry_value)*sqrt(1.0/(PI*disk.Dj(1)*domain.t))*exp(-x*x);
		}
		disk.Fj[1] = (   laplace_val*disk.dl2*l2
		               + bndry_laplace[A]*disk.inner_bndry_value*disk.dl
		               - bndry_laplace[C]*disk.Fj[2] )/( bndry_laplace[A] + bndry_laplace[C] );

		// fixed gradient
		disk.Fj[0] = disk.Fj[1] - disk.inner_bndry_value*disk.dl;

	} else if( DIRICHLET == disk.inner_bndry_type ){

		// fixed laplacian
		laplace_val = 0.0;
		if( disk.inner_bndry_laplacian == SELF_SIM && domain.t > 0.0){  // for Rafikov, 2012 solt'n
			double x = disk.dl/sqrt(4.0*disk.D0*domain.t);
			laplace_val = (1.0-disk.inner_bndry_value)*sqrt(1.0/(PI*disk.Dj(1)*domain.t))*exp(-x*x);
		}
		disk.Fj[1] = (   laplace_val*disk.dl2*l2
		               - bndry_laplace[A]*disk.Fj[0]
		               - bndry_laplace[C]*disk.Fj[2] )/bndry_laplace[B];

	} else {
		cerr << "ERROR --- Inner Bndry Type Improperly Specified as " 
		     << disk.inner_bndry_type << endl;
		return EXIT_FAILURE;
	} // end inner BC if/else
	
	if( NEUMANN == disk.outer_bndry_type ){		// OUTER BNDRY SWITCH

		tmp = disk.dl*pow(lambda,N-2);
		laplace_val = 0.0;
	
		// fixed laplacian
		disk.Fj[N-2] = (   tmp*(tmp*laplace_val - bndry_laplace[C]*disk.outer_bndry_value)
		                 - bndry_laplace[A]*disk.Fj[N-3] )/( bndry_laplace[B] + bndry_laplace[C] );

		// fixed gradient
		disk.Fj[N-1] = disk.Fj[N-2] + tmp*disk.outer_bndry_value;
		
	} else if( DIRICHLET == disk.outer_bndry_type ){
		
		// fixed laplacian
		tmp = disk.dl*pow(lambda,N-2);
		laplace_val = 0.0;
		disk.Fj[N-2] = (   disk.dl2*pow(lambda,2*N-4)*laplace_val - bndry_laplace[A]*disk.Fj[N-3]
		                 - bndry_laplace[C]*disk.Fj[N-1] )/bndry_laplace[B]; 

	} else {
		cerr << "ERROR --- Outer Bndry Type Improperly Specified as " 
			<< disk.outer_bndry_type << endl;
		return EXIT_FAILURE;
	} // end outer BC if/els

	// Check for negative
	for( int j = 0 ; j < N ; j++ ){
		if( disk.Fj[j] < 0.0 ){
			if( disk.density_floor < 0.0 ){
				cerr << "ERROR IN CN SOLVER: Density negative @ j = " << j << endl
					<< "	>> t = " << domain.t << ", tStart = " << domain.tStart << ", dt= " << domain.dt << endl;
				status = EXIT_FAILURE;
			} else {
				disk.Fj[j] = disk.density_floor;  // if floor enabled
				cout << "WARNING IN CN SOLVER: Density negative @ j = " << j << endl 
					<< "	>> t = " << domain.t << ", tStart = " << domain.tStart << ", dt = " << domain.dt << endl
					<< "		Density Floor of " << disk.density_floor << " activated" << endl;
			} // end floor if/else
		}// end negative density if
	} // end j for
	
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
double udSolver::Mdot( const problemDomain &domain,
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

	tmp -= disk.Fj[j]*secondary.torque(disk,disk.l[j],domain.M)/disk.Dj(j);
	
	return tmp;
} // end mDOt

#endif