#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "mr.h"

#ifndef FD_SOLVER
#define FD_SOLVER

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

	double laplace_coeffs[STENCIL_SIZE];  // ""                           2nd deriv
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
	double tmp0,lambda=disk.lambda,dFoDdt,FoD;

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
			dFoDdt = tmp0*laplace_coeffs[k]*disk.Fj[j-CNTR+k];

		// upwind differencing for the advective term
		if( disk.l[j] < secondary.l_a ){ 
			dFoDdt -= (   disk.Fj[j+1]*secondary.torque(disk,l[j+1],domain.M)/disk.Dj(j+1)
			            - disk.Fj[j  ]*secondary.torque(disk,l[j  ],domain.M)/disk.Dj(j  )
			          )/( disk.l[j+1]-disk.l[j]);
		} else {
			dFoDdt -= (   disk.Fj[j  ]*secondary.torque(disk,l[j  ],domain.M)/disk.Dj(j  )
			            - disk.Fj[j-1]*secondary.torque(disk,l[j-1],domain.M)/disk.Dj(j-1)
			          )/( disk.l[j]-disk.l[j-1]);
		} // end upwind if/else

		// find new Fj/Dj
		FoD = disk.Fj[j]/disk.Dj[j] + domain.dt*dFoDdt;

		// Find Fj via Dj
		if( disk.nd == 0 )
			disk.Fj[j] = FoD*disk.D0*pow(disk.l[j],disk.np);
		else 
			disk.Fj[j] = pow(FoD*disk.D0*pow(disk.l[j],disk.np),1.0/(1.0-nd));
	} // end j for

	/*
	 * --------- Update boundary conditions ...
	 */
	double lp1 = lambda+1.0,la2=lambda*lambda,
		grad_const = 0.0, laplace_const = 0.0, laplace_val = 0.0;
	tmp1 = lambda*lambda+lambda+1.0, tmp2 = lambda*lambda-lambda-1.0;
	double r2,r1,c,l1,l2;

	// --------------- INNER BOUNDARY

	// Fixed Laplacian
	laplace_const = 2.0/(la2*disk.dl2*lp1);
	laplace_val = 0.0;
	if( disk.inner_bndry_laplacian == SELF_SIM && domain.t > 0.0){
		double x = disk.dl/sqrt(4.0*disk.D0*domain.t);
		laplace_val = (1.0-disk.inner_bndry_value)*sqrt(1.0/(PI*disk.Dj(1)*domain.t))*exp(-x*x);
	}

	l1 =  laplace_const*pow(lambda,3)*(lambda+2.0)/tmp1;
	r1 =  laplace_const*(la2+lambda-1.0)/lambda;
	r2 = -laplace_const*(lambda-1.0)/lambda/tmp1;
	c  = -(l1+r1+r2);
	disk.Fj[1] = (laplace_val-r2*disk.Fj[3]-r1*disk.Fj[2]-l1*disk.Fj[0])/c;

	if( NEUMANN == disk.inner_bndry_type ){			// Fixed Gradient
		grad_const = 1.0/(disk.dl*lambda);
		r2 = -grad_const/lp1;
		r1 =  grad_const*lp1;
		c  = -(g2+g1);
		disk.Fj[0] = (disk.inner_bndry_value-r2*disk.Fj[2]-r1*disk.Fj[1])/c;
	} else if( DIRICHLET == disk.inner_bndry_type ){ // Fixed Value
		disk.Fj[0] = disk.inner_bndry_value;
	} else { // none specified
		cerr << "ERROR --- Inner Bndry Type Improperly Specified as " 
			<< disk.inner_bndry_type << endl;
		return EXIT_FAILURE;
	} // end outer BC if/else
	
	// --------------- OUTER BOUNDARY

	// Laplacian
	laplace_const = 2.0*lambda*pow(lambda,-2.0*(N-2.0))/lp1/disk.dl2;
	laplace_val = 0.0;
	l2 =  laplace_const*la2*la2*(lambda-1.0)/tmp1;
	l1 = -laplace_const*lambda*tmp2;
	r1 =  laplace_const*(2.0*lambda+1.0)/tmp1;
	c  = -(l2+l1+r1);
	disk.Fj[N-2] = (laplace_val-l2*disk.Fj[N-4]-l1*disk.Fj[N-3]-r1*disk.Fj[N-1])/c;
 
	if( NEUMANN == disk.outer_bndry_type ){			// fixed gradient
		grad_const = pow(lambda,-(N-1.0))/tmp1/disk.dl;
		l2 = -grad_const*la2*la2/lp1;
		l1 =  grad_const*lp1;
		c  = -(M[N-1][L2]+M[N-1][L1]);
		disk.Fj[N-1] = (disk.outer_bndry_value-l2*disk.Fj[N-3]-l1*disk.Fj[N-2])/c;
	} else if( DIRICHLET == disk.outer_bndry_type){  // Constant Value
		disk.Fj[N-1] = disk.outer_bndry_value;
	} else { // none specified
		cerr << "ERROR --- Outer Bndry Type Improperly Specified as " 
			<< disk.outer_bndry_type << endl;
		return EXIT_FAILURE;
	} // end outer BC if/else

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
