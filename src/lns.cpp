/*
 * 	LNS
 *
 *		July 2012
 *
 * 		Finds steady state solution from Liu & Schapiro 2010 for the case
 *		of a binary black hole who's orbital separation remains fixed.
 *
 *		The solutions is more-or-less analytic, but involves solving a pair
 *		of coupled 1st order ODEs. We use an RK4/5 solver to accomplish this.
 *
 *		The quantity rStar determines where in the integration we switch
 *		between two different formulations of the ODEs ... each poorly behaved
 *		in certain regimes. It can be fed in via command line at runtime.
 *
 *		The integration works from the OUTER disk bounds inward.
 */
#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "readWrite.h"
#include "mr.h"
#include "rk4.h"
#include "initialize.h"

// see main for explanation of variables
double r_out,r_isco,nv,a,s1,s2,s18,g,dhdr;
double nu(double r){ return pow(r/r_out,nv); } // viscosity

/*
 *	ff 
 *		Dimensionless torque fcn from L&S 2010
 */
double ff(double s)
{
	static double g_star, Dp, vbar,s12,ss;

	g_star = g*s18;
	if( s < s1 )
		g_star = -g*pow(s,8.0);

	ss = s*s;
	s12 = s1*s1;
	Dp = max( fabs(ss - s12) , ss*dhdr );
	vbar = pow(s,nv*2.0);

	return g_star*pow(Dp,-4.0)/vbar;
} // end f


/*
 * 	FNG and FNH
 *
 *		1st order coupled ODEs from LNS 2010 we solve for analytic 
 *	steady state
 */
void FnG( const double x , vDoub_i &F , vDoub_o &dFdx ){
	dFdx[0] = ff(1.0-x);
	dFdx[1] = F[0];
} // end FnG

void FnH( const double x , vDoub_i &F, vDoub_o &dFdx ){
	dFdx[0] = ff(1.0-x);
	dFdx[1] = exp( F[0] - F[1] );
}// end HnG

int main(int argc , char **argv ){

	int status = EXIT_SUCCESS;

	problemDomain domain;
	gasDisk disk;
	secondaryBH secondary;
  
	// read in from parameter file, params.in
  if( EXIT_SUCCESS != (status = readParams(domain,disk,secondary)))
    return status;

	domain.problemType = RAMPED; // just needs to not look for some file

	// intialize grid
  if(EXIT_SUCCESS != (status = initialize(0,NULL,domain,disk,secondary)))
    return status;

	r_out  = 300.0*300.0; // FIXME lMax*lMax/M;    // outer bounds
	r_isco = disk.lMin*disk.lMin/domain.M;         // ISCO (inner bounds)
	nv     = .5*disk.np+1;                         // viscosity power law const
	a      = secondary.l_a*secondary.l_a/domain.M; // position of secondary BH
	s1 = sqrt(a/r_out);                            // LNS dimensionless radius-esque thing
	s2 = sqrt(r_isco/r_out);                       // """, at the inner bnds
	s18 = pow(s1,8);                
	g = 2.0/3.0*secondary.f*secondary.q*secondary.q
	      *sqrt(r_out);                            // measure of secondary strength (~q^2)
	dhdr = disk.dhdr;

	// xStar determines where switch ODE formulation (matters, a lot!)
	double xStar = 0.5;
	if( argc > 1 ){
		sscanf(argv[1],"%lf",&xStar);
	} else {
		xStar = 0.5*(1.0 + sqrt(a/r_out));
	}// end xStar if 
	cerr << "xStar = " << xStar << endl;

	/*
	 * setup LNS equation for F and G in terms and solve
	 * in terms of x = 1 - x, with ode45 (RK adaptive 5th 
	 * order)
	 */
	double xMin = 0.0,                        // s_start = 1.0
		xMax = 1.0 - sqrt(r_isco/r_out);        // s_end = s2 = sqrt(r_isco/r_out)
	size_t maxIters = (size_t)1E5, nVar = 2;  // some solver vars
	vDoub FInit(nVar,0.0);                    // ICs
	mDoub FOut_FnG(maxIters,nVar+1);          // output from ode45 solver
	mDoub FOut_FnH(maxIters,nVar+1);          // output from ode45 solver

	double relTol = 1E-7;
	double absTol = 1E-11;

	// integrate from s_start -> s_star
	size_t nIter_FnG = 0;
	if( 0 == ( nIter_FnG = ode45(FInit,FOut_FnG,xMin,xStar,0.1,FnG,relTol,absTol))){
	  cerr << "ERROR IN y ... ODE45 failed to integrate" << endl;
	  return EXIT_FAILURE;
	} // end err if

	// next integrate from s_star -> s2
	size_t nIter_FnH = 0;
	FInit[0] = FOut_FnG[nIter_FnG-1][1];
	FInit[1] = log(FOut_FnG[nIter_FnG-1][2]);
	if( 0 == ( nIter_FnH = ode45(FInit,FOut_FnH,xStar,xMax,0.1,FnH,relTol,absTol))){
		cerr << "ERROR IN y ... ODE45 failed to integrate" << endl;
		return EXIT_FAILURE;
	}

	/*
	 * transform to sigma(r) ...
	 */
	size_t totalIters = nIter_FnG + nIter_FnH;
	vDoub r(totalIters);
	vDoub sigma(totalIters);
	double mdot = exp(-FOut_FnH[nIter_FnH-1][2]),s=0;		// mdot = exp{-H(s2)}

	// solution for F and H region (above rStar)
	for( size_t i = 0, j = totalIters - 1 ; i != nIter_FnG ; ++i, --j ){
		s = 1.0 - FOut_FnG[i][0];
		r[j] = s*s*r_out;
		sigma[j] = exp(-FOut_FnG[i][1])*(1.0-FOut_FnG[i][2]*mdot)/s/nu(r[j]);
	}// end FnH loop

	// solution for F and H region (below rStar)
	for( size_t i = 0, j = nIter_FnH-1 ; i != nIter_FnH ; ++i , --j ){
		s = 1.0 - FOut_FnH[i][0];
		r[j] = s*s*r_out;
		sigma[j] = exp(-FOut_FnH[i][1])*(1.0-mdot*exp(FOut_FnH[i][2]))/s/nu(r[j]);
	} // end FnG loop



	/*
   *		Convert solution from r/s -> l/FJ, interpolate
   *		onto our simulation grid and print result
   */	

	// convert r and sigma into l and FJ
	double tmp,l_tmp;
	for( size_t i = 0 ; i != r.size() ; ++i ){
		l_tmp = sqrt(r[i]*domain.M);
		tmp = r[i]*omega_k(l_tmp,domain.M);
		sigma[i] = 3.0*PI*nu(r[i])*r[i]*r[i]*omega_k(l_tmp,domain.M)*sigma[i];
		r[i] = l_tmp;
	}

	// interpolate solution onto grid
	disk.Fj[0  ] = sigma[0  ]; disk.Fj[disk.N-1] = sigma[totalIters - 1];
	size_t j = 0; double m;
	for( size_t i = 1 ; i != disk.N-1 ; ++i ){ // at every grid pt
		while( r[j] < disk.l[i] ) ++j; // advance rk sltn until we bracket grid pt
		m = (sigma[j]-sigma[j-1])/(r[j]-r[j-1]);
		disk.Fj[i] = (disk.l[i]-r[j])*m + sigma[j];
	}

	// print solution
	for( size_t i = 0 ; i != disk.N ; ++i )
		cout << disk.l[i] << "\t" << disk.Fj[i] << endl;

	return status;
} // end main
