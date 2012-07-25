#include <iostream>
#include <cmath>
#include <fstream>

#include "mr.h"
#include "rk4.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

inline double max( double x, double y){ return ( x > y ? x : y ); }

double r_out  = 1.0E5,  // outer bnds
	r_isco = 6.0,         // ISCO
	hbar   = 0.1,         // dh/dr
	nv     = 0.5,         // visc const.
	q      = 0.08,         // binary mass ratio
	fconst = 0.01,        // calibration factor
	a      = 100.0;				// position of secondary

const size_t N = 1000;            // number of grid pts

double s1 = sqrt(a/r_out),
	s2 = sqrt(r_isco/r_out),
	s18 = pow(s1,8);

double nu(double r){ return pow(r/r_out,nv); }

double g = 2.0/3.0*fconst*q*q*sqrt(r_out);	//FIXME

double f(double s)
{
	static double g_star, Dp, vbar,s12,ss;

	g_star = g*s18;
	if( s < s1 )
		g_star = -g*pow(s,8.0);

	ss = s*s;
	s12 = s1*s1;
	Dp = max( fabs(ss - s12) , ss*hbar );
	vbar = pow(s,nv*2.0);

	return g_star*pow(Dp,-4.0)/vbar;
} // end f

void FnG( const double x , vDoub_i &F , vDoub_o &dFdx ){
	dFdx[0] = f(1.0-x);
	dFdx[1] = F[0];
} // end FnG

void FnH( const double x , vDoub_i &F, vDoub_o &dFdx ){
	dFdx[0] = f(1.0-x);
	dFdx[1] = exp( F[0] - F[1] );
}// end HnG

int main(int argc , char **argv ){

	double xStar = 0.5;
	if( argc > 1 ){
		sscanf(argv[1],"%lf",&xStar);
	} // end xStar if 

	cerr << "xStar = " << xStar << endl;

	/*
	 * setup LNS equation for F and G in terms and solve
	 * in terms of x = 1 - x, with ode45 (RK adaptive 5th 
	 * order)
	 */
	double xMin = 0.0,                  // s_start = 1.0
		xMax = 1.0 - sqrt(r_isco/r_out);  // s_end = s2 = sqrt(r_isco/r_out)
	size_t maxIters = 1E5, nVar = 2;    // some solver vars
	vDoub FInit(nVar,0.0);              // ICs
	mDoub FOut_FnG(maxIters,nVar+1);    // output from ode45 solver
	mDoub FOut_FnH(maxIters,nVar+1);    // output from ode45 solver

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
	 * Print r and sigma
	 */
	for( size_t i = 0 ; i != r.size() ; ++i )
		cout << r[i] << "\t" << sigma[i] << endl;

	return EXIT_SUCCESS;
} // end main
