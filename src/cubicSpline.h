/*
 *	CUBIC SPLINE
 *
 *		3rd order spline interpolation: given vectors x and f
 *		of equal length, this routine calculates 2nd derivs
 * 		necessary for spline interpolation upon instantiation.
 *		
 *		Call interp (with an x value) to 
 */
#include "mr.h"

#ifndef INC_CUBESPLINE
#define INC_CUBESPLINE

class cubicSpline
{
public:
	cubicSpline(vDoub_i &x,vDoub_i &f);   // constructor
	double interp(double x);              // evaluate at x
private:
	size_t jLast;
	vDoub_i &xx, &ff;
	const size_t n;
	vDoub f2;
}; // end cubic spline



/*
 *	CONSTRUCTOR 
 *
 *		Initializes vector to store 2nd derivs
 *		at grid points and then fills it
 */
cubicSpline::cubicSpline(vDoub_i &x,vDoub_i &f)
	: xx(x),ff(f),n(x.size()),f2(x.size()),jLast(0)
{
	double p,sig;
	vDoub u(n-1);

	// tri-diag decomp y2 and u used as tmp storage
	f2[0] = u[0] = 0.0;
	for(size_t i=1 ; i < n-1 ; ++i ){
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*f2[i-1]+2.0;
		f2[i] = (sig-1.0)/p;
		u[i] = (f[i+1]-f[i])/(x[i+1]-x[i]) - (f[i]-f[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}// end i for	
	f2[n-1] = 0.0;

	// back substitution to obtain 2nd derivs, stored in f2
	for(ptrdiff_t i=n-2;i>=0;i--)  f2[i] = f2[i]*f2[i+1] + u[i];

}// end constructor for cubic spline



/*
 *	INTERP
 *
 *		Given x value, finds the interpolated
 *		function value there. Stores index point
 *		of original vector (xx and ff) where the last
 *		call was made to speed up this process
 */
double cubicSpline::interp(double x) 
{
#ifdef _CHECKBNDS_
	if( x < xx[0] || x > xx[n-1] )
		throw("x beyond vector domain in cubic spline interp");
#endif

	// find grid points bracketing x
	size_t jl = jLast, jh;
	while( xx[jl]   > x && jl   > 0    ) --jl;
	while( xx[jl+1] < x && jl+1 <  n-1 ) ++jl;
	jLast = jl; // update memory of last eval pt
	jh = jl+1;  // grid point above x
	if( xx[jl] == x ) return ff[jl];

	// some coeffs we need to find
	double f,h,b,a;
	h = xx[jh] - xx[jl];
	if(h==0.0) throw("Bad input to cubic spline interp");
	a = (xx[jh]-x)/h;
	b = (x-xx[jl])/h;
	
	// evaluate and return cubic polynomial
	return a*ff[jl]+b*ff[jh]
		+ ((a*a*a-a)*f2[jl]+(b*b*b-b)*f2[jh])*(h*h)/6.0;
}// end cubic spline interp

#endif
