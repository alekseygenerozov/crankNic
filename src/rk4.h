#include <iostream>
#include "mr.h"

using std::cout;
using std::cerr;
using std::endl;

#ifndef INC_RK4
#define INC_RK4

/*
 *		RK4
 *
 *			A standard 4-step runge kutta stepper
 */
void rk4(	vDoub_i &f,                 // current values
          vDoub_i &dfdx,              // current derivatives
          const double x,             // current x pos
          const double h,             // step size
          vDoub_o &fnew,              // solution vector
          void deriv( const double,   // ODE fcn
                      vDoub_i &, 
                      vDoub_o &))      
{
	size_t n = f.size();
	vDoub dfm(n), dft(n), ft(n);
	double hh = h*0.5,
		h6 = h/6.0,
		xh = x+hh;
	
	for( size_t i = 0 ; i != n ; ++i )	// 1st 1/2 guess
		ft[i] = f[i] + hh*dfdx[i];	
	deriv(xh,ft,dft);										// derives @ 1st 1/2 guess

	for( size_t i = 0 ; i != n ; ++i )	// 2nd 1/2 guess
		ft[i] = f[i] + hh*dft[i];
	deriv(xh,ft,dfm);										// derivs @ 2nd 1/2 guess

	for( size_t i = 0 ; i != n ; ++i ){	// guess
		ft[i] = f[i] + h*dfm[i];
		dfm[i] += dft[i];	// accumulate 1/2 guess's
	}// end i for
	deriv(x+h,ft,dft);                  // derivs @ guess

	for( size_t i = 0 ; i != n ; ++i )
		fnew[i] = f[i] + h6*(dfdx[i] + dft[i] + 2.0*dfm[i] );
}// end rk4

/*
 *  ODE 45
 *
 *  Implements RK4 with adaptive stepsize and a 
 *  fifth order correction term.
 */
size_t ode45( vDoub_i &f_init,             // initial values
              mDoub_o &fOut,               // solution
              const double xMin,           // starting x position
              const double xMax,           // final x position
              double h,                    // initial stepsize
              void deriv( const double,    // fcn of interest
                          vDoub_i &, vDoub_o &),
              const double rel_tol=1E-7,   // error tolerance
              const double abs_tol=1E-11 )
{
  size_t nvar = f_init.size();
	vDoub fCurr = f_init, fFull(nvar),fHalf(nvar),fHalf2(nvar),
		dfdx(nvar),dfdxHalf(nvar);
	double x=xMin,hh=0.5*h,err=0.0,tmp=0;

	// copy intial f_init into fOut
	fOut[0][0] = xMin;
	for(size_t j=0;j!=nvar;++j) fOut[0][j+1] = f_init[j];

	// solve ...
	int nIters = 1;
	for( size_t i = 1 ; nIters < fOut.nrows() ; ++i ){

		// check we're not over-shooting
		if( ( x + h ) > xMax )
			h = (xMax - x);
		
		// find current derivative
		deriv(x,fCurr,dfdx);

		// take a full step
		rk4(fCurr,dfdx,x,h,fFull,deriv);

		// take two half steps
		hh = 0.5*h;
		rk4(fCurr,dfdx,x,hh,fHalf,deriv);
		deriv(x+hh,fHalf,dfdxHalf);
		rk4(fHalf,dfdxHalf,x+hh,hh,fHalf2,deriv);

		// Calculate error
		err = 0.0;
		for( size_t j = 0 ; j != nvar ; ++j ){
			tmp = (fHalf2[j] - fFull[j])/( abs_tol + rel_tol*fHalf2[j]);
			err += tmp*tmp;
		}// end j for
		err /= nvar;

		// use error to accept / reject step, resize
		if( err < 1.0 ){

			x += h;

			for( size_t j = 0 ; j != nvar ; ++j )
				fCurr[j] = fHalf2[j] + (fHalf2[j] - fFull[j])/15.0;
			h = .9 * h * pow( err , -1.0/5.0 );

			// copy result into fOut
			fOut[nIters][0] = x;
			for(size_t j=0;j!=nvar;++j) fOut[nIters][j+1] = fCurr[j];
		
			++nIters;
	
		} else {
			h /= 2.0;
		}

		if( x >= xMax )
			return nIters;
	}// end i for
	return 0;
}// end ode45

#endif
