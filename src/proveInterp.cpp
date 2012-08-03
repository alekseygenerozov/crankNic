/*
 *	PROVE INTERP
 *
 *		Small program to check out performance of interpolator 
 *		for the solver in any state. Could come in handy later,
 *		haha. To run:
 *
 *		g++ proveInterp.cpp -o Interp
 *		./Interp -r outputfiles/T020.dat > interp.dat
 *
 *		Where you can replace "outputfiles/T020.dat" with any dat file
 */

#include "global.h"
#include "readWrite.h"
#include "mr.h"
#include "cubicSpline.h"
#include "initialize.h"

int main(int argc , char **argv ){

	int status = EXIT_SUCCESS;

	/*
	 * Read in parameter file and update global variables
	 */

	// read in from parameter file, params.in
	if( EXIT_SUCCESS != (status = readParams()))
		return status;

	// intialize grid
	double l[N], FJ[N];            // not used until after ODE solver
	int fCnt; double t;
	if(EXIT_SUCCESS != (status = initialize(argc,argv,l,FJ,fCnt,t)))
		return status;

	// write globals to params.out
	writeParams();

	/*
	 * Perform interpolation on data and print results to stdout
	 */
	vDoub xx(N,l), ff(N,FJ);
	cubicSpline cSpline(xx,ff);

	/*
	 * Print result, at higher res ...
	 */
	int Ni = 3000;
	double xMin = sqrt(6)+.1, xMax = 30.0,
		dx = (xMax - xMin) / ( Ni - 1.0 ), x;
	for( size_t j = 0 ; j != Ni ; ++j ){
		x = xMin + j * dx;
		cout << x << "\t" << cSpline.interp(x) << endl;
	}

	return status;
} // end main
