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

#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "readWrite.h"
#include "mr.h"
#include "cubicSpline.h"
#include "initialize.h"

int main(int argc , char **argv ){

	int status = EXIT_SUCCESS;

	// pull integration bounds from command line
	double a = 6.0, b = 15.0;
	if( argc > 2 ){
		sscanf(argv[1],"%lf",&a);
		sscanf(argv[2],"%lf",&b);
		if( a > b ) {
			cerr << "ERROR IN PROVE QUAD -- integration bounds out of order" << endl;
			return EXIT_FAILURE;
		} // end error if
	} // end argc if


	problemDomain domain;
	gasDisk disk;
	secondaryBH secondary;

	// read in from parameter file, params.in
	if( EXIT_SUCCESS != (status = readParams(domain,disk,secondary)))
		return status;

	// intialize grid
	if(EXIT_SUCCESS != (status = initialize(0,NULL,domain,disk,secondary)))
		return status;

	// hand-code constant value for Fj (remove it from problem)
	for( size_t j = 0 ; j != disk.N ; ++j ) disk.Fj[j] = 1.0;

	// set DJ to simply 1.0
	disk.nd = 0.0;
	disk.np = 0.0;
	disk.D0 = 1.0;

	// write globals to params.out
	writeParams(domain,disk,secondary);

	// Perform interpolation on data
	cubicSpline cSpline(disk.l,disk.Fj);
	
	// print interpolation to cout
	{
		int N = 1000;
		double lMin = 6.0, lMax = 17.0, dl = (lMax - lMin)/(N-1.0),l;
		for( int i = 0 ; i != N ; ++i ){
			l = lMin + i*dl;
			cout << l << "\t" << cSpline.interp(l) << "\t" << secondary.torque(disk,l,domain.M)<< endl;
		}// end i for
	}// end scope

	// perform quadrature
	secondary.moveSecondary(disk,domain.dt,domain.M);

	return status;
} // end main
