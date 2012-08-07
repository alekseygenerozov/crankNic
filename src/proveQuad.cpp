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

#define _QUAD_CHECK_ 	// switches torque to cleaner, square fcn

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
	double a = 7.0, b = 16.0;
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

	// hand-code function for Fj:
	{
		double ll;
		for( size_t j = 0 ; j != disk.N ; ++j ){
			ll = disk.l[j];
			disk.Fj[j] = (ll-8.0)*(ll-10.0)*(ll-12.0)*(ll-14.0) + 50.0;
		} // end j for
	}// end ll scope

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
	cerr << secondary.gaussTorqueInt(cSpline,disk,a,b,domain.M) << endl;

	return status;
} // end main
