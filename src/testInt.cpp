#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "disk.h"
#include "initialize.h"
#include "integrator.h"
#include "readWrite.h"
#include "torques.h"

int main(){

	int status = EXIT_SUCCESS;

	// read in from parameter file, params.in
	if( EXIT_SUCCESS != (status = readParams()))
		return status;

	// Create arrays for data
	double r[N];				// radial position
	double ff[N];				// test function
	double FF;					// integral


	// Intialize r and ff (begign for ff)
	if(EXIT_SUCCESS != (status = initialize(r,ff)))
		return status;

	// print Parameters we'll use
	if(EXIT_SUCCESS != (status = writeParams()))
		return status;

/*
	// ---------------------------------------- integral tests ...
	double real,perDiff;

  // test function f(r) = r
  for( int j = 0 ; j < N ; j++ )
    ff[j] = r[j];               // f(r) = r
	FF = simpsonInt(ff);
	real = .5*(rMax*rMax - rMin*rMin);
	perDiff = (FF-real)/(real)*100.0;
	fprintf(stderr,"f(x) = x\n\t>>>     percent diff: %g\n",perDiff);

	// test function f(r) = -(r-1)(r-3)
  for( int j = 0 ; j < N ; j++ )
    ff[j] = -(r[j]-1.0)*(r[j]-3.0);
  FF = simpsonInt(ff);
  real =  -1.0/3.0*(pow(rMax,3)-pow(rMin,3))
          +2.0/1.0*(pow(rMax,2)-pow(rMin,2))
          -3.0/1.0*(pow(rMax,1)-pow(rMin,1));
  perDiff = (FF-real)/(real)*100.0;
  fprintf(stderr,"f(x) = -(x-1)(x-3)\n\t>>>     percent diff: %g\n",perDiff);
	
	// test function f(r) = -r(r-1)(r-2)(r-3)
	for( int j = 0 ; j < N ; j++ )
		ff[j] = -1.0*r[j]*(r[j]-1.0)*(r[j]-2.0)*(r[j]-3.0);
	FF = simpsonInt(ff);
	real =  -1.0/5.0*(pow(rMax,5)-pow(rMin,5))
					+3.0/2.0*(pow(rMax,4)-pow(rMin,4))
					-11./3.0*(pow(rMax,3)-pow(rMin,3))
					+3.00000*(pow(rMax,2)-pow(rMin,2));
	perDiff = (FF-real)/(real)*100.0;
	fprintf(stderr,"f(x) = -x(x-1)(x-2)(x-3)\n\t>>>     percent diff: %g\n",perDiff);

	// test torque profile:
	for( int j = 0 ; j < N ; j++ )
		ff[j] = tidalTorque(r[j],a,h(r[j]));
	FF = simpsonInt(ff);
	fprintf(stderr,"TORQUE PROFILE:\n\t>>>		%g\n",FF);
*/

/*
	// Test new Bode for torque profile
	for( int j=0;j<N;j++){
		ff[j] = 1.0;
	}
	FF = bodeInterpIntWithTorque(r,ff);
	fprintf(stderr,"TORQUE PROFILE:\n\t>>>		%g\n",FF);
*/

	// Test Gaussian Quadrature version:
//	fprintf(stderr,"\t\t>> Result of 16-point Gauss Quadrature scheme: %g\n",gaussIntWithTorque(16));	
//	fprintf(stderr,"\t\t>> Result of 48-point Gauss Quadrature scheme: %g\n",gaussIntWithTorque(48));	
//	fprintf(stderr,"\t\t>> Result of 96-point Gauss Quadrature scheme: %g\n",gaussIntWithTorque(96));	
//	fprintf(stderr,"\t\t>> Result of 100-point Gauss Quadrature scheme: %g\n",gaussIntWithTorque(100));	

	// test speed of this thing:
	for( int i = 0 ; i < 1000000 ; i++ )
		gaussIntWithTorque(96);

  // print last file:
  if(EXIT_SUCCESS != (status = writeStandard(-1,r,ff,a,0.0)))
    return status;

	return status;
} // end main
