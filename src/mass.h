#include "disk.h"
#include "cnSolver.h"

#ifndef INC_MASS
#define INC_MASS

/*
 *	MASS INTEGRAL
 *
 *		Integrates total disk mass, using Simpson's
 *		1/3rd rule on a logarithmic grid
 *
 */
double massIntegral(const double *l, const double *Fj, const int jMin = 0 , const int jMax = N-1 ){

	double result = 0.0,
		evenCoeff = (2.0*lambda-1.0)/(6.0*lambda*lambda*lambda) + (2.0-lambda)/6.0,
		oddCoeff  = (lambda+1.0)*(lambda+1.0)/(6.0*lambda*lambda);

	result += (2.0-lambda)/6.0*dmdl(Fj[jMin],l[jMin]);
	for(int j=jMin+1;j<jMax;j++)
		result += (((j-jMin)%2==0)?evenCoeff:oddCoeff)*dmdl(Fj[j],l[j])*pow(lambda,j);
	result += (2.0*lambda-1.0)/6.0*pow(lambda,N-4)*dmdl(Fj[jMax],l[jMax]);

	return result*dl*(lambda+1.0);
}// end massIntegral

/*
 *	MASS ANNULUS
 *
 *		Returns the mass contained within two radii 
 *
 */
bool massAnnulus( 
			const double *l, 
			const double *Fj, 
			const double alMin, 
			const double alMax, 
			double &totalMass,
			int &jMin, 
			int &jMax )
{
	
	if( alMin > alMax ) return false;

	jMin=0;
	jMax=N-1;

	// find indices that bracket region ...
	bool minFound = false;
	for( int j = 0 ; j < N ; ++j ){
		if( !minFound && l[j] > alMin ){
			jMin = j - 1;
			minFound = true;
		} // end min if
		if( l[j] > alMax ){
			jMax = j;
			break;
		}
	} // end j for

	if( (jMax-jMin)%2 == 0 ) jMin++;		// for simpson's to work

	totalMass =	massIntegral(l,Fj,jMin,jMax);
	
	return true;
}// end massAnnulus


/*
 *		WRITE MASS
 *
 *		Appends time and total mass to file mass.out
 */
int writeMass(const double *l, const double *Fj,const double t , cnSolver &solver ){
	
	static bool good = true;
	double totalMass;
	int jMin, jMax;
	
	if(!good) return EXIT_FAILURE;
	
	ofstream fout("mass.out",ofstream::app);
	if(!fout){
		cerr << "WARNING IN MAIN ... Output file mass.out did not open properly." << endl;
		good = false;
		return EXIT_FAILURE;
	} // end err if

	massAnnulus(l,Fj,l_a-4,l_a+4,totalMass,jMin,jMax);

	fout << t << "\t" << massIntegral(l,Fj) << "\t" << totalMass << "\t"
	     << solver.Mdot(l,Fj,jMin) << "\t" << solver.Mdot(l,Fj,jMax) << endl;

	fout.close();
	
	return EXIT_SUCCESS;	
}// end writeMass

#endif
