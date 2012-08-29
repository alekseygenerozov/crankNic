#include "problemDomain.h"
#include "gasDisk.h"
#include "secondaryBH.h"
#include "udSolver.h"

#ifndef INC_MASS
#define INC_MASS

/*
 *	MASS INTEGRAL
 *
 *		Integrates total disk mass, using Simpson's
 *		1/3rd rule on a logarithmic grid
 *
 */
double massIntegral( const gasDisk &disk,
                     int jMin = 0, 
                     int jMax = UNSET )
{
	if( jMax == UNSET ) jMax = disk.N - 1;

	double result = 0.0, lambda = disk.lambda,
		evenCoeff = (2.0*lambda-1.0)/(6.0*lambda*lambda*lambda) + (2.0-lambda)/6.0,
		oddCoeff  = (lambda+1.0)*(lambda+1.0)/(6.0*lambda*lambda);

	result += (2.0-lambda)/6.0*disk.dmdl(jMin);
	for(int j=jMin+1;j<jMax;j++)
		result += (((j-jMin)%2==0)?evenCoeff:oddCoeff)*disk.dmdl(j)*pow(lambda,j);
	result += (2.0*lambda-1.0)/6.0*pow(lambda,(int)disk.N-4)*disk.dmdl(jMax);

	return result*disk.dl*(lambda+1.0);
}// end massIntegral

/*
 *	MASS ANNULUS
 *
 *		Returns the mass contained within two radii 
 *
 */
bool massAnnulus(
			const gasDisk &disk,
			const double alMin, 
			const double alMax, 
			double &totalMass,
			int &jMin, 
			int &jMax )
{
	
	if( alMin > alMax ) return false;

	jMin=0;
	jMax=disk.N-1;

	// find indices that bracket region ...
	bool minFound = false;
	for( size_t j = 0 ; j < disk.N ; ++j ){
		if( !minFound && disk.l[j] > alMin ){
			jMin = j - 1;
			minFound = true;
		} // end min if
		if( disk.l[j] > alMax ){
			jMax = j;
			break;
		}
	} // end j for

	if( (jMax-jMin)%2 == 0 ) jMin++;		// for simpson's to work

	totalMass =	massIntegral(disk,jMin,jMax);
	
	return true;
}// end massAnnulus


/*
 *		WRITE MASS
 *
 *		Appends time and total mass to file mass.out
 */
int writeMass( const problemDomain &domain, 
               const gasDisk &disk,
               const secondaryBH &secondary,
               const udSolver &solver )
{
	
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

	massAnnulus(disk,secondary.l_a-4,secondary.l_a+4,totalMass,jMin,jMax);

	fout << domain.t << "\t" << massIntegral(disk) << "\t" << totalMass << "\t"
	     << solver.Mdot(domain,disk,secondary,jMin) << "\t" << solver.Mdot(domain,disk,secondary,jMax) << endl;

	fout.close();
	
	return EXIT_SUCCESS;	
}// end writeMass

#endif
