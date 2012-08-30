/*
 *	UNIT FACTORS
 *
 *		Designed to allow conversion from code units to normal units.
 *
 *		In addition the following dimensionless constants are used
 *		in the program to perform calculations of temperature/pressure/etc ...
 *
 *				GAMMA == Pressure Const / ( c * Omega const * sigma const )
 *				ETA   == fcn of fundamental consts, sigma and M
 *        KS    == kappa (opacity) * sigma constant
 *
 *		Finally, this program has a routine to print to units.out the unit factors
 */

#ifndef INC_UNITS
#define INC_UNITS

#include <fstream>

using std::endl;

class unitFactors {
public:
	unitFactors();	
	void writeOut() const;

	double gamma,M7,s5,mu,kappa,GM,c,mp,kb,sb,
		r,omega,H,nu,t,Lambda,sigma,P,T,l,FJ,DJ,eta,ks;

}; // end unitFactors

/*
 *	DEFAULT CONSTRUCTOR
 */
unitFactors::unitFactors()
{	
	gamma = 1.0E8;

	// PROBLEM PARAMETERS
	M7    = 1.0;               // mass of primary in 1E7 Msun
	s5    = 1.0;               // typical surface density in 1E5 g/cm^2
	mu    = 0.5;               // mean molecular weight
	kappa = .35;               // opacity, cm^2 / g

	// UNIVERSAL CONSTANTS
	GM    = 1.32732643E33*M7;  // cm^3/s^2
	c     = 2.99792458E10;     // cm/s
	mp    = 1.67262158E-24;    // grams
	kb    = 1.3806503E-16;     // ergs/K
	sb    = 5.670373E-5;       // stefan-boltzman, cgs 

	// NORMALIZATION FACTORS
	r     = GM/(c*c);                 // cm
	omega = c*c*c/GM;                 // s
	H     = r;                        // cm
	
	nu     = GM/c;                    // cm^2/s
	t      = GM/(c*c*c);              // s
	Lambda = c*c;                     // cm^2/s^2

	sigma = 1.0E5*s5;                 // grams / cm^2
	P     = c*c*c*c*sigma/(gamma*GM); // dynes/cm^2 
	T     = 2.0*mu*mp*c*c/(gamma*kb); // Kelvin

	l  = GM/c;                        // cm^2/s
	FJ = GM*GM/(c*c)*sigma;           // ergs
	DJ = GM*c;                        // cm^4/s^3

	// DIMENSIONLESS CONSTANTS
	eta   = 64.0*sb*mu*mu*mu*mu*mp*mp*mp*mp*c*c*c
                        *GM/(3.0*gamma*gamma*gamma*kb*kb*kb*kb*sigma);
	ks    = kappa*sigma;    // <-- and ^ used in temp/press calcs
} // end constructor

/*
 *    WRITE MASS
 *
 *    Appends time and total mass to file mass.out
 */
void unitFactors::writeOut() const 
{


  std::ofstream fout("units.out");
	fout << "gamma = " << gamma << endl;
	fout << endl;

	fout << "// PROBLEM PARAMETERS" << endl;
	fout << "M7    = " << M7 << endl;
	fout << "s5    = " << s5 << endl;
	fout << "mu    = " << mu << endl;
	fout << "kappa = " << kappa << endl;
	fout << endl;

	fout << "// UNIVERSAL CONSTANTS" << endl;
	fout << "GM = " << GM << endl; 
	fout << "c  = " << c << endl;
	fout << "mp = " << mp << endl;
	fout << "kb = " << kb << endl;
	fout << "sb = " << sb << endl;
	fout << endl;

  fout << "// NORMALIZATION FACTORS" << endl;
	fout << "r     = " << r << endl;
	fout << "omega = " << omega << endl;
	fout << "H     = " << H << endl;
	fout << endl;
	fout << "nu     = " << nu << endl;
	fout << "t      = " << t << endl;
	fout << "Lambda = " << Lambda << endl;
	fout << endl;
	fout << "sigma = " << sigma << endl;
	fout << "P     = " << P << endl;
	fout << "T     = " << T << endl;
	fout << endl;
	fout << "l  = " << l << endl;
	fout << "FJ = " << FJ << endl;
	fout << "DJ = " << DJ << endl;
	fout << endl;

	fout << "// DIMENSIONLESS CONSTANTS" << endl;
	fout << "eta = " << eta << endl;
	fout << "ks  = " << ks << endl;
	fout << endl;

  fout.close();

}// end writeMass

#endif
