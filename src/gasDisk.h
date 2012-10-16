#include "problemDomain.h"

#ifndef INC_DISK
#define INC_DISK

// BOUNDARY CONDITIONS
const int NEUMANN = 0;
const int DIRICHLET = 1;
const int ZERO = 0;
const int SELF_SIM = 1;
const int UNSET = -1;

// DISK PHYSICS OPTIONS
const int PWR_LAW   = 0;
const int BETA_DISK = 1;

class gasDisk {
public:
	gasDisk();
	void resize(size_t N, double lMin, double lMax, double l_a);
	void manualResize(size_t N, double lMin, double lMax, double lambda);
	void buildGrid();
	double dmdl(const size_t j) const {return Fj[j]/DJ[j];};
	
	// RESOLUTION PARAMETERS
	size_t N;         // Size of Simulation
	double lambda; // stretch factor of log grid
	double lMax;   // outer boundary
	double lMin;   // inner boundary
	double dl;     // cell size
	double dl2;    // cell size squared
	int STENCIL;   // For Gradient derivative term

	// PHYSICAL PARAMETERS
	double D0;          // Diffsn cnst (dflt: t=tau for const visc case)
	double nd;          // diffsn Fj power-law index
	double np;          // diffsn l pwr-law indx (dflt = Om_k, const visc)
	double dhdr;        // r/h for disk scale height
	double alpha;       // viscosity param
	int visc_model;	    // pwr law (default) or beta disk

	vDoub Fj,l,T,H,DJ;

	int outer_bndry_type, inner_bndry_type, outer_bndry_laplacian, inner_bndry_laplacian;
	double outer_bndry_value, inner_bndry_value;
	
	double density_floor;	// force density positive
}; // end disk



/*
 *	DEFAULT CONSTRUCTOR
 */
gasDisk::gasDisk()
	: N(500), lambda(1.0), lMax(2.0), lMin(0.1), STENCIL(0), D0(1.0/16.0),
	  nd(0.0), np(-2.0), dhdr(0.1), visc_model(PWR_LAW), alpha(0.1),
		density_floor(UNSET),
	  outer_bndry_type(NEUMANN),   inner_bndry_type(NEUMANN),
	  outer_bndry_laplacian(ZERO), inner_bndry_laplacian(ZERO),
	  outer_bndry_value(UNSET),    inner_bndry_value(UNSET)
{
	buildGrid();
}// end default constructor



/*
 *	RESIZE
 *
 *		Given # of grid cells, domain bounds and position of secondary,
 *		builds the radial grid with resolution optimized at the secondary
 */
void gasDisk::resize(size_t N, double lMin, double lMax, double l_a)
{
	if( N < 0 ) throw("Gridsize must be positive");
	this->N = N;
	if( lMin < 0 || lMin > lMax ) throw("lMin must be < lMax, positive");
	this->lMin = lMin;
	this->lMax = lMax;
	
	// optimize lambda ... 
	int n = (int)1E5;
	double lamMin = 1.0 - 1E-5, lamMax = 1.05, optLam = lamMin, res, lam, 
		dLam = (lamMax - lamMin)/(n-1.0), optRes = 1E99;
	for( int i = 0 ; i != n ; ++i ){
		lam = i*dLam + lamMin;
		res = (lam-1.0)*(l_a-lMin+(lMax-lMin)/(pow(lam,N)-1.0));
		if( res < optRes ){
			res = optRes;
			optLam = lam;
		} // end optRes if
	}// end i for
	this->lambda = lam;
	
	buildGrid();
}// end resize



/*
 *	MANUAL RESIZE
 *
 *		Resive grid, setting the logarithmic stretch factor by hand
 */
void gasDisk::manualResize(size_t N, double lMin, double lMax, double lambda)
{
	if( N < 0 ) throw("Gridsize must be positive");
	this->N = N;
	if( lMin < 0 || lMin > lMax ) throw("lMin must be < lMax, positive");
	this->lMin = lMin;
	this->lMax = lMax;
	if( lambda < 1.0 || lambda > 1.5 ) throw("lambda must be in [ 1.0 , 1.5 ]");
	this->lambda = lambda;

	buildGrid();
}// end manual_resize



/*
 *	BUILD GRID
 *
 *		Sets up radial grid points and resets FJ to all zeros
 */
void gasDisk::buildGrid()
{
	l.resize(N);
	Fj.resize(N);   // clears current values
	T.resize(N);		// ""
	H.resize(N);    // ""
	DJ.resize(N);		// "" 
	
	// calcualte innermost grid cell size
	if( lambda == 1.0 ){
		dl = (lMax-lMin)/(N-1.0);
	} else {
		dl = (lMax-lMin)*(lambda-1.0)/(pow(lambda,N-1)-1.0);
	}
	dl2  = dl*dl;
	cerr << "dl = " << dl << endl;

	// setup grid
	for( size_t j = 0 ; j < N ; j++ ){
			if( lambda == 1.0 ){
				l[j] = lMin + j*dl;
			} else {
				if( j == 0 )
					l[j] = lMin;
				else
					l[j] = lMin + dl*(pow(lambda,j)-1.0)/(lambda-1.0);
			}// end lambda if/else
	}// end j for
}// end buildGrid

#endif
