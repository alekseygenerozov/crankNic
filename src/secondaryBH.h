#include "problemDomain.h"
#include "quadCoeffs.h"
#include "cubicSpline.h"

#ifndef INC_SECONDARY
#define INC_SECONDARY

const int STATIC  = 0;
const int DYNAMIC = 1;

class secondaryBH {
public:
	secondaryBH();
	double torque(const gasDisk&,const double,const double) const;
	void moveSecondary(const gasDisk&,const double,const double);
	double gaussTorqueInt(cubicSpline&,const gasDisk&,const double,const double,const double);
	int writeOut(const problemDomain&) const;
	double q;      // BH mass ratio
	double f;      // numerical parameter (see Armitage & Natarajan)
	double l_a;    // binary separation
	static const double sigma_conversion_factor = 1.61939E-6;	// S ~ 10^5 g/cm^2 @ behind secondary
	int position;  // static or dynamic (default static)
	int GW_loss;   // include GW wave hardening (default OFF)
}; // end secondaryBH



/*
 *	DEFAULT CONSTRUCTOR
 */
secondaryBH::secondaryBH()
	: q(0.1), f(0.01), l_a(10.0), position(STATIC) , GW_loss(OFF) {}


/*
 *  TORQUE
 *
 *		Smoothed torque density of secondary on disk, due to
 *		to linblad resonances.
 */
double secondaryBH::torque( const gasDisk &disk, const double l , const double M) const
{
#ifdef _QUAD_CHECK_
	if( l < 7.0 || l > 16.0 ) return 0.0;
	return 10.0;
#endif

	if( q == 0.0 ){
		return 0.0;
	}// end simple cas if
	double	tmp1 = f*q*q*M*M*0.5,
					l2 = l*l,
					la2 = l_a*l_a,
					lh2 = disk.h(l,M)*M,
					tmp2;

	if( l2 < la2 - lh2 ){
		tmp2 = l2/(l2-la2);
		return -tmp1*tmp2*tmp2*tmp2*tmp2/l2;
	}
	if( l2 < la2 ){
		tmp2 = l2/lh2;
		return -tmp1*tmp2*tmp2*tmp2/lh2;
	}
	if( l2 < la2 + lh2 ){
		tmp2 = la2/lh2;
		return tmp1*tmp2*tmp2*tmp2*tmp2/l2;
	}
	tmp2 = la2/(l2-la2);
	return 	tmp1*tmp2*tmp2*tmp2*tmp2/l2;
}// end tidal torque

/*
 *  GAUSS TORQUE INT
 *
 *    Performs Gauss Quadrature of FJ*Lambda/DJ 
 *    in region a < l < b. Needs an interpolation
 *    object for FJ, from which it computes DJ.
 *
 *    GQ coefficients are stored in quadCoeffs.h
 */
double secondaryBH::gaussTorqueInt( cubicSpline &FJ_cSpline,
                                    const gasDisk& disk,
                                    const double a,
                                    const double b ,
                                    const double M )
{
	if( b <= a ) return 0.0;

/*
// ---------------- FIXME
	int N = 10000;
	double ds = (b-a)/(N-1.0),s,f;
  double result = 0.0;

	f = FJ_cSpline.interp(a);
  result += f*torque(disk,a,M)/disk.Dj(f,a);
  for(int i= 1 ; i<N-1 ; i++ ){
		s = a + i*ds;
		f = FJ_cSpline.interp(s);
    result += (((i)%2==0)?2.0:4.0)*f*torque(disk,s,M)/disk.Dj(f,s);
	}
	f = FJ_cSpline.interp(b);
  result += f*torque(disk,b,M)/disk.Dj(f,b);

  return result*ds/3.0;
// ---------------- FIXME
*/


	static const size_t N_GAUSS = 48;
	double bMao2 = 0.5*(b-a), bPao2 = 0.5*(b+a),
		sum = 0.0,x,Ftmp;
	
	// first half of integral (splitting like this is optimal 
	// for spline fcn)
	for( size_t i = 0 ; i != N_GAUSS ; ++i ){
		x = bPao2 + bMao2*gaussQuad_x_96[i];
		Ftmp = FJ_cSpline.interp(x);
		sum += torque(disk,x,M)*Ftmp/disk.Dj(Ftmp,x)*gaussQuad_w_96[i];
	} // end i for

	// second half of integral
	for( size_t i = 0 ; i != N_GAUSS ; ++i ){
		x = bPao2 - bMao2*gaussQuad_x_96[i];
		Ftmp = FJ_cSpline.interp(x);
		sum += torque(disk,x,M)*Ftmp/disk.Dj(Ftmp,x)*gaussQuad_w_96[i];
	} // end i for

	return sum*bMao2;

}// end gaussTorqueInt



/*
 *  MOVE SECONDARY
 *
 *    Performs "back-reaction" of disk-secondary torque
 *    on the secondary blackhole, hardening the binary.
 *
 *    Integrates product of FJ and Torque density
 *    profile, then uses this to find dl_a/dt and updates
 *    l_a for next time step.
 */
void secondaryBH::moveSecondary( const gasDisk &disk, const double dt , const double M)
{
	if( position == STATIC || q == 0.0 ) return;  // do nothing

	const double LL = sqrt(disk.h(l_a,M));

	// interpolate data
	cubicSpline FJ_cSpline(disk.l,disk.Fj);

	// perform integrations in four regions, using 
	// 96-pt gaussian quadrature scheme
	double dldt = 0.0;
	dldt += gaussTorqueInt(FJ_cSpline,disk,disk.lMin,l_a-LL,M); // ( ISCO < l < l_a - LL )
	dldt += gaussTorqueInt(FJ_cSpline,disk,l_a+LL,disk.lMax,M); // ( l_a + LL < l < l_out )
	dldt += gaussTorqueInt(FJ_cSpline,disk,l_a-LL,l_a,M);       // ( l_a - LL < l < l_a )
	dldt += gaussTorqueInt(FJ_cSpline,disk,l_a,l_a+LL,M);       // ( l_a < l < l_a + LL )

	dldt *= (-1.0/(M*q)) * sigma_conversion_factor;

	if( GW_loss == ON )
		dldt -= 32.0*pow(M/l_a,7.0)*q/((1.0+q)*(1.0+q)*5.0);

	// update binary position
	l_a += dldt*dt;

} // end moveSecondary



/*
 *	WRITE OUT
 */
int secondaryBH::writeOut(const problemDomain &domain) const
{
	static bool first = true, good = true;
	if(!good) return EXIT_FAILURE;
	
	ofstream fout("secondary.out",first?ofstream::trunc:ofstream::app);
	if(!fout){
		good = false;
		cerr << "WARNING IN SECONDARY::writeOut ... Output file mass.out did not open properly." << endl;
		return EXIT_FAILURE;
	} // end err if
	
	if(first){
		fout << "#t	l_a	q" << endl;
		first = false;
	}

	fout << domain.t << "\t" << l_a << "\t" << q << endl;

  fout.close();

	return EXIT_SUCCESS;
}// end writeOut

#endif
