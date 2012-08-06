#ifndef INC_SECONDARY
#define INC_SECONDARY

const int STATIC  = 0;
const int DYNAMIC = 1;

class secondaryBH {
public:
	secondaryBH();
	double torque(const gasDisk &disk, const double l, const double M) const;
	
	double q;      // BH mass ratio
	double f;      // numerical parameter (see Armitage & Natarajan)
	double l_a;    // binary separation
	int position;  // static or dynamic (default static)
}; // end secondaryBH



/*
 *	DEFAULT CONSTRUCTOR
 */
secondaryBH::secondaryBH()
	: q(0.1), f(0.01), l_a(10.0), position(STATIC) {}


/*
 *  TORQUE
 *
 *		Smoothed torque density of secondary on disk, due to
 *		to linblad resonances.
 */
double secondaryBH::torque( const gasDisk &disk, const double l , const double M) const
{

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

#endif
