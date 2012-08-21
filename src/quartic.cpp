#include <math.h>
#include <iostream>

using std::cout;
using std::endl;

double quartic( double A, double B, double C ){
	
	double t0 = pow( sqrt(3)*sqrt( 256.0*A*A*A*C*C*C + 27.0*A*A*B*B*B*B)
	                 + 9.0*A*B*B , 1.0/3.0 );
	double t1 = pow( 18.0 , 1.0/3.0 ) * A;
	double t2 = 4.0 * pow( 2.0/3.0 , 1.0/3.0 ) * C;
	double t3 = sqrt( t0/t1 - t2/t0 );
	
	return .5*t3*(1.0+sqrt(2.0*B/(A*pow(t3,1.5))-1.0));
}// end quartic

int main(){

	cout << quartic(1.0,1.0,1.0) << endl

	return 0;
} // end main
