#ifndef INC_POLYSOLVE
#define INC_POLYSOLVE

/*
 *  QUARTIC 
 *
 *    Finds solutions to the equation:
 *
 *      Ax^4 - Bx - C = 0.0
 *
 *        !!!!! Assumes A,B > 0.0 !!!!!
 *
 *    Used by the Temperature solver
 */
double quartic( double A, double B, double C ) const
{
  static const double c0=sqrt(3.0),c1=pow(18.0,1.0/3.0),c2=4.0*pow(2.0/3.0,1.0/3.0);
  static double t0,t1,t2,t3,A2,B2,tmp;

  A2 = A*A;
  B2 = B*B;

  t0 = pow(c0*sqrt( A2*(256.0*A*C*C*C + 27.0*B2*B2))
                   + 9.0*a*B2 , 1.0/3.0 );
  t1 = c1*A;
  t2 = c2*C;
  t3 = sqrt( t0/t1 - t2/t0 );

  return .5*t3*(1.0+sqrt(2.0*B/(A*t3*t3*t3)-1.0));
}// end quartic



/*
 *	QUADRATIC
 *
 *		Finds positive solution to equation
 *
 *			Ax^2 + Bx - C = 0.0
 *
 *			!!! ASSUMES A > 0 and B,C < 0 !!!
 *
 */
void quadratic( double a, double b, double c ){
  return (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
}

#endif
