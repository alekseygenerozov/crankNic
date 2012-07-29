#ifndef INC_MASS
#define INC_MASS

/*
 *	MASS INTEGRAL
 *
 *		Integrates total disk mass, using Simpson's
 *		1/3rd rule on a logarithmic grid
 *
 */
double massIntegral(double *l, double *Fj){

  double result = 0.0,
    evenCoeff = (2.0*lambda-1.0)/(6.0*lambda*lambda*lambda) + (2.0-lambda)/6.0,
    oddCoeff  = (lambda+1.0)*(lambda+1.0)/(6.0*lambda*lambda);

  result += (2.0-lambda)/6.0*dmdl(Fj[0],l[0]);
  for(int j=1;j<N-1;j++)
      result += ((j%2==0)?evenCoeff:oddCoeff)*dmdl(Fj[j],l[i])*pow(lambda,j);
  result += (2.0*lambda-1.0)/6.0*pow(lambda,N-4)*dmdl(Fj[N-1],l[N-1]);

  return result*dr*(lambda+1.0);
}// end massIntegral

#endif
