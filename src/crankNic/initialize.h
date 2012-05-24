#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"

int initialize( double *r, double *sigma ){

	/*
	 * PROBLEM 1
	 *
	 *		Delta function response for S&S Disk,
	 *		constant viscosity, zero torque
	 *
	 *				>> Can be compared to analytic
	 *					 expression
	 */
	if( problemType == 1 ){
		// We intialize from file ...
		FILE* fp = fopen("analytic_T0.dat","r");
		if(!fp){
			fprintf(stderr,"ERROR IN INITIALIZE.H --- Failed to Open IC file:\n");
			return EXIT_FAILURE;
		} // end error if
		double tmp1,tmp2;
		for( int i = 0 ; i < N ; i++ ){
			fscanf(fp,"%lg",&tmp1);
			fscanf(fp,"%lg",&tmp2);
			r[i] = tmp1;
			sigma[i] = tmp2;
		}// end i for	
		fclose(fp);
	} // end delta-function initialize

	/*
	 * PROBLEM 2
	 *
	 *		Static Secondary, constant viscosity
	 *			Following Liu & Schapiro
	 */
	else if( problemType == 2 ){

		for( int i = 0 ; i < N ; i++ ){
			if( lambda == 1.0 ){
				r[i] = rMin + i*dr;
			} else {
		    if( i == 0 )
					r[i] = rMin;
		    else
					r[i] = rMin + dr*( pow(lambda,i) - 1.0 )/(lambda-1.0);
			}// end lambda if/else
			
			sigma[i] = (r[i]-rMin)/(rMax-rMin);
		}// end i for

/*	EVENTUAL IMPLEMENTATION  -------------------------------------

		double g,g_star,s,ds,s1=sqrt(a/rMax),s2=sqrt(rMin/rMax),ff,dp;
		g = 2.0/3.0*f*q*q*sqrt(M)*sqrt(rMax)/nu(rMax);
		double F[N],G_H[N];
		F[N-1] = 0.0;
		G_H[N-1] = 0.0;

		bool firstH = true;
		for( int i = N-2 ; i >= 0 ; i-- ){
			s = sqrt(r[i]/rMax);
			ds = s - sqrt(r[i+1]/rMax);
			dp = max(fabs(s*s-s1*s1),s*s*dhdr);
			g_star = (s>s1?g*pow(s1,8):-g*pow(s,8));
			ff = g_star/nu(r[i])*pow(1.0/dp,4);

			F[i] = F[i+1] - ff*ds;
			if( s > s1)
				G_H[i] = G_H[i+1] - exp(F[i+1])*ds;
			else
				if( firstH ){
					G_H[i] = log(G_H[i+1]) - ds*exp(F[i+1]-log(G_H[i+1]));
					firstH = false;
				} else
					G_H[i] = G_H[i+1] - ds*exp(F[i+1]-G_H[i+1]);

			fprintf(stderr,"i,r,F,G = %d\t%g\t%g\t%g\n",i,r[i],F[i],G_H[i]);

		} // end i for	

		double G0 = exp(G_H[0]);
		for( int i = 0 ; i < N ; i++ ){
			s = sqrt(r[i]/rMax);

			if( s < s1 )
				sigma[i] = exp(-F[i])/nu(r[i])*(1.0-exp(G_H[i]-G_H[0]))/s;
			else
				sigma[i] = exp(-F[i])/nu(r[i])*(1.0-G_H[i]/G0)/s;
		}// end i for
----------------------------------------------------------------------*/

		sigma[0] = 1E-8;	
	} // end linTorq test problem


	/*
	 *	PROBLEM 3 -- Square Pulse
	 *		To test the equation
	 *			u_t = c/r * u_x	
	 *					(i.e. the advective term)
	 */
	else if( problemType == 3 ){
		for( int i = 0 ; i < N ; i++ ){
			if( lambda == 1.0 ){
				r[i] = rMin + i*dr;
			} else {
				if( i == 0 )
					r[i] = rMin;
				else
					r[i] = rMin + dr*( pow(lambda,i) - 1.0 )/(lambda-1.0);
			}// end lambda if/else

			if( r[i] > (rMax-rMin)*0.7 && r[i] < (rMax-rMin)*(0.8) )
				sigma[i] = 1.0;
			else
				sigma[i] = 0.1;
		}// end i for
	}// end square test problem

	// if not explicitly set, normalize torque at outer boundary
	if(-1.0==nu0)
		nu0 = (n_v==0?1.0:pow(1.0/r[N-3],n_v));

	// If Dirichlet BVs and no value set, use IC file:
	if(-1.0==outer_bndry_value)
		outer_bndry_value = sigma[N-1];
	if(-1.0==inner_bndry_value)
		inner_bndry_value = sigma[0];

	return EXIT_SUCCESS;

} // end main
