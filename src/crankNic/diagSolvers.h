#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
 *  Tridiagonal matrix solver (from wikipedia)
 */
void solveMatrix (int n, double *a, double *b, double *c, double *v, double *x)
{

	/*
	 * n - number of equations
	 * a - sub-diagonal -- indexed from 1..n-1
	 * b - the main diagonal
	 * c - sup-diagonal -- indexed from 0..n-2
	 * v - right part
	 * x - the answer
	 */

	for (int i = 1; i < n; i++){
		double m = a[i]/b[i-1];
		b[i] = b[i] - m*c[i-1];
		v[i] = v[i] - m*v[i-1];
	}

	x[n-1] = v[n-1]/b[n-1];

	for (int i = n - 2; i >= 0; i--)
		x[i]=(v[i]-c[i]*x[i+1])/b[i];
	
} // end solve matrix
