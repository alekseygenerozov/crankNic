#include <stdio.h>
#include <stdlib.h>
#include "nr3.h"
#include "banded.h"

void printMat( MatDoub a ){
	int nr=a.nrows(),nc=a.ncols(),i,j;
	fprintf(stderr,"	[	\n");
	for(i=0;i<nr;i++){
		fprintf(stderr,"\n		");
		for(j=0;j<nr;j++)
			fprintf(stderr,"%.2g  ",a[i][j]);	
	}// end i for
	fprintf(stderr,"\n	]\n");
} // end printMat

void printVec( VecDoub a){
	int nr=a.size(),i;
	fprintf(stderr,"  [ ");
	for(i=0;i<nr;i++)
		fprintf(stderr,"%.2g  ",a[i]);
	fprintf(stderr,"]\n");
} // end printMat

int main(){

// We test our banded matrix solver on some small
// examples before using in our program

// TEST 1:  10x10 matrix with 1's on a 3-thick diagonal band
//          dotted into a vector of ones 

	int nr=10,nc=nr,i,j;
	MatDoub m(nr,nc);	// Identity matrix
	VecDoub x(nr,nc);	// solution
	VecDoub b(nr,nc);	// RHS
	for(i=0;i<nr;i++){
		b[i] = i + 1.0;
		for(j=0;j<nr;j++) 
			m[i][j] = (i==j)?1.0:0;
	} // end i for

	fprintf(stderr,"\n\n  b =");
	printVec(b);

	fprintf(stderr,"\n\n  m =");
	printMat(m);

	Bandec B(m,1,1);
	B.solve(b,x);

	fprintf(stderr,"\n\n x = m.b =");
	printVec(x);

	return 0;
}// end main
