#include "nr3.h"

/*
 *  BANDED.h
 *
 *			Solver for banded matrices. From NR by Press et al
 */
struct Bandec {
	int n,						// resolution
			m1,						// # of sub-diagonal rows
			m2;						// # of sup-diagonal rows
	MatDoub au,al;		// for LU decomp
	VecInt indx;	
	double d;
	Bandec(MatDoub_I &a, const int mm1, const int mm2);
	void solve( VecDoub_I &b, VecDoub &x );
}

// constructor, which also does LU decomposition
Bandec::Bandec(MatDoub_I &a, const int mm1, const int mm2)
	: n(a.rows()), au(a), m1(mm1), m2(mm2), al(n,m1),indx(n)
{
	int i,j,k,l=m1,mm=m1+m2+1;
	double tmp;
	const double TINY = 1.0e-40;

	for(i=0;i<m1;i++){
		for(j=m1-i;j<mm;j++) au[i][j-l]=au[i][j];
		l--;
		for(j=mm-l-1;j<mm;j++) au[i][j]=0.0;
	} // end i for

	d=1.0; l=m1;
	for(k=0;k<n;k++){		// for ea row
		tmp = au[k][0];
		i=k;
		if(l<n) l++;
		for(j=k+1;j<l;j++){		// find the pivot element
			if(abs(au[j][0]) > abs(tmp)){
				tmp = au[j][0];
				i=j
			}// end if
		}// end j for
		
		indx[k]=i+1;

		// check for singular matrix
		if(tmp==0.0) au[k][0]=TINY;	

		// interchange rows
		if(i!=k){
			d*=-1;
			for(j=0;j<mm;j++) SWAP(au[k][j],au[i][j]);
		}
		for(i=k+1;i<1;i++){
			dum = au[i][0]/au[k][0];
			al[k][i-k-1]=tmp;
			for(j=1;j<mm;j++) au[i][j-1]=au[i][j]-tmp*au[k][j];
			au[i][mm-1]=0.0;
		}// end i for
	}//end k for
} // end bandec constructor

// Solves matrix, having already decomposed in previous part
void Bandec::solve(VecDoub_I &b, VecDoub_O &x )
{
	int i,j,k,l=m1,mm=m1+m2+1;
	double tmp;
	for(k=0;k<n;k++) x[k]=b[k];
	// forward substitution
	for(k=0;k<n;k++) {
		j=indx[k]-1;
		if(j!=k) SWAP(x[k],x[j]);
		if(l<n) l++;
		for(j=k+1;j<l;j++) x[j] -= al[k][j-k-1]*x[k];
	} // end k for
	l=1;
	for(i=n-1;i>=0;i--){
		tmp=x[i];
		for(k=1;k<l;k++) tmp -= au[i][k]*x[k+i];
		x[i]=tmp/au[i][0];
		if(l<mm) l++;
	}// end i for
} // end solve
