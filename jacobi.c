/*************************************************************
*                                                            *
* jacobi -- compute eigenvalues and eigenvectors of a real   *
*           symmetric matrix by the Jacobi method.           *
*                                                            *
* Arguments                                                  *
*     a[0...n][0...n]  a real symmetric matrix               *
*     n                size of a                             *
*     d[0...n]         eigenvalues of a                      *
*     v[0...n][0...n]  a matrix whose columns contain the    *
*                      normalized eigenvectors of a          *
*     nrot             number of Jacobi rotations            *
*                                                            *
*                                                            *
*      Press, W.H. et al. 1992, Numerical Recipes in C,      *
*      2nd ed. (Cambridge University Press).                 *
*written out by Charlotte Deane                              *
23rd Jan 1998                                                *
*                                                            *
*************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jacobi.h"

int jacobi(double a[MAXDIM][MAXDIM], int n, double D[],
	   double V[MAXDIM][MAXDIM], int nrot)
{
  void eigsort(double D[], double V[MAXDIM][MAXDIM], int n);

  int j,iq,ip,i,k;
  double tresh;      /* In the first three sweeps, the pq rotaion 
			is carried out only if |a[ip][jq]| > tresh */
  double g;
  double theta,tau,t,sm,s,h,c;

  if (n > MAXDIM) {
    fprintf(stderr, "Error in routine jacobi\n");
    fprintf(stderr, "The size of the matrix must be smaller than MAXDIM.\n");
    return -1;
  }

  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) {
      V[ip][iq]=0.0;
    }
    V[ip][ip]=1.0;
  }

  for (ip=0; ip<n; ip++) {
    D[ip] = a[ip][ip];
  }
  nrot=0;

  for (i=0; i<MAXSWEEP; i++) {

/*
for(k=0;k<n;k++){
for(j=0;j<n;j++){
printf ("%f a k%d j%d \n",a[k][j],k,j);
}
}
*/


    sm=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++)
	sm += fabs(a[ip][iq]);
    }

    if (sm < CONV) {
      eigsort(D, V, n);
      return 0;
    }

    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;

    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++) {
	g = 100.0 * CONV_INV * fabs(a[ip][iq]);
	if (i > 4 && fabs(D[ip]) > g && fabs(D[iq]) > g)
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h = D[iq]-D[ip];
	  if (fabs(h) > g) 
	    t = (a[ip][iq])/h;
	  else {
	    theta = 0.5 * h/(a[ip][iq]);
	    t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt(1 + t*t);
	  s = t*c;
	  tau = s/(1.0+c);
	  h = t * a[ip][iq];
	  D[ip] -= h;
	  D[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=0; j<=ip-1; j++) {
	    ROTATE(a,j,ip,j,iq)
	  }
	  for (j=ip+1; j<=iq-1; j++) {
	    ROTATE(a,ip,j,j,iq)
	  }
	  for (j=iq+1; j<n; j++) {
	    ROTATE(a,ip,j,iq,j)
	  }
	  for (j=0; j<n; j++) {
	    ROTATE(V,j,ip,j,iq)
	  }
	  ++(nrot);
	}
      }
    }
  }
  fprintf(stderr, "Too many iterations in routine jacobi\n");
  return -1;
}
#undef ROTATE

void eigsort(double D[], double V[MAXDIM][MAXDIM], int n) 
{
  int i, j, k;
  double p;

  for (i=0; i<n-1; i++) {
    k = i;
    p = D[k];
    for (j=i+1; j<n; j++) {
      if (D[j] >= p) {
	k = j;
	p = D[k];
      }
    }
    if (k != i) {
      D[k] = D[i];
      D[i] = p;
      for (j=0; j<n; j++) {
	p = V[j][i];
	V[j][i] = V[j][k];
	V[j][k] = p;
      }
    }
  }
}
