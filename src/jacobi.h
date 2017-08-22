#define MAXDIM    4    /* maximum size of the matrix */
#define MAXSWEEP  100   /* maximum number of iterations */
#define CONV       1.0E-10
#define CONV_INV   1.0E10
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
extern int jacobi(double a[MAXDIM][MAXDIM], int n, double d[],
		  double v[MAXDIM][MAXDIM], int nrot);
