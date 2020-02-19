#include <R.h>
 #include <math.h>
 void g_Boehm_JProteomeRes2014_310df0pf ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (100.0*x[3+i**k]+200.0*x[1+i**k]*p[0])/(x[3+i**k]+x[0+i**k]*p[0]+2.0*x[1+i**k]*p[0]) ;
y[1+i**l] = -(100.0*x[3+i**k]-200.0*x[4+i**k]*(p[0]-1.0))/((x[2+i**k]*(p[0]-1.0)-x[3+i**k])+2.0*x[4+i**k]*(p[0]-1.0)) ;
y[2+i**l] = (100.0*x[3+i**k]+100.0*x[0+i**k]*p[0]+200.0*x[1+i**k]*p[0])/(2.0*x[3+i**k]+x[0+i**k]*p[0]+2.0*x[1+i**k]*p[0]-x[2+i**k]*(p[0]-1.0)-2.0*x[4+i**k]*(p[0]-1.0)) ; 
}
}