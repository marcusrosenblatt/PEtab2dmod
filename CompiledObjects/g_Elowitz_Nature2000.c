#include <R.h>
 #include <math.h>
 void g_Elowitz_Nature2000_flofp2gb ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0]+x[7+i**k]*p[1] ; 
}
}