#include <R.h>
 #include <math.h>
 void g_Fujita_SciSignal2010_kdyawnlv ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0]*(x[5+i**k]+x[7+i**k]) ;
y[1+i**l] = p[1]*(x[2+i**k]+x[4+i**k]) ;
y[2+i**l] = x[8+i**k]*p[2] ; 
}
}