#include <R.h>
 #include <math.h>
 void g_Elowitz_Nature2000_deriv_t3wiaz6k ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (p[1])*(x[16+i**k]) ;
y[1+i**l] = (p[1])*(x[25+i**k]) ;
y[2+i**l] = (p[1])*(x[34+i**k]) ;
y[3+i**l] = (p[1])*(x[43+i**k]) ;
y[4+i**l] = (p[1])*(x[52+i**k]) ;
y[5+i**l] = (p[1])*(x[61+i**k]) ;
y[6+i**l] = (p[1])*(x[70+i**k]) ;
y[7+i**l] = (p[1])*(x[79+i**k]) ;
y[8+i**l] = (p[1])*(x[88+i**k]) ;
y[9+i**l] = (p[1])*(x[97+i**k])+1.0 ;
y[10+i**l] = (p[1])*(x[106+i**k])+x[7+i**k] ;
y[11+i**l] = (p[1])*(x[115+i**k]) ;
y[12+i**l] = (p[1])*(x[124+i**k]) ;
y[13+i**l] = (p[1])*(x[133+i**k]) ;
y[14+i**l] = (p[1])*(x[142+i**k]) ;
y[15+i**l] = (p[1])*(x[151+i**k]) ;
y[16+i**l] = (p[1])*(x[160+i**k]) ;
y[17+i**l] = (p[1])*(x[169+i**k]) ;
y[18+i**l] = (p[1])*(x[178+i**k]) ;
y[19+i**l] = (p[1])*(x[187+i**k]) ;
y[20+i**l] = (p[1])*(x[196+i**k]) ; 
}
}