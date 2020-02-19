#include <R.h>
 #include <math.h>
 void g_Zheng_PNAS2012_50ltj0b7 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = x[0+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[1+i**l] = x[1+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[2+i**l] = x[4+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[3+i**l] = x[8+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[4+i**l] = x[2+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[5+i**l] = x[3+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[6+i**l] = x[6+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[7+i**l] = x[14+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[8+i**l] = x[5+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[9+i**l] = x[7+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[10+i**l] = x[10+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[11+i**l] = x[11+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[12+i**l] = x[9+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[13+i**l] = x[13+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ;
y[14+i**l] = x[12+i**k]/(x[0+i**k]+x[1+i**k]+x[4+i**k]+x[8+i**k]+x[2+i**k]+x[3+i**k]+x[6+i**k]+x[14+i**k]+x[5+i**k]+x[7+i**k]+x[10+i**k]+x[11+i**k]+x[9+i**k]+x[13+i**k]+x[12+i**k]) ; 
}
}