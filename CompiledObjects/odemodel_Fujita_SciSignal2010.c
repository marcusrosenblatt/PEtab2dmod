/** Code auto-generated by cOde 1.0.0 **/
#include <R.h> 
 #include <math.h> 

static double parms[26];
static double forc[0];
static double cons[0];
static double eventcounter[2];
static double range[2];

#define nGridpoints 2 
#define nSplines 0 
#define precision 1e-05 

#define EGF_step parms[0] 
 #define EGF_ramp parms[1] 
 #define reaction_1_k1 parms[2] 
 #define reaction_1_k2 parms[3] 
 #define EGFR_turnover parms[4] 
 #define reaction_9_k1 parms[5] 
 #define reaction_2_k1 parms[6] 
 #define reaction_2_k2 parms[7] 
 #define reaction_3_k1 parms[8] 
 #define reaction_4_k1 parms[9] 
 #define reaction_7_k1 parms[10] 
 #define reaction_5_k1 parms[11] 
 #define reaction_5_k2 parms[12] 
 #define reaction_6_k1 parms[13] 
 #define reaction_8_k1 parms[14] 
 #define EGF_impulse parms[15] 
 #define y0_0 parms[16] 
 #define y1_0 parms[17] 
 #define y2_0 parms[18] 
 #define y3_0 parms[19] 
 #define y4_0 parms[20] 
 #define y5_0 parms[21] 
 #define y6_0 parms[22] 
 #define y7_0 parms[23] 
 #define y8_0 parms[24] 
 #define y9_0 parms[25] 
#define tmin range[0]
#define tmax range[1]


void odemodel_Fujita_SciSignal2010_initmod(void (* odeparms)(int *, double *)) {
	 int N=26;
	 odeparms(&N, parms);
	 for(int i=0; i<2; ++i) eventcounter[i] = 0;
}

void odemodel_Fujita_SciSignal2010_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void odemodel_Fujita_SciSignal2010_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = -1.0*(((EGF_step+y[9]+EGF_ramp*time/3600.0)*y[0]*reaction_1_k1-y[1]*reaction_1_k2))-1.0*(y[0]*EGFR_turnover)+1.0*(68190.0*EGFR_turnover);
 	 ydot[1] = 1.0*(((EGF_step+y[9]+EGF_ramp*time/3600.0)*y[0]*reaction_1_k1-y[1]*reaction_1_k2))-1.0*(y[1]*reaction_9_k1);
 	 ydot[2] = -1.0*((y[3]*y[2]*reaction_2_k1-y[4]*reaction_2_k2))+1.0*(y[4]*reaction_3_k1)-1.0*(y[2]*reaction_4_k1)+1.0*(y[1]*reaction_9_k1);
 	 ydot[3] = -1.0*((y[3]*y[2]*reaction_2_k1-y[4]*reaction_2_k2))+1.0*(y[5]*reaction_7_k1);
 	 ydot[4] = 1.0*((y[3]*y[2]*reaction_2_k1-y[4]*reaction_2_k2))-1.0*(y[4]*reaction_3_k1);
 	 ydot[5] = 1.0*(y[4]*reaction_3_k1)-1.0*((y[6]*y[5]*reaction_5_k1-y[7]*reaction_5_k2))+1.0*(y[7]*reaction_6_k1)-1.0*(y[5]*reaction_7_k1);
 	 ydot[6] = -1.0*((y[6]*y[5]*reaction_5_k1-y[7]*reaction_5_k2))+1.0*(y[8]*reaction_8_k1);
 	 ydot[7] = 1.0*((y[6]*y[5]*reaction_5_k1-y[7]*reaction_5_k2))-1.0*(y[7]*reaction_6_k1);
 	 ydot[8] = 1.0*(y[7]*reaction_6_k1)-1.0*(y[8]*reaction_8_k1);
 	 ydot[9] = 1.0*(0.0);

}

/** Event function **/
void odemodel_Fujita_SciSignal2010_myevent(int *n, double *t, double *y) {

	 double time = *t;

	 if(*t == 0 & eventcounter[0] == 0) {
		y[9] = EGF_impulse;
		eventcounter[0] = eventcounter[0] + 1.;
	 }

	 if(*t == 60.0 & eventcounter[1] == 0) {
		y[9] = 0;
		eventcounter[1] = eventcounter[1] + 1.;
	 }


}
