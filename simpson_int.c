//#include <math.h>
#define EPS_sim 1.0e-4
#define JMAX 20
#include "trapezoidal_sim.c"


double qsimp(double (*func)(double ,double []), double a, double b,double vararg[])
/*Returns the integral of the function func from a to b. The parameters EPS can be set to the
desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
number of steps. Integration is performed by Simpson¡¯s rule.*/
{
double trapzd_sim(double (*func)(double, double []), double a, double b, int n, double vararg[]);
void nrerror(char error_text[]);
int j;
double s,st,ost=0.0,os=0.0;
for (j=1;j<=JMAX;j++) {
st=trapzd_sim(func,a,b,j,vararg);
s=(4.0*st-ost)/3.0;   // Compare equation (4.2.4), above.
if (j > 5) // Avoid spurious early convergence.
if (fabs(s-os) < EPS_sim*fabs(os) ||
(s == 0.0 && os == 0.0)) return s;
os=s;
ost=st;
}
nrerror("Too many steps in routine qsimp");
return 0.0;  //Never get here.
}
