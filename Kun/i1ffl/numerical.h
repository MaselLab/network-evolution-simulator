#include "RngStream.h" 
#include "netsim.h"

/* uniformally distributed random number generator from Numerical Recipes. 
 * Replaced by the uniform random number generator from RngStream */
//extern float ran1(long *seed);

/* Normally distributed random number with mean zero and variance 1 */
extern float gasdev(RngStream);
// poisson distributed
//extern float poidev(float xm, RngStream);

/* Returns an exponentially distributed, positive, random deviate of unit mean */
extern float expdev(RngStream);

// binomially distributed
//extern float bnldev(float pp, int n, RngStream);

/* cumulative probability function of norma distribution*/
extern float pnorm(float);
//extern float qz(float, float, float, float);

/* quantile function of normal distribution*/
extern float qnorm7(float, float);
//extern float qnorm4(float, float);

/*Newton-Raphson root-finding method*/
extern float rtsafe(void (*funcd)(float, int, float, float*, float*, float*, float*, float*, int), 
		    int, float, float *, float*, float*, float x1, float x2, float xacc, int);

/*find max or min*/
extern void find_max(float*, int, int, float*, int*);
extern void find_x(float*, int, int, float, float*, int);
/*Wilcoxon two-sample test*/
extern float Wilcoxon_test(Genotype*, Genotype*, int, int, float *, float *, float *);

/*support functions of Wilcoxon test*/
void quick_sort(float *, int *, int);
void quick_sort2(float *, int);
