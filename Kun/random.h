#include "RngStream.h"
/* random number generator from Numerical Recipes */
//extern float ran1(long *seed);
/* Normally distributed random number with mean zero and variance 1 */
extern float gasdev(RngStream);
// poisson distributed
extern float poidev(float xm, RngStream);
/* Returns an exponentially distributed, positive, random deviate of unit mean */
extern float expdev(RngStream);
// binomially distributed
extern float bnldev(float pp, int n, RngStream);
/* quantile function of normal distribution*/
extern float qz(float, float, float, float);

