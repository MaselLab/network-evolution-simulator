/* random number generator from Numerical Recipes */
extern float ran1(long *seed);
/* Normally distributed random number with mean zero and variance 1 */
extern float gasdev(long *seed);
// poisson distributed
extern float poidev(float xm, long *seed);
/* Returns an exponentially distributed, positive, random deviate of unit mean */
extern float expdev(long *seed);
// binomially distributed
extern float bnldev(float pp, int n, long *seed);


