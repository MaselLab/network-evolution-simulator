#ifndef NUMERICAL
#define NUMERICAL

#include "netsim.h"

/* Normally distributed random number with mean zero and variance 1 */
float gasdev(RngStream);

/* Returns an exponentially distributed, positive, random deviate of unit mean */
float expdev(RngStream);

/*Newton-Raphson root-finding method*/
float rtsafe(void (*funcd)(float, int, float, float*, float*, float*, float*, float*), 
		    int, float, float *, float*, float*, float, float);

#endif

