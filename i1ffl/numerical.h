#include "RngStream.h" 

/* Normally distributed random number with mean zero and variance 1 */
extern float gasdev(RngStream);

/* Returns an exponentially distributed, positive, random deviate of unit mean */
extern float expdev(RngStream);

/*Newton-Raphson root-finding method*/
extern float rtsafe(void (*funcd)(float, int, float, float*, float*, float*, float*, float*, int), 
		    int, float, float *, float*, float*, float, float, int);

/*find max or min*/
extern void find_max(float*, int, int, float*, int*, float);
extern void find_x(float*, int, int, float, float*, int);
