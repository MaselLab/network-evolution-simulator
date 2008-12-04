#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"

long seed = 28121; 

int main(int argc, char *argv[])
{
  int i;

  for (i=0; i<10; i++) {
    printf("%d ran1=%g\n", i, ran1(&seed));
  }
}

/* function prototypes for reference:

   float ran1(long *seed)
   float gasdev(long *seed)
   float gammln(float xx)
   float poidev(float xm, long *seed)
   float expdev(long *seed)
   float bnldev(float pp, int n, long *seed) */



