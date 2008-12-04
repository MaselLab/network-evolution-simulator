#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"

long seed = 28121; 

int main(int argc, char *argv[])
{
  int i;
  int r;

  for (i=0; i<100; i++) {
    printf("%d ran1=%g\n", i, ran1(&seed));
    printf("%d gasdev=%g\n", i, gasdev(&seed));
    r = rand()%10;
    printf("%d rand(): %d\n", i, r);
  }
}

/* function prototypes for reference:

   float ran1(long *seed)
   float gasdev(long *seed)
   float gammln(float xx)
   float poidev(float xm, long *seed)
   float expdev(long *seed)
   float bnldev(float pp, int n, long *seed) */



