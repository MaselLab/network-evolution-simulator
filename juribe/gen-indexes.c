#include <stdio.h>
#include <math.h>

#define TFBS 4

/* this simply prints the binary representation */
void printBinaryRepresentation(int decimal) {
   int binary = 0;
   int place = 0;
   int num = decimal;

   while (num != 0)
   {
      binary = binary + (num%2 * pow(10, place));
      num = num /2;
      place = place + 1;
   }
   printf("%.*d (%*d)", TFBS, binary, TFBS-1, decimal);
}

int countBits(unsigned int v) {

  unsigned int c; // c accumulates the total bits set in v
  for (c = 0; v; c++) {
      v &= v - 1; // clear the least significant bit set
  }
  return (c);
}

/* this gets the bit at position pos */
int getBit(unsigned int val, int pos) {
  return (val >> pos) & 1;
}

/* this set the bit at position pos */
int setBit(unsigned int val, int pos) {
  unsigned int mask = 1 << pos;
  return mask | val;
}

int main(int argc, char *argv[])
{
  unsigned int N = pow(2, TFBS);
  unsigned int i = 0;

  /* loop through all 2^TFBS(=N) rows */
  while (i < N) {
    unsigned int bitCount = countBits(i); 
    unsigned int zeros;
    int lastpos = -1;

    printBinaryRepresentation(i);
    printf(" (bits = %d) ", bitCount);

    /* there are TFBS-bitCount zero entries, create a new binary number for each */
    for (zeros = 0; zeros < TFBS-bitCount; zeros++) {
      unsigned int konpos;
      int b;
      konpos = i;

      /* go through the binary number bitwise */
      for (b = 0; b < TFBS; b++) {
	/* get the last untoggled bit location where we have a zero */
	if ((getBit(konpos, b) == 0) && b > lastpos) {
	  lastpos = b;
	  /* printf(" bitVal: 0 at %d ", b); */
	  konpos = setBit(konpos, b);
	  break;
	}
      }
      printf("[");    
      printBinaryRepresentation(konpos);
      printf("] ");    
    }
    printf("\n");
    i++;
  }
}
