#include <stdio.h>
#include <math.h>

//#define TFBS 5
#define MAX_TFBS 100

/* this simply prints the binary representation */
void printBinaryRepresentation(int decimal, unsigned int TFBS) {
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

/* the following 3 functions aren't used any more as they have been absorbed 
 * into the main() program */

/* counts the number of "ones" in the binary represesentation of the number */
int countBits(unsigned int v) {

  unsigned int c; // c accumulates the total bits set in v
  for (c = 0; v; c++) {
      v &= v - 1; // clear the least significant bit set
  }
  return (c);
}

/* this gets the bit at position pos */
int getBit(unsigned int val, int pos) {
  return (val >> pos) & 1; // >> means val shifts right pos times. Why the & 1??
}

/* this set the bit at position pos */
int setBit(unsigned int val, int pos) {
  unsigned int mask = 1 << pos;// << means 1 shifts left pos times
  return mask | val; // not sure why it returns mask or val, is val changed at all?
}

int main(int argc, char *argv[])
{
  int debug = 0;
  unsigned long int TFBS;

  for (TFBS=0; TFBS < MAX_TFBS; TFBS++) {

  unsigned long int N = pow(2, TFBS);
  unsigned long int N_squared = N*N;
  unsigned long int i = 0;
  unsigned long int total_kon = 0;

  /* loop through all 2^TFBS(=N) rows */
  while (i < N) {
    unsigned int bitCount = countBits(i); 
    unsigned int p;
    
    if (debug) {
      printBinaryRepresentation(i, TFBS);
      printf(" (bits = %d) ", bitCount);
    }
    unsigned long int konpos;
    
    /* there are at most TFBS possible states accessible from current state
       since we only bind at most one new TF, go through and check if it is
       valid */
    for (p = 0; p < TFBS; p++) 
      /* if there is not already a 1 in the p-position */
      if(!(i & (1 << p))) {
	konpos = i | (1 << p) ;  /* then add a one there, effectively doing what the old setBit function did */
	if (debug) {
	  printf("[");    
	  printBinaryRepresentation(konpos, TFBS);
	  printf("] ");   
	} 
	total_kon++;
      }
    if (debug) {
      printf("\n");
    }
    i++;
  }
  printf("%lu %lu %lu %g\n", TFBS, total_kon, N_squared, (double)total_kon/N_squared);
  fflush(stdout);
  //system("date");
  }
  //system("PAUSE");
}

/* 
   Local Variables:
   compile-command: "gcc -o count-kons count-kons.c -lm"
   End: 
*/
