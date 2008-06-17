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
/*counts the number of "ones" in the binary represesentation of the number*/
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
  unsigned int mask = 1 << pos;/* << means 0 shifts left pos times--NOT WORKING.
                                    for Koff values. Ones not being changed to 0.
  printf("val: %d",mask); */
  return mask ^ val;   /* here you need to keep the mask 1, but use XOR (exclusive OR) */
}

int main(int argc, char *argv[])
{
  unsigned int N = pow(2, TFBS);
  unsigned int i = 0;

  /* loop through all 2^TFBS(=N) rows */
  while (i < N) {
    unsigned int bitCount = countBits(i); 
    unsigned int b;

    printBinaryRepresentation(i);
    printf(" (bits = %d) ", bitCount);
    unsigned int koffpos;

 
      for (b = 0; b < TFBS; b++) 
          if((i & (1 << b))){
              koffpos = i ^ (1 << b);//inclusive or exclusive??
              printf("[");    
              printBinaryRepresentation(koffpos);
              printf("] ");    
           }
           printf("\n");
           i++;
  }
  system("PAUSE");
}
