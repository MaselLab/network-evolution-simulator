#include <stdio.h>
#include <math.h>

#define TFBS 3

/*this should print all the non-zero Kon values in each column based on the given 
hinderances array*/

static unsigned int hinderances[TFBS] = {2, 1, 0};

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
     //printf("%.*d (%*d)", TFBS, binary, TFBS-1, decimal);
     printf("%d ", decimal); 
}

//this function checks if there are consecutive ones
int checkConfigThree(int num) {
    while (num!=3 && num!= 0) {
        num=(num >> 1);
    }
    if (num==3)
       return 1;
    else
       return 0;
}

//this function makes an array of the impossible states 
void makeImpStates(int BS, unsigned int hind[TFBS], unsigned int check[TFBS], int *size) {
     int num[BS];
     int i,j;
     j=0;
     
     for (i=0; i<BS-1; i++) {
        num[i]= hind[i] | hind[i+1];
     }
 
     for (i=0; i<BS-1; i++){ 
        if (checkConfigThree(num[i])) {
           check[j]=num[i];
           j++;
        }
     }
      *size=j;
}


int main(int argc, char *argv[])
{
  unsigned int N = pow(2, TFBS);
  unsigned int col = 0;

  /* loop through all 2^TFBS(=N) columns */
  
  //determine impossible states
  unsigned int p;
  unsigned int check[TFBS];
  int j, b;
  j=0;
  makeImpStates(TFBS, hinderances, check, &j);
  
  while (col < N) {
    b=0;
    //check that starting state is possible
    while (b<j && ((check[b] & col) != check[b])) {
       b++;
    }
    //if starting state is valid
    if(b==j){
      printBinaryRepresentation(col);
    
    unsigned int konpos;
    unsigned int hinderance_mask;
    
    for (p = 0; p < TFBS; p++)
      if ((col & (1 << p))) {
    	konpos = col ^ (1 << p) ;  /* then add a one there, effectively doing what the old setBit function did */
        hinderance_mask = hinderances[p];
        if (!(konpos & hinderance_mask)) {
           printf("[");    
	       printBinaryRepresentation(konpos);
	       printf("] ");    
        }	
      }
    printf("\n");
    }
    //if starting state not valid
   else{
      printBinaryRepresentation(col); 
      printf("\n");
    }
    col++;
  }    
  system("PAUSE");
}
