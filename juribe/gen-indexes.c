#include <stdio.h>
#include <math.h>

#define TFBS 2
#define DIM 4

static float *matrix[DIM][DIM] = {{NULL, NULL, NULL, NULL}, 
			 {NULL, NULL, NULL, NULL}, 
			 {NULL, NULL, NULL, NULL}, 
			 {NULL, NULL, NULL, NULL}};
static float kon[TFBS] = {0.1, 0.4};
static float koff[5] = {0.15, 0.25, 0.35, 0.45, .55};

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
     printf(" %d ", decimal);
     
}

void print_matrix(float *matrix[DIM][DIM]) {
  int i, j;
  for (i=0; i < DIM; i++) {
    for (j=0; j < DIM; j++) {
      if (matrix[i][j] != NULL)
	printf("%.4g ", *matrix[i][j]);
      else
	printf("NULL ");
    }
    printf("\n");
  }
}

void printKPos(int TFBSites, float *matrix[DIM][DIM], float kon[TFBS], float koff[5]){
  unsigned int N = pow(2, TFBSites);
  unsigned int col = 0;
  int i,j;

  /* loop through all 2^TFBS(=N) columns */
  while (col < N) {
     
    unsigned int p;
    printBinaryRepresentation(col);
    unsigned int rowOn;
     unsigned int rowOff;

    //simplify this code
    for (p = 0; p < TFBSites; p++) 
       if ((col & (1 << p))){
           rowOn = col ^ (1 << p);
           matrix[rowOn][col] = &(kon[p]);
           printf("[");    
           printBinaryRepresentation(rowOn);
           printf("] "); 
       } else {
           rowOff = col | (1 << p) ;  /* then add a one there, effectively doing what the old setBit function did */
    	   matrix[rowOff][col] = &(koff[p]);
           printf("[");    
           printBinaryRepresentation(rowOff);
           printf("] "); 
      } 
      printf("\n");
      col++;//make this incrementation conditional--inly increment when col not all zeros
  } 
      
}  


int main(int argc, char *argv[])
{
   printKPos(TFBS, matrix, kon, koff);
  
  printf("initial matrix:\n");
  print_matrix(matrix);

  /* modify kon */
  kon[0] = 0.0789;
  kon[1] = 0.9999;

  printf("\nmodified matrix:\n");

  /* this shows the modified kon values */
  print_matrix(matrix);
  
  kon[0] = 0.034;
  kon[1] = 0.333;
  
  printf("\nmodified matrix:\n");
  print_matrix(matrix);
  
  
  system("PAUSE");
}
