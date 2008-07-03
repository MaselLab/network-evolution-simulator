#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TFBS 3
#define DIM 8
#define SIZE 0 //number of elements in hinderance array

static float kon[TFBS] = {0.1, 0.4, 0.3};
static float koff[5] = {0.15, 0.25, 0.35, 0.45, .55};
static float diag[DIM] = {0,0,0,0,0,0,0,0};
static int hinderances[TFBS][2] = {};//no hinderance

static float Xvector[DIM] = {0, .3, .3, .2, 0, .4, 0, 0};


struct Rowtype {
    int rownum;
    float *kval;
};
    
struct Ttype {
  int col;
  struct Rowtype *row;
  int rowCount;
};

/* this returns a 0 if state e is not possible and 1 if state e is
   possible based on hind array. This is a new function from last time
   to deal with all special cases of more than two binding sites
   hindering each other. The hinderance array holds ordered pairs
   ie. hinderances[TFBS][2] ={{0,2}, {1,2}} where the first number
   represents the first binding site that is hindered and the second
   number is the number of binding sites involved in the hinderance.
   So {0,2} means 0 is hindered and that there are two binding sites
   involved, 0 and 1. So the function just checks to make sure that
   only ONE of the number of involved positions is bound at one
   time. SIZE is the number of hinderance occurences. */
int ifPossible(unsigned int e, int hind[TFBS][2]){
    int possible, i, m , b, d;
    unsigned int s, w;
    possible = 1;
    i = SIZE-1;
   
    if (e>0) {
       while (possible==1 && i>=0) {
          d = hind[i][1];
          s = (e >> (TFBS - hind[i][0] - hind[i][1]));
          int p= (int) pow(2,d);
          w = s % p;
          if (w>0) {
            m=d;
            b=(1 << d);
            while (m>=0 && w!=b) {
               m--;
               b = b/2;
            }
            if (m<0) {
               possible = 0;
            }
          }
          i--;
       }
    }
    return (possible);   
}

//populates the diagonal entries
void diagonal(int col, float diag[DIM], struct Ttype *arrayT, int m, int n){
    int x;
    diag[col]=1;
          for (x=0; x<m; x++) {
             diag[col] -= *(arrayT[n].row[x].kval);
          }
      arrayT[n].row[m].rownum = col;
      arrayT[n].row[m].kval = &(diag[col]);   
} 
     
//creates the Tarray structure
void createStruct(int TFBSites,  struct Ttype *arrayT, float kon[TFBS], float koff[5], int *colCount, int hind[TFBS][2], float s[DIM]){
  
  unsigned int N = pow(2, TFBSites);
  unsigned int col = 0;
  
  unsigned int check[TFBS];
  int j, b, n, m, x;
 
  j=0;
  n=0;
 
  /* loop through all 2^TFBS(=N) columns */
  while (col < N) {
 
     if (ifPossible(col, hind)==1) {
        
        arrayT[n].col = col;
        arrayT[n].row = malloc((TFBS+2)*sizeof(struct Rowtype));
          
       int p;
       unsigned int rowOn, rowOff, row;

       m=0;
       for (p = 0; p < TFBSites; p++) {
          if ((col & (1 << p))){
              rowOn = col ^ (1 << p);              
              if (ifPossible(rowOn, hind)==1) {
                 arrayT[n].row[m].rownum = rowOn;
                 arrayT[n].row[m].kval = &(kon[p]);
                 m++;  
              }	
          } else  {
              rowOff = col | (1 << p) ; 
              if (ifPossible(rowOff, hind)==1) {
                 arrayT[n].row[m].rownum = rowOff;
                 arrayT[n].row[m].kval = &(koff[p]);
                 m++;   
              }	 
          } 
       }
      diagonal(col, s, arrayT, m, n);
      m++;
      arrayT[n].rowCount = m;
      n++;                                
   }
    col++;
  }       
  *colCount = n;       
} 

//multiplies by Xvector
void multiplyT(struct Ttype *arrayT, float *newX, int colCount)
{  
    int c, n, j;  
    float sum, k;
    float storeX[DIM];
   
    for (j=0; j<DIM; j++) {  //copy and store newX array and initialize newX array to 0
        storeX[j] = newX[j];
        newX[j] = 0;
    }
    
    c=0;  //goes through Tarray
    while (c < colCount) {
       sum=0;  // next part goes through Rowtype array
       n = 0;
       while (n < arrayT[c].rowCount) {
         k = *(arrayT[c].row[n].kval);
         sum += storeX[(arrayT[c].row[n].rownum)]* k ;/*multiplies element in storeX with corresponding element in SmallT then adds*/
         n++;
       }
      newX[(arrayT[c].col)] = sum;  //store sum of each column in corresponding newX place
      c++;
    } 
} 

//prints the Tarray
void print_arrayT(struct Ttype *arrayT, int colCount){
     int p, q;  
     
     p=0;
    while (p < colCount) {
       printf( "Column: %d\n", arrayT[p].col); 
       q=0;
       while (q < arrayT[p].rowCount) {
    	  printf( "Row%d: %d\n",q, arrayT[p].row[q].rownum);
	      printf( "Value%d: %.2f\n",q, *arrayT[p].row[q].kval);     
    	  q++;
       }
       p++; 
       printf("\n");   
    }
    printf("\n");
} 

//prints the Xvector
void print_vector(float *newX, int length){
    int d; 
    for (d=0; d<length; d++) {
        printf("%.2f ", newX[d]); 
    }
      printf("\n"); 
      printf("\n"); 
}     


int main(int argc, char *argv[])
{
    int d, colCount, p, q; 
    struct Ttype *arrayT;
    float *newX;
    
    newX = calloc(pow(2,TFBS), sizeof(float));
    arrayT = malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));
    
    for (d=0; d<pow(2,TFBS); d++) {
         newX[d] = Xvector[d];
     }    
    
  createStruct(TFBS, arrayT, kon, koff, &colCount, hinderances, diag);
  
  printf("Initial Tarray\n");
  print_arrayT(arrayT, colCount);
  
  kon[0] = 0.36;
  kon[1] = 0.26;
  kon[2] = 0.36;
  
  for (d=0; d<colCount;d++) { 
     diagonal(arrayT[d].col, diag ,arrayT, (arrayT[d].rowCount-1), d);
 }    
     
  printf("Changed Tarray\n");
  print_arrayT(arrayT, colCount); 
  
  printf("Original Xvector\n");
  print_vector(newX, pow(2, TFBS));
  
  for (d=0; d<40; d++) { 
  multiplyT(arrayT, newX, colCount);
  } 
  
  printf("Converged Xvector\n"); 
  print_vector(newX, pow(2, TFBS));

  for (d=0; d<colCount; d++) {
       free(arrayT[d].row);
  }    
  free(arrayT);
  free(newX);
  
  system("PAUSE");
}
