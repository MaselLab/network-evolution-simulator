#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TFBS 7
#define DIM 12
#define SIZE 2 //number of elements in hinderance array

static float kon[TFBS] = {1, 2, 3,4,5,6,7};
static float koff[3] = {1.5, 2.5, 3.5};//need to add cooperativity--does it matter left or right? 
static float hammDist[7] = {2,2,2,2,2,2,1};
static float diag[DIM] = {0,0,0,0,0,0,0,0,0,0,0,0};
static int hinderances[TFBS][2] = {{0,5},{2,5}};

static float Xvector[DIM] = {0, .3, .3, .2, 0, .4, 0, 0,1,1,1,1};

FILE *sparseMatrix1;


struct Coltype {
    int colnum;
    float *kval;
};
    
struct Ttype {
  int row;
  struct Coltype *col;
  int colCount;
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
void diagonal(int row, float diag[DIM], struct Ttype *arrayT, int m, int n){
    int x;
    printf("diag[%d]=%d\n", row, diag[row]);
    diag[row]=0;
          for (x=0; x<m; x++) {
             diag[row] -= *(arrayT[n].col[x].kval);
          }
      arrayT[n].col[m].colnum = row;
      arrayT[n].col[m].kval = &(diag[row]); 
      printf("arrayT[n].col[m].kval=%.2f\n",*arrayT[n].col[m].kval);  
} 
     
//creates the Tarray structure
void createStruct(int TFBSites,  struct Ttype *arrayT, float kon[TFBS], float koff[3], int *rowCount, int hind[TFBS][2], float s[DIM]){
  
  unsigned int N = pow(2, TFBSites);
  unsigned int row = 0;
  
  unsigned int check[TFBS];
  int j, b, n, m, x;
 
  j=0;
  n=0;
 
  /* loop through all 2^TFBS(=N) rows */
  while (row < N) {
 
     if (ifPossible(row, hind)==1) {
        
        arrayT[n].row = row;
        arrayT[n].col = malloc((TFBS+2)*sizeof(struct Coltype));
          
       int p;
       unsigned int colOn, colOff, col;

       m=0;
       for (p = 0; p < TFBSites; p++) {
          if (!(row & (1 << p))){
              colOn = row | (1 << p);              
              if (ifPossible(colOn, hind)==1) {
                 arrayT[n].col[m].colnum = colOn;
                 arrayT[n].col[m].kval = &(kon[p]);
                 m++;  
              }	
          } else  {
              colOff = row ^ (1 << p) ; 
              if (ifPossible(colOff, hind)==1) {
                 arrayT[n].col[m].colnum = colOff;
                 int hamDist;
                 hamDist = hammDist[(TFBS-1)-p];
                 //printf("hamDist = %d\n", hamDist);
                 arrayT[n].col[m].kval = &(koff[hamDist]);
                 m++;   
              }	 
          } 
       }
      diagonal(row, s, arrayT, m, n);
      m++;
      arrayT[n].colCount = m;
      n++;                                
   }
    row++;
  }       
  *rowCount = n; 
} 


//prints the Tarray
void print_arrayT(struct Ttype *arrayT, int rowCount){
   
     int p, q;  
     p=0;
     
     printf( "%d,   %d,   %.2f\n\n", arrayT[6].row, arrayT[6].col[3].colnum, *arrayT[6].col[3].kval);
     double fine;
     fine = *arrayT[6].col[3].kval;
    while (p < rowCount) {
       printf( "Row%d: %d\n",p, arrayT[p].row); 
       q=0;
       while (q < arrayT[p].colCount) {
    	  printf( "Column%d: %d\n",q, arrayT[p].col[q].colnum);
	      printf( "Value%d: %.2f\n",q, *arrayT[p].col[q].kval);     
    	  q++;
       }
       p++; 
       printf("\n");   
    }
    printf("\n");
     
    // system("PAUSE");
     p=0;
     if (sparseMatrix1 = fopen("sparseMatrix1.txt", "w")){
        //ask alex about this line of code...
       fprintf(sparseMatrix1, "1Row  Column  Value\n\n");
       //printf( "%d,   %d,   %.2f\n\n", arrayT[6].row, arrayT[6].col[3].colnum, *arrayT[6].col[3].kval);
       fprintf(sparseMatrix1, "%d,   %d,   %.2f\n\n", arrayT[6].row, arrayT[6].col[3].colnum, *arrayT[6].col[3].kval);
   // system("PAUSE");
    
     while (p < rowCount) {
       q=0;
       while (q < arrayT[p].colCount) { 
          fprintf(sparseMatrix1, " %d  ", arrayT[p].row+1);
    	  fprintf(sparseMatrix1, " %d    ", arrayT[p].col[q].colnum+1);
	      fprintf(sparseMatrix1, " %.2f \n", *arrayT[p].col[q].kval);  
          //printf( "\n\nValue%d: %.2f\n",q, *arrayT[p].col[q].kval);    
    	  q++;
  	  
       }
       p++; 
       fprintf(sparseMatrix1, "\n");   
    }
    fprintf(sparseMatrix1, "\n");
    //fclose(sparseMatrix1);
  }else{
    printf("Error printing to File");
  }  
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
    int d, rowCount, p, q; 
    struct Ttype *arrayT;
    float *newX;
    
    
    newX = calloc(pow(2,TFBS), sizeof(float));
    arrayT = malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));
    rowCount=0;
    //sparseMatrix1 = fopen("sparseMatrix1.txt", "w");
     
  createStruct(TFBS, arrayT, kon, koff, &rowCount, hinderances, diag);
  
  printf("Initial Tarray\n");
  //system("PAUSE");
  //printf( "%d,   %d,   %.2f\n\n", arrayT[6].row, arrayT[6].col[3].colnum, *arrayT[6].col[3].kval);
  print_arrayT(arrayT, rowCount);
  

  fclose(sparseMatrix1);
  
  for (d=0; d<rowCount; d++) {
       free(arrayT[d].col);
  }    
  free(arrayT);
  free(newX);
  
  system("PAUSE");
}
