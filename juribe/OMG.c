#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TFBS 2
#define DIM 4

static float Tmatrix[4][4]={{1, 2, 3, 0}, {4, 5, 0, 0}, {7, 0, 8, 0}, {0, 0, 0, 0}};                            
static float Xvector[4] = {0, 1, 1, 1};

static float Fmatrix[8][8]={{.2, .2, .1, 0, .4, 0, 0, 0}, {.1, .2, 0, .3, 0, .2, 0, 0}, {.2, 0, .2, .1, 0, 0, .3, 0}, {0, .2, .2, .1, 0, 0, 0, .2},
                            {.2, 0, 0, 0, .3, .2, .1, 0}, {0, .2, 0, 0, .3, .2, 0, .2}, {0, 0, .2, 0, .2, 0, .3, .2}, {0, 0, 0, .2, 0, .2, .1, .4}};                            
static float Fvector[8] = {.2, 0, .1, .3, .2, 0, .2, .1};


struct Rowtype {
    int rownum;
    float kval;
};
    
struct Ttype {
  int col;
  struct Rowtype *row;
  int rowCount;
};  


void compressT(float matrixT[DIM][DIM], struct Ttype *arrayT, int *colCount)
{
    int i, j, n, m, check, zeros;
     
     for(i=0; i<DIM; i++){  //this just prints the original T matrix
         for(j=0; j<DIM; j++){
             printf("%.1f ", matrixT[i][j]);
         }
         printf("\n");
     } 
     printf("\n");   
       
    n=0;  //n, m move through the Ttype array and Rowtype array
    i=0;  //i, j move through T matrix
    while (i<pow(2,TFBS)) {
        zeros=0;
        for (check=0; check<pow(2,TFBS); check++) {
            if(matrixT[check][i] != 0)
                zeros++;
        }
        while(zeros == 0){
           i++;
           for(check=0; check<pow(2,TFBS); check++) {
              if(matrixT[check][i] != 0)
                zeros++;
           }
        }                
       if(i<pow(2,TFBS)){
         arrayT[n].col=i;
         arrayT[n].row = malloc((zeros+1)*sizeof(struct Rowtype));

         m=0;  //m moves through the Rowtype array inside the Ttype structure
         for(j=0; j<pow(2,TFBS); j++) { //j moves through the rows of the T array
             if(matrixT[j][i] != 0) {
                arrayT[n].row[m].rownum = j;  //put in row number
                arrayT[n].row[m].kval = matrixT[j][i];  //put in value at that place
                m++;  //go to next element in Rowtype array
             }
         }
         arrayT[n].rowCount = m;
         n++;  //if the entire column is non-zero, then go to next place in Ttype array
         i++;
       }
   }    
   *colCount = n;
} 
 
void multiplyT(struct Ttype *arrayT, float *newX, int colCount)
{  
    int c, n, j;  
    float sum;
    float storeX[DIM];
   
    for(j=0; j<DIM; j++){  //copy and store newX array and initialize newX array to 0
        storeX[j] = newX[j];
        newX[j] = 0;
    }
    
    c=0;  //goes through Tarray
    while (c < colCount) {
       sum=0;  // next part goes through Rowtype array
       n = 0;
       while (n < arrayT[c].rowCount) {
         sum += storeX[(arrayT[c].row[n].rownum)]* arrayT[c].row[n].kval;/*multiplies element in storeX with corresponding element in SmallT then adds*/
         n++;
       }
      newX[(arrayT[c].col)] = sum;  //store sum of each column in corresponding newX place
      c++;
    } 
}    

int main(){
    int d;
    int colCount;
    struct Ttype *arrayT;
    float *newX;
    newX = calloc(pow(2,TFBS), sizeof(float));
    arrayT = malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));//=4 * sizeof()
    
    for(d=0; d<pow(2,TFBS); d++){
         newX[d] = Xvector[d];
     }    
    
    compressT(Tmatrix, arrayT, &colCount);
   
    int p,q;
    p=0;
    while (p < colCount) {
       printf( "Column: %d\n", arrayT[p].col); 
       q=0;
       while(q < arrayT[p].rowCount) {
    	  printf( "Row%d: %d\n",q, arrayT[p].row[q].rownum);
	      printf( "Value%d: %.2f\n",q, arrayT[p].row[q].kval);     
    	  q++;
       }
       p++; 
       printf("\n");   
    }
    printf("\n");
    
    multiplyT(arrayT, newX, colCount);
    int f;
    for(f=0; f<pow(2,TFBS); f++){
        printf("%.2f ", newX[f]); 
    }
    printf("\n"); 
    
    multiplyT(arrayT, newX, colCount);
    multiplyT(arrayT, newX, colCount);
    
    for(f=0; f<pow(2,TFBS); f++){
        printf("%.2f ", newX[f]); 
    }
    printf("\n");  
   
    int j;
    for (j=0; j<colCount; j++){
        free(arrayT[j].row);
    }    
    free(arrayT);
    free(newX);
                                         
    system("PAUSE"); 
}
      
