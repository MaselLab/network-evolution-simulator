#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TFBS 2
#define DIM 4

static int Tmatrix[4][4]={{1, 2, 0, 0}, {4, 5, 0, 6}, {7, 0, 0, 9},
                            {0, 10, 0, 12}};
                            
static double Xvector[4] = {0, 1, 1, 1};

struct Rowtype {
    int rownum;
    float kval;
};
    
struct Ttype {
  int col;
  struct Rowtype *row;
  int rowCount;
};  
         

void compressT(int matrixT[DIM][DIM], struct Ttype *SmallT, int *colCount)
{
    int i, j, n, m, check, zeros;
    
     
     for(i=0; i<DIM; i++){  //this just prints the original T matrix
         for(j=0; j<DIM; j++){
             printf("%d ", matrixT[i][j]);
         }
         printf("\n");
     } 
     printf("\n");   
       
   
    n=0;//n, m move through the Ttype array and Rowtype array
    i=0;//i, j move through T matrix
    while (i<pow(2,TFBS)) {
        zeros=0;
        for(check=0; check<pow(2,TFBS); check++) {
            if(Tmatrix[check][i]!=0)
                zeros++;
        }
        if(zeros==0){
           i++;
           for(check=0; check<pow(2,TFBS); check++) {
              if(Tmatrix[check][i]!=0)
                zeros++;
           }
        }                
      
       SmallT[n].col=i;
       SmallT[n].row=malloc((zeros+1)*sizeof(struct Rowtype));

       m=0;  //m moves through the Rowtype array inside the Ttype structure
       for(j=0; j<pow(2,TFBS); j++) { //j moves through the rows of the T array
           if(Tmatrix[j][i] != 0) {
              SmallT[n].row[m].rownum = j;  //put in row number
              SmallT[n].row[m].kval=Tmatrix[j][i];  //put in value at that place
              m++;  //go to next element in Rowtype array
           }
       }
       SmallT[n].rowCount = m;
       
       n++;  //if the entire column is non-zero, then go to next place in Ttype array
       i++;
    }
    *colCount = n;
} 
 
void multiplyT(struct Ttype *SmallT, float *newX, int colCount)
{  
    int c, n, j;  
    float sum;
    int storeX[DIM];
   
    for(j=0; j<DIM; j++){  //copy and store newX array and initialize newX array to 0
        storeX[j]=newX[j];
        newX[j]=0;
    }
    
    c=0;  //goes through Tarray
    while (c < colCount) {
       sum=0;  // next part goes through Rowtype array
       n = 0;
       while (n < SmallT[c].rowCount) {
         sum+= storeX[(SmallT[c].row[n].rownum)]* SmallT[c].row[n].kval;/*multiplies element in storeX with corresponding element in SmallT then adds*/
         n++;
       }
      newX[(SmallT[c].col)]=sum;  //store sum of each column in corresponding newX place
      c++;
    } 
}    

int main(){
    int d;
    int colCount;
    struct Ttype *SmallT;
    float *newX;
    newX=calloc(4, sizeof(float));
    
    for(d=0; d<4; d++){
         newX[d]=Xvector[d];
     }    

    SmallT=malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));//=4 * sizeof()

    compressT(Tmatrix, SmallT, &colCount);
   
    int p,q;
    p=0;
    while (p < colCount) {
       printf( "Column: %d\n", SmallT[p].col); 
       q=0;
       while(q < SmallT[p].rowCount) {
    	  printf( "Row%d: %d\n",q, SmallT[p].row[q].rownum);
	      printf( "Value%d: %.2f\n",q, SmallT[p].row[q].kval);     
    	  q++;
       }
       p++;    
    }
    printf("\n");
    
    multiplyT(SmallT, newX, colCount);
    int f;
    for(f=0; f<4; f++){
        printf("%.2f ", newX[f]); 
    }
    printf("\n"); 
    
    multiplyT(SmallT, newX, colCount);
    for(f=0; f<4; f++){
    printf("%.2f ", newX[f]); 
    }
    printf("\n");  
   
    free(SmallT);
    free(newX);
                                         
    system("PAUSE"); 
}
      
