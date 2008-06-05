#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* THIS IS THE ONE YOU NEED TO WORK ON!!!*/
/*Fix compressT.. in another file*/
#define TFBS 2
#define DIM 4

static int Tmatrix[4][4]={{1, 2, 0, 0}, {4, 5, 0, 6}, {7, 0, 0, 9},
                            {0, 10, 0, 12}};
                            
static double Xvector[4] = {0, 1, 1, 1};

struct Rowtype{
    int rownum;
    float kval;
};
    
struct Ttype{
    int col;
    struct Rowtype **row;//TFBS+1 is the maximum # of non zero entries in any column
    //does the row have to be a pointer?? or a double pointer **??
};  
         

void compressT(int matrixT[DIM][DIM], struct Ttype **SmallT)
{
    int i, j, k, m, check, zeros;
    
     
     for(i=0; i<DIM; i++){  //this just prints the original T matrix
         for(j=0; j<DIM; j++){
             printf("%d ", matrixT[i][j]);
         }
         printf("\n");
     } 
     printf("\n");   
       
   
    k=0;//k moves through the Ttype array
    i=0;//i moves through the columns of the T array
    while (i<pow(2,TFBS)) {
        printf("i= %d\n", i);
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
       printf("zeros= %d\n", zeros);
       printf("i= %d\n", i); 
       printf("k: %d\n", k); 
        
       printf("before SmallT[k] malloc\n");
       system("PAUSE");
       SmallT[k]->row=malloc((zeros+1)*sizeof(struct Rowtype));
      
       SmallT[k]->col=i;

       m=0;  //m moves through the Rowtype array inside the Ttype structure
       for(j=0; j<pow(2,TFBS); j++) { //j moves through the rows of the T array
           printf("i=%d, j=%d, Tmatrix[j][i]=%d, k=%d, m=%d\n", i, j, Tmatrix[j][i], k, m);
           if(Tmatrix[j][i] != 0) {
              system("PAUSE");
              SmallT[k]->row[m]->rownum = j;  //put in row number
              SmallT[k]->row[m]->kval=Tmatrix[j][i];  //put in value at that place
              printf("after row[m] assignment, k=%d, m=%d, smallT row=%d, kval=%g\n", 
                k, m, SmallT[k]->row[m]->rownum, SmallT[k]->row[m]->kval);
              m++;  //go to next element in Rowtype array
             
              system("PAUSE");
           }
       }
       printf("before NULL ptr assignment\n");
       system("PAUSE");
       //SmallT[k]->row[m]= NULL;
       printf("after NULL ptr assignment\n");
       //Need to assign NULL pointer at the end of the Rowtype array... not sure how to do this
       //SmallT[k]->row[m]=NULL;
       
       k++;  //if the entire column is non-zero, then go to next place in Ttype array
       i++;
       printf("SmallT, i = %d, k = %d\n", i, k);
    }
    SmallT[k] = NULL; 
    printf("k: %d\n", k); 
    printf("m: %d\n", m);
    printf("\n");
 
} 
 
void multiplyT(struct Ttype **SmallT, float *newX)
{  
    int c, n, j;  
    float sum;
    int storeX[DIM];
   
    for(j=0; j<DIM; j++){  //copy and store newX array and initialize newX array to 0
        storeX[j]=newX[j];
        newX[j]=0;
    }
    
   c=0;  //goes through Tarray
   while (SmallT[c]!=NULL){  /*goes through Tarray of structures (outer)*/
       sum=0;  // next part goes through Rowtype array
       for(n=0; n<3; n++){  /* number of non-zero elements in each column, should be NULL while*/
           sum+= storeX[(SmallT[c]->row[n]->rownum)]* SmallT[c]->row[n]->kval;/*multiplies element in storeX with corresponding element in SmallT then adds*/
       }
       printf("%d %.2f\n", (SmallT[c]->col), sum);
       newX[(SmallT[c]->col)]=sum;  //store sum of each column in corresponding newX place
       c++;
   } 
}    

int main(){
    int d;
    struct Ttype **SmallT;
    float *newX;
    newX=calloc(4, sizeof(float));
    
    //how do we intialize newX to Xvector without specifically assigning each element...
    for(d=0; d<4; d++){
         newX[d]=Xvector[d];
     }    

    SmallT=malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));//=4 * sizeof()

    compressT(Tmatrix, SmallT);
   
    int p,q;
    p=0;
    while (SmallT[p]!=NULL){
        printf( "Column: %d\n", SmallT[p]->col); 
        for(q=0;q<3;q++){
           printf( "Row%d: %d\n",q, SmallT[p]->row[q]->rownum);
           printf( "Value%d: %.2f\n",q, SmallT[p]->row[q]->kval);     
        }
        p++;    
    }
  // printf("Here again"); 
  // printf("\n");
   multiplyT(SmallT, newX);
   
   int f;
   for(f=0; f<4; f++){
       printf("%.2f ", newX[f]); 
   } 
  /* printf("\n"); 
   multiplyT(SmallT, newX);
   
   for(f=0; f<4; f++){
       printf("%.2f ", newX[f]); 
   } 
   printf("\n"); */
   

    free(SmallT);
    free(newX);
    
                                         
    system("PAUSE"); }
      
