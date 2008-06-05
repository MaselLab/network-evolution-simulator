#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* THIS IS THE ONE YOU NEED TO WORK ON!!!*/
/*Fix compressT.. in another file*/
#define TFBS 2
#define DIM 4

static int Tmatrix[4][4]={{1, 2, 0, 0}, {4, 5, 0, 6}, {0, 0, 0, 9},
                            {0, 10, 0, 12}};
                            
static double Xvector[4] = {0, 1, 1, 1};

struct Rowtype{
    int rownum;
    float kval;
};
    
struct Ttype{
    int col;
    struct Rowtype row[(TFBS+1)];//TFBS+1 is the maximum # of non zero entries in any column
    //does the row have to be a pointer?? or a double pointer **??
};  
         

void compressT(int matrixT[DIM][DIM], struct Ttype **SmallT)
{
    int i, j, k, m, check, zeros;
    
     
     for(i=0; i<DIM; i++){//this just prints the original T matrix
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
        }}   
           
                
                          
       printf("zeros= %d\n", zeros);
       printf("i= %d\n", i); 
       printf("k: %d\n", k);   
       SmallT[k] = malloc(zeros*sizeof(struct Rowtype)+1);
       SmallT[k]->col=i;

       m=0;  //m moves through the Rowtype array inside the Ttype structure
       for(j=0; j<pow(2,TFBS); j++) { //j moves through the rows of the T array
           if(Tmatrix[j][i] != 0) {
              SmallT[k]->row[m].rownum=j;  //put in row number
              SmallT[k]->row[m].kval=Tmatrix[j][i];  //put in value at that place
              m++;
              //printf("%d\n",SmallT[k]->row[m].rownum);   //go to next element in Rowtype array
           }
       }
       
       //if the entire column is non-zero, then go to next place in Ttype array
       k++;
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
   
    for(j=0; j<DIM; j++){//copy and store newX array and initialize newX array to 0
        storeX[j]=newX[j];
        newX[j]=0;
    }
    //goes through Tarray
    c=0;
   while (SmallT[c]!=NULL){/*number of elements in Tarray, k-- for some reason this is an issue...have
   to adjust based on number of significant rows*/
       sum=0;
       //goes through Rowtype array
       for(n=0; n<3; n++){/*3 is the number of non-zero elements in each column, which is m but its
       specific to each column so changes based on column #*/
           sum+= storeX[(SmallT[c]->row[n].rownum)]* SmallT[c]->row[n].kval;/*multiplies element in
           storeX with corresponding element in SmallT then adds*/
       }
       printf("%d %.2f\n", (SmallT[c]->col), sum);
       newX[(SmallT[c]->col)]=sum;//store sum of each column in corresponding newX place
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
    

    SmallT=malloc(pow(2,TFBS)*sizeof(struct Ttype)+1);//=4 * sizeof()

    compressT(Tmatrix, SmallT);
    //printf("Here");
    
   /* int count, s;//gives the number of significant elements in T-array- same as k in compressT
    count=1;
    s=1;
    while(SmallT[s]->col!=0){
        count++;
        s++;
    }
    printf("count: %d\n", count);*/
    // printf("here2"); 
    //SmallT=realloc(SmallT, sizeof(struct Ttype)*(count+2));//works bc same as malloc= 4* sizeof()
    /*SmallT=realloc(SmallT, sizeof(struct Ttype)*(count)); should work but doesnt...
    does realloc only work to make memeory larger??
    /* realloc?? not sure how to reallocate memory to make it smaller. i think in the compressT
    function its setting all unused positions to 0 and places them at end of array, so when I try
    to realloc(SmallT) it does not actually make the chunk of memory allocated to it smaller.*/
   
    int p,q;
    for (p=0; p<3; p++){
        printf( "Column: %d\n", SmallT[p]->col); 
        for(q=0;q<3;q++){
           printf( "Row%d: %d\n",q, SmallT[p]->row[q].rownum);
           printf( "Value%d: %.2f\n",q, SmallT[p]->row[q].kval);     
        }    
    }
   // printf("Here again"); 
  /* multiplyT(SmallT, newX);
   
   int f;
   for(f=0; f<4; f++){
       printf("%.2f ", newX[f]); 
   } 
   printf("\n"); 
   multiplyT(SmallT, newX);
   
   for(f=0; f<4; f++){
       printf("%.2f ", newX[f]); 
   } 
   printf("\n"); */
   

    free(SmallT);
    free(newX);
    
                                         
    system("PAUSE"); }
      
