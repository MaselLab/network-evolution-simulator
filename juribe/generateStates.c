#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TFBS 14
#define SIZE 1
#define HIND_LENGTH 15

static int startPos[TFBS]= {1,2,3,3,10,10,11,11,11,13,16,16,17,17};
static float kon[TFBS] = {1, 2, 3,4,5,6,7,8,9,10,11,12,13,14};
static float koff[5] = {1.5, 2.5, 3.5, 4.5, 5.5};
static int hammDist[TFBS] = {2,2,2,2,2,2,1,2,1,2,2,2,1,2};
static float diag[TFBS]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//static int startPos[TFBS]= {1,2,3,3,10,10,11,11,11,13,16,16,17,17, 20,20,21,24,24,25,25,25};

/*static int startPos[TFBS]= {1,2,3,3,10,10,11,11,11,13,16,16,17,17,20,20,21,24,24,25,25,25,
                           26,26,26,27,28,28,29,32,36,37,37,41,44,46,48,51,53,53,57,57,57,
                           59,59,61,62,62,63,66,67,69,74,78,79,79,80,80,81,82,84,84,85,86,86,
                           88,88,89,90,90,91,91,92,92,93,94,100,100,100,103,103,104,104,110,
                           116,118,119,119,123,123,125,126,127};*/
                           
 struct Coltype {
    int colnum;
    float *kval;
};
    
struct Ttype {
  int row;
  struct Coltype *col;
  int colCount;
};                          
                           
int convertToDecimal(int bits[TFBS]){
     int n = TFBS-1;
     int count =0;
     int record[TFBS];
     int rec=0;
     while( n>=0 ){
         if(bits[n]==0){
            count++;
         }else{
             record[rec]=(TFBS-1)-n;
             rec++;
         }
         n--;
     }
     int i;
     int decimal=0;
     for(i=0; i<rec;i++){
        decimal+= pow(2,record[i]);
     }
     return decimal;
}  
 
                          
int isHindered(int bindSite, int bits[TFBS]){
    int s = bindSite - 1;
    int count=0;
    //printf("%d\n", s);
    //printf("%d\n", startPos[bindSite]);
    int check = startPos[bindSite]-HIND_LENGTH;
    if(check<0){
        check=0;}
    //printf("%d\n", check);
    while(startPos[s] <= startPos[bindSite] && startPos[s] > check){
       if(bits[s] == 1){
          return 1;
       } else {
          s--;
          count++;
       }
    }
    //printf("%d\n",count);
    return 0;
  }
  
  void configure(int bindSite, int bits[TFBS], int *numStates, int *statesArray){

     if(bindSite<TFBS-1){
        bits[bindSite] = 0;
        configure(bindSite+1,bits,numStates,statesArray);
        if(!isHindered(bindSite, bits)){
           bits[bindSite] = 1;
           configure(bindSite+1,bits,numStates,statesArray);
        }
     } else {
        bits[TFBS-1] = 0;
        statesArray[(*numStates)] = convertToDecimal(bits);
        (*numStates)++;
        int i;
        for(i=0; i<TFBS; i++){
           printf("%d", bits[i]);
        }
        printf("\n");
        if(!isHindered(bindSite, bits)){
           bits[TFBS-1]=1;
           convertToDecimal(bits);
           statesArray[(*numStates)] = convertToDecimal(bits);
           (*numStates)++;
           for(i=0; i<TFBS; i++){
              printf("%d", bits[i]);
           }
           printf("\n");
        }
      } 
      
  }
  
  void diagonal(int row, float diag[TFBS], struct Ttype *arrayT, int m, int n){
    int x;
    diag[row]=0;
          for (x=0; x<m; x++) {
             diag[row] -= *(arrayT[n].col[x].kval);
          }
      arrayT[n].col[m].colnum = row;
      arrayT[n].col[m].kval = &(diag[row]);   
} 
  
  void transitions(int size, int *viableStates, int TFBSites, struct Ttype *arrayT, float kon[TFBS], float koff[5],int hammDist[TFBS]){
       int i, p,j,m;
       int n=0;
       for( i=0;i<size; i++){
          arrayT[n].row = i;
          arrayT[n].col = malloc((TFBS+2)*sizeof(struct Coltype));
         //printf("viableStates:%d, row num:%d\n",viableStates[i], i);
        m=0;
        for(p=0;p<TFBSites;p++){
          int col = viableStates[i];
          
          if(col& (1<<p)){
            col = viableStates[i] ^ (1<<p);
           if( col!=viableStates[i]){
             for(j=0;j<size; j++){
               if(col==viableStates[j]){
                 arrayT[n].col[m].colnum = j;
                 int a = hammDist[p];
                 arrayT[n].col[m].kval = &(koff[a]);
                 m++;
                 // printf("    col=%d, j=%d, p=%d\n", col, j, p);
                // printf(" %d  %d  \n", i, j);
               }
             }
           }
         }else{
          
          int col = viableStates[i] | (1<<p);
          if(col!=0 && col!=viableStates[i]){
            for( j=0;j<size; j++){
              if(col==viableStates[j]){
                 arrayT[n].col[m].colnum = j;
                 arrayT[n].col[m].kval = &(kon[p]);
                 m++;
                 //printf("    col=%d, j=%d, p=%d\n", col, j, p);
                // printf(" %d  %d  \n", i, j);
               }
             }
           }
          }
       }
       diagonal(i,diag, arrayT, m, n);
       m++;
       arrayT[n].colCount = m;
       n++;
     }
 
  }
  
void decimalToBinary(int decimal) {
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
     printf("%d ", binary);    
}

void print_arrayT(struct Ttype *arrayT, int size){
     int p, q;  
    
    /*p=0; 
    printf("\n"); 
    while (p < size) {
           
       printf( "Row%d: %d\n",p, arrayT[p].row); 
       q=0;
       while (q < arrayT[p].colCount) {
    	  printf( "col%d: %d\n",q, arrayT[p].col[q].colnum);
	      printf( "Value%d: %.2f\n",q, *arrayT[p].col[q].kval);     
    	  q++;
       }
       p++; 
       printf("\n");   
    }
    printf("\n");*/
    printf("\n"); 
    p=0;  
    while (p < size) {
       q=0;
       while (q < arrayT[p].colCount) {
          printf( "%d  %d  %.2f\n",p,q, *arrayT[p].col[q].kval); 
    	  //printf( "col%d: %d\n",q, arrayT[p].col[q].colnum);
	      //printf( "Value%d: %.2f\n",q, *arrayT[p].col[q].kval);     
    	  q++;
       }
       p++; 
       printf("\n");   
    }
    printf("\n");
} 

            


  int main(int argc, char *argv[]){
      
     int *bits = calloc(TFBS, sizeof(int));
     int *viableStates;
     struct Ttype *arrayT;
    
     viableStates = malloc(30*sizeof(int));
     arrayT = malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));
     int array = 0;
    
     configure(0,bits,&array,viableStates);
     system("PAUSE");
     //printf("viableStates=%d\n", viableStates[14]);
     //printf("array=%d\n", array);
     
     transitions(array,viableStates,TFBS,arrayT,kon,koff,hammDist);
     print_arrayT(arrayT,array);
     
     
     int d;
     for (d=0; d<array; d++) {
       free(arrayT[d].col);
  }    
  free(arrayT);
  free(viableStates);
  free(bits);
     system("PAUSE");}
  