
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TFBS 18
#define SIZE 1
#define HIND_LENGTH 15

static int startPos[TFBS]= {0,0,6,6,7,7,7,7,8,8,12,12,14,14,17,17,20,20};
//static float kon[TFBS] = {1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
static float koff[5] = {1.5, 2.5, 3.5, 4.5, 5.5};
static int hammDist[TFBS] = {2,2,2,2,2,2,1,2,1,2,2,2,1,2,2,1,2,1};
//static float diag[TFBS]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

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
                           
int convertToDecimal(int *bits, int TFBSites){
     int n = TFBSites-1;
     int count =0;
     int record[TFBSites];
     int rec=0;
     while( n>=0 ){
         if(bits[n]==0){
            count++;
         }else{
             record[rec]=(TFBSites-1)-n;
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
 
                          
int isHindered(int bindSite, int *bits){
    int s = bindSite-1;
    int count=0;
    //printf("%d\n", s);
    //printf("%d\n", startPos[bindSite]);
    int check = startPos[bindSite]-HIND_LENGTH;
    if(check<0){
        check=0;}
   // printf("%d\n", check);
    while(s>=0 && startPos[s] <= startPos[bindSite] && startPos[s] >= check){
       if(bits[s] == 1){
                  //printf("hindered");
          return 1;
       } else {
          s--;
          count++;
       }
    }
    //printf("%d\n",count);
    return 0;
  }
  
  void configure(int bindSite, int *bits, int *numStates, int *statesArray, int TFBSites){

     if(bindSite<TFBSites-1){
        bits[bindSite] = 0;
        configure(bindSite+1,bits,numStates,statesArray,TFBSites);
        if(!isHindered(bindSite, bits)){
           bits[bindSite] = 1;
           configure(bindSite+1,bits,numStates,statesArray,TFBSites);
        }
     } else {
        bits[TFBSites-1] = 0;
        statesArray[(*numStates)] = convertToDecimal(bits, TFBSites);
        (*numStates)++;
        int i;
        for(i=0; i<TFBSites; i++){
           printf("%d", bits[i]);
        }
        printf("\n");
        if(!isHindered(bindSite, bits)){
           bits[TFBSites-1]=1;
           convertToDecimal(bits, TFBSites);
           statesArray[(*numStates)] = convertToDecimal(bits, TFBSites);
           (*numStates)++;
           for(i=0; i<TFBSites; i++){
              printf("%d", bits[i]);
           }
           printf("\n");
        }
      } 
     // system("PAUSE");
  }
  
  void diagonal(int row, float *diag, struct Ttype *arrayT, int m, int n){
    int x;
    diag[row]=0;
          for (x=0; x<m; x++) {
             diag[row] -= *(arrayT[n].col[x].kval);
          }
      arrayT[n].col[m].colnum = row;
      arrayT[n].col[m].kval = &(diag[row]);   
} 
  
  void transitions(int size, int *viableStates, int TFBSites, struct Ttype *arrayT, float kon[], float koff[5],int hammDist[TFBS], float *diag){
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
             printf( "%d  %d   %d  %.2f\n",p,q, arrayT[p].col[q].colnum,  *arrayT[p].col[q].kval); 
          //printf( "%d  %d  %.2f\n",p,q, *arrayT[p].col[q].kval); 
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
     float *diag; 
     diag = malloc(TFBS*sizeof(float));
     /*float *kon;
     kon = malloc(TFBS*sizeof(float));*/
     float kon[TFBS];
     int i;
     for(i=0; i<TFBS; i++){
              kon[i]=i+1;
     }
     
    //kon[TFBS]={1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
     viableStates = malloc(50*sizeof(int));
     arrayT = malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));
     int array = 0;
    
     configure(0,bits,&array,viableStates,TFBS);
     system("PAUSE");
     //printf("viableStates=%d\n", viableStates[14]);
     //printf("array=%d\n", array);
     
     transitions(array,viableStates,TFBS,arrayT,kon,koff,hammDist, diag);
     print_arrayT(arrayT,array);
     
     int test[12]={1,0,1,0,0,0,0,0,0,0,0,0};
     printf("%d", isHindered(1,test)); 
     
     
     int d;
     for (d=0; d<array; d++) {
       free(arrayT[d].col);
  }    
  free(arrayT);
  free(viableStates);
  free(bits);
  free(diag);
  //free(kon);
     system("PAUSE");}
  
