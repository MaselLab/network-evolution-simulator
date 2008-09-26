#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TFBS 12
#define SIZE 1
#define HIND_LENGTH 15

static int startPos[TFBS]= {1,2,3,3,10,10,11,11,11,13,16,16};

/*static int startPos[TFBS]= {1,2,3,3,10,10,11,11,11,13,16,16,17,17, 20,20,21,24,24,25,25,25};

/*static int startPos[TFBS]= {1,2,3,3,10,10,11,11,11,13,16,16,17,17,20,20,21,24,24,25,25,25,
                           26,26,26,27,28,28,29,32,36,37,37,41,44,46,48,51,53,53,57,57,57,
                           59,59,61,62,62,63,66,67,69,74,78,79,79,80,80,81,82,84,84,85,86,86,
                           88,88,89,90,90,91,91,92,92,93,94,100,100,100,103,103,104,104,110,
                           116,118,119,119,123,123,125,126,127};*/
 void convertToDecimal(int bits[TFBS]){
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
     printf("%d\n", decimal);
}  

void rowTransitions(int row){
         
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
  
  void configure(int bindSite, int bits[TFBS]){

     if(bindSite<TFBS-1){
        bits[bindSite] = 0;
        configure(bindSite+1,bits);
        if(!isHindered(bindSite, bits)){
           bits[bindSite] = 1;
           configure(bindSite+1,bits);
        }
     } else {
        bits[TFBS-1] = 0;
        convertToDecimal(bits);
        int i;
        for(i=0; i<TFBS; i++){
           printf("%d", bits[i]);
        }
        printf("\n");
        if(!isHindered(bindSite, bits)){
           bits[TFBS-1]=1;
           convertToDecimal(bits);
           for(i=0; i<TFBS; i++){
              printf("%d", bits[i]);
           }
           printf("\n");
        }
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


            


  int main(int argc, char *argv[]){
      
     int bits[TFBS] = {0,0,0,0,0,0,0,0,0,0,0,0};
     int *viableStates;
     viableStates = malloc(20*sizeof(int));
     //printf("boolean: %d\n",isHindered(5, bits));
     //int numStates =0;
     configure(0,bits);
     //printf("numStates = %d\n", numStates);
     //printf("hello world\n");
     int testBits[TFBS] = {1,0,0,0,0,0,0,0,0,0,0,1};
    // convertToDecimal(testBits);
     
     unsigned int test = 2;
     if(test & (0)){
       printf("test=%d\n", (test&2049));
     }
     
     system("PAUSE");}
  
