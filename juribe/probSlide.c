#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM 15
 
static unsigned long statesA[NUM] = {0,1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192};
static float probs[NUM] = {.0009,.0417,.0417,.0138,.0138,.1184,.1184,.0417,.0417,.1885,.1885,.0064,.0064,.0889,.0889};
static int prev[8] = {1,0,0,0,0,0,0,0};

void convertFile( char *fileName, float *vect, int size){
     int n;
     FILE *file;
     file = fopen(fileName, "r");
     if(file==NULL){
        printf("Unable to open file!");
     } else {
         for(n=0; n<size; n++){
            int success = fscanf(file, "  %f", &(vect[n]));
            if (!success) {
               printf("The file ran out of input. Check the size parameter to convertFile.\n");              
            }
         }
         
         fclose(file);
     }
}

void probSlide(int *statesArray, float *prob, float *outcome, int size, float *previous){
     //separate into 00, 01, 10
     int y;
     int ones, twos, none;
     float *onesArray, *twosArray, *noneArray;
     noneArray=malloc(50*sizeof(float));
     onesArray= malloc(20*sizeof(float));
     twosArray = malloc(20*sizeof(float));
     ones =0;
     twos=0;
     none =0;
     for(y=0; y<size; y++){
         printf("%d\n",statesArray[y]);    
           
        if(((statesArray[y])>>1)%2==1 ){
           twosArray[twos] = prob[y];
           twos++;
        }else if((statesArray[y]-1)%2==0){
              onesArray[ones] = prob[y];
              ones++;
        } else{
               noneArray[none] = prob[y];
               none++;
        }
     }        
      printf("\n none=%d, ones=%d, twos=%d\n\nonesArray\n",none, ones, twos); 
      
      for(y=0; y<ones; y++){
         printf("%.4f\n", onesArray[y]);
         }
         printf("twosArray\n");
      for(y=0; y<twos; y++){
         printf("%.4f\n", twosArray[y]);
         }
         printf("noneArray\n");
      for(y=0; y<none; y++){
         printf("%.4f\n", noneArray[y]);
         }  
      //Add up probabilities
      
      float sum1, sum2, sum0;
      sum1 = 0;
      sum2 = 0;
      sum0 = 0;
      for(y=0; y<ones; y++){
         sum1 += onesArray[y];
         }
      for(y=0; y<twos; y++){
         sum2 += twosArray[y];
         }
      for(y=0; y<none; y++){
         sum0 += noneArray[y];
         }
         //normalize to previous
      outcome[0] = sum0*previous[0];
      outcome[1] = sum1*previous[0];
      outcome[2] = sum2*previous[0];
      printf("\n%.4f %.4f %.4f\n", sum0, sum1, sum2); 
       
     int g;   
     for(g=0; g<3; g++){
        previous[g] = outcome[g];
        }
        
      
         
}

int main(int argc, char *argv[]){
   int *states;
   states = malloc(20*sizeof(unsigned long));
   float *prob;
   float *outcome;
   float *previous;
   previous = malloc(10*sizeof(float));
   
   //previous[0]=1;
    
    int i;
    for(i=0; i<10; i++){
       previous[i]=prev[i];
    }
    for(i=0; i<NUM; i++){
       states[i]=statesA[i];
    }
    
    for(i=0; i<NUM; i++){
      printf("%d\n", states[i]);
    }
    system("PAUSE");
    prob = malloc(20*sizeof(float));
    for(i=0; i<NUM; i++){
       prob[i]=probs[i];
    }
    system("PAUSE");
     for(i=0; i<NUM; i++){
       printf("%.4f\n",prob[i]);
    }
    outcome = malloc(20*sizeof(float));
    
    float *vector;
    vector = malloc(16*sizeof(float));
    system("PAUSE");
    
    convertFile("b.txt", vector, 16);
    int n;
    
    for(n=0; n<16; n++){
            printf( "%f\n", (vector[n]));
    }
    probSlide(states, prob, outcome, NUM, previous);
    
    /*printf("\n");
    int k;
    for (k=0; k<3; k++){
        printf("%.4f\n", outcome[k]);
    }
    printf("\n");
    for (k=0; k<3; k++){
        printf("%.4f\n", previous[k]);
    }*/
    system("PAUSE");
}
    
    
