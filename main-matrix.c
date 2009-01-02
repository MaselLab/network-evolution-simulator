#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <errno.h>

#include "random.h"
#include "lib.h"
#include "netsim.h"



int array_size = 1000;

FILE *sparseMatrixV1;
FILE *statesV1;
FILE *columnV1;
FILE *leftEdgePositions;
FILE *bVector;

typedef int (*compfn)(const void*, const void*);

int intcmp(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

int compare(struct AllTFBindingSites  *elem1, struct AllTFBindingSites *elem2){
    if(elem1->leftEdgePos < elem2->leftEdgePos)
        return -1;
    else if(elem1->leftEdgePos > elem2->leftEdgePos)
         return 1;
    else
         return 0;
}

 struct Rowtype {
    long rownum;
    float *kval;
};
    
struct Ttype {
  int col;
  struct Rowtype *row;
  int rowCount;
};  

unsigned long convertToDecimal(int *bits, int TFBS){
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
     unsigned long decimal=0;
     for(i=0; i<rec;i++){
        decimal+= (long)pow(2,record[i]);
     }
     return decimal;
}              

int isHindered(int bindSite, int *bits, int *startPos){
    int s = bindSite-1;
    int count=0;
    //printf("%d\n", s);
    //printf("%d\n", startPos[bindSite]);
    int check = startPos[bindSite]-(HIND_LENGTH);
    if(check<0){
        check=0;}
    //printf("%d\n", check);
    while(s>=0 && startPos[s] <= startPos[bindSite] && startPos[s] >= check){
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
  
void configure(int bindSite, int *bits, int *numStates, unsigned long *statesArray, int TFBS, int *startPos){

     if(bindSite<TFBS-1){
        bits[bindSite] = 0;
        configure(bindSite+1,bits,numStates,statesArray,TFBS,startPos);
        if(!isHindered(bindSite, bits,startPos)){
           bits[bindSite] = 1;
           configure(bindSite+1,bits,numStates,statesArray,TFBS,startPos);
        }
     } else {
        bits[TFBS-1] = 0;
        statesArray[(*numStates)] = convertToDecimal(bits, TFBS);
        (*numStates)++;
       int i;
        for(i=0; i<TFBS; i++){
          printf("%d", bits[i]);
          }
        printf("\n");
        if(!isHindered(bindSite, bits,startPos)){
           bits[TFBS-1]=1;
           //convertToDecimal(bits, TFBS);
           //printf("%d\n", TFBS);
           statesArray[(*numStates)] = convertToDecimal(bits, TFBS);
           //printf("numstates=%d, %lu convert=%lu\n", (*numStates), statesArray[(*numStates)], convertToDecimal(bits, TFBS));
           (*numStates)++;
           for(i=0; i<TFBS; i++){
              printf("%d", bits[i]);
           }
           printf("\n");
        }
      } 
          //system("PAUSE");
  }

void diagonal(int col, float *diag, struct Ttype *arrayT, int m, int n){
      
    int x;
    float value=0;
    diag[col]=0;
          for (x=0; x<m; x++) {
             value -= *(arrayT[n].row[x].kval);
          }
      //printf("%d   %.2f  %.2f\n", col, diag[col], value);
      diag[col]=value;
      arrayT[n].row[m].rownum = col;
      arrayT[n].row[m].kval = &(diag[col]); 
      
} 
  
void transitions(int startSite,int size, unsigned long *viableStates, int TFBSites, struct Ttype *arrayT, float kon[], float koff[5],int *hammDist, float *diag, int *TFon){
       statesV1 = fopen("statesV1.txt", "w");
       columnV1 = fopen("columnV1.txt","w");
     if ((statesV1 = fopen("statesV1.txt", "w"))) {
     if ((columnV1 = fopen("columnV1.txt", "w"))) {
      
       int i, p, j, m, a;
       unsigned long row;
       int n=0;
       printf("size=%d TFBS= %d\n", size, TFBSites);
       for( i=0;i<size; i++){
          arrayT[n].col = i;
          arrayT[n].row = malloc((array_size)*sizeof(struct Rowtype));
         printf("viableStates:%lu, col num:%d\n",viableStates[i], i);
       fprintf(statesV1,"%lu \n", viableStates[i]);
         fprintf(columnV1,"%d\n",i);
        m=0;
        row = viableStates[i];
        for(p=0;p<TFBSites;p++){
     
          if(!(row & (1L<<p))){     
              //need to do decimal to binray here!!!!!!!!  
             row = viableStates[i] | (1L<<p);
             if(row!=0 && row!=viableStates[i]){
                for(j=0;j<size; j++){
                   if(row==viableStates[j]){
                      arrayT[n].row[m].rownum = j;
                      a = hammDist[p];
                      arrayT[n].row[m].kval = &(koff[a]);
                      m++;
                     // printf("rownum=%d, j=%d\n", arrayT[n].row[m].rownum, j);
                      //printf("koff=%f\n", koff[a]);
                     //printf("    col=%d, rownum=%d, j=%d, p=%d\n",row, &(arrayT[n].row[m].rownum),j, p);
                     //printf(" %d  %d  \n", i, j);
                   }
                 }
              }
          }else{
              row = viableStates[i] ^ (1L<<p);
              if(row != viableStates[i]){     
                 for( j=0;j<size; j++){   
                    if(row==viableStates[j]){
                    
                      arrayT[n].row[m].rownum = j;
                      arrayT[n].row[m].kval = &(kon[p]);
                      m++;
                      //printf("rownum=%d, j=%d\n", arrayT[n].row[m].rownum, j);
                      //printf("    col=%d, rownum=%d, j=%d, p=%d\n", row, &(arrayT[n].row[m].rownum),j, p);
                     //printf(" %d  %d  \n", i, j);
                    }
                  }
              }//printf("\n");
          }
       }
       diagonal(i,diag, arrayT, m, n);
       m++;
       arrayT[n].rowCount = m;
      // printf("rowCount=%d\n",  arrayT[n].rowCount);
       n++;
       }
     }}
     fclose(statesV1);
  }
  
void print_arrayT(struct Ttype *arrayT, int size, unsigned long *viableStates){
     int p, q;  
     printf("%d\n",size); 
     //printf("CHECK:%d  %d | %d  %d  %d\n", 64,arrayT[64].row[0].rownum, viableStates[64], viableStates[arrayT[64].row[0].rownum],  *arrayT[64].row[0].kval ); 
   //  system("PAUSE");
    p=0;  
    while (p < size) {
       q=0;
       //printf("rowCount=%d, q=%d\n",  arrayT[p].rowCount,q);
       while (q < arrayT[p].rowCount) {
          //if(arrayT[p].row[q].rownum!=4){
         //printf("rownum=%d\n", &(arrayT[p].row[q].rownum));
          printf( "%d  %lu | %lu  %lu  %.2f\n",p,(arrayT[p].row[q].rownum), viableStates[p], viableStates[arrayT[p].row[q].rownum],*arrayT[p].row[q].kval); 
    	 // }else{
              //printf( "%d  %d | %d  %d  %d\n",p,arrayT[p].row[q].rownum, viableStates[p], viableStates[arrayT[p].row[q].rownum],  1);   
          //}//printf( "col%d: %d\n",q, arrayT[p].col[q].colnum);
	     // printf( "Value%d: %.2f\n",q, *arrayT[p].row[q].kval);     
    	  q++;
       }
       p++; 
       //printf("%d\n", size);  
      // system("PAUSE"); 
    }
    //system("PAUSE"):
    printf("\n");
} 

void print_arrayT_MATLAB(struct Ttype *arrayT, int size, unsigned long *viableStates){
     int p, q;  
     printf("\n Col   Row      kval\n"); 
    p=0;  
    while (p < size) {
       q=0;
       while (q < arrayT[p].rowCount) {
           if(arrayT[p].row[q].rownum!=4){
           printf(" %d      %lu      %f\n", p,arrayT[p].row[q].rownum,  *arrayT[p].row[q].kval); 
           }else{
                printf(" %d      %lu      %d\n", p,arrayT[p].row[q].rownum,  1);  
           }
           q++;
       }
       p++;
    }
    int count=0;
       while(count<size){
          printf(" %d     %d    %d\n", count, 4, 1);
          count++;
          }
    //system("PAUSE");
   
    sparseMatrixV1 = fopen("sparseMatrixV1.txt", "w");
    if ((sparseMatrixV1 = fopen("sparseMatrixV1.txt", "w"))) {
        //ask alex about this line of code...
       //fprintf(sparseMatrixV1, "1Row  Column  Value\n\n");
        //printf( "1Row  Column  Value\n\n");
        int p=0;  
      while (p < size) {
      int q=0;
       while (q < arrayT[p].rowCount) {
       if(arrayT[p].row[q].rownum!=4){
       fprintf(sparseMatrixV1, "%ld,   %d,   %f\n" ,arrayT[p].row[q].rownum +1,p+1,  *arrayT[p].row[q].kval);
      }else{} 
      // printf( "%d,   %d,   %.2f\n\n", arrayT[6].row, arrayT[6].col[3].colnum, *arrayT[6].col[3].kval);
       q++;
       }
       p++;
      }
     int count=0;
       while(count<size){
          fprintf(sparseMatrixV1," %d     %d    %d\n", 5, count+1, 1);
          count++;
          }
      }
    fclose(sparseMatrixV1);
}

void print_vector_MATLAB(int size){
     int i;
      bVector = fopen("bVector.txt", "w");
      if ((bVector = fopen("bVector.txt", "w"))) {
       for(i=0; i<size; i++){
         if(i!=4){
            fprintf(bVector, "%d\n", 0);
         }else{
            fprintf(bVector, "%d\n", 1);
         }
       }
       fclose(bVector);
     }
}

int TestByteOrder()
{
   short int word = 0x0001;
   char *byte = (char *) &word;
   return byte[0];
   //return(byte[0] ? LITTLE_ENDIAN : BIG_ENDIAN);
}
  
  
void convertFile1( char *fileName, unsigned long *vect, int size){
     int n;
     FILE *file;
     file = fopen(fileName, "r");
     if(file==NULL){
        printf("Unable to open file!");
     } else {
         for(n=0; n<size; n++){
            int success = fscanf(file, "  %lu", &(vect[n]));
            if (!success) {
               printf("The file ran out of input. Check the size parameter to convertFile.\n");              
            }
         }
         
         fclose(file);
     }
}

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

void probSlide(unsigned long *statesArray, float *prob, float *outcome, int size, float *previous){
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
     system("PAUSE");
     for(y=0; y<size; y++){
         printf("%lu\n",statesArray[y]);    
     //}    
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

int nextStartSite(int startSite, int *startPos){
    int next;
    int basePairNum;
    basePairNum = startPos[startSite];
    printf("startSite = %d\n basePairNum = %d\n",startSite, basePairNum);
    next = basePairNum;
    while(next==basePairNum){
       next++;
    }
    return next;
}

void populateFinal(float zero, float ones, float twos, float *final, int nextSite){
     final[0] = zero;
     final[nextSite] = ones;
     final[nextSite+1] = twos;
}

int main(int argc, char *argv[])
{
    /*NOTE: To make following code compatible with all platfroms and structures need to implement doe
    for 4 types: 
        - 32bit Little-Endian
        - 32bit Big-Endian
        - 64bit Little-Endian
        - 64bit Big-Endian
    Current code wirtten for 64bit Big-Endian
    */
    
    //Check Endianess of system
    int testEndian = TestByteOrder();
    if(testEndian==1){printf("%d = LittleEndian\n\n", testEndian);}
    else{printf("%d = BigEndian--WARNING: This code will not work correctly on BigEndian systems!\n\n", testEndian);}
    
    //Check platform of system
    int testPlatform = sizeof(long);
    if(testPlatform == 4){printf("32-bit : Number of Binding Sites cannot exceed 32sites\n\n");}
    else{printf("64-bit\n\n");}
    
  FILE *fpkdis;
  char fperrors_name[80];
  int i, j;
  Genotype indiv;
  float initProteinConc[NGENES], kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int curr_seed;
  int TFBS;
  
  TFBS = 10;
  verbose = 0;

  /* change to get a different genotype */
  dummyrun = 4;
  current_ploidy = 1;

  /* parse command-line options */
  while ((c = getopt (argc, argv, "hvd:r:p:w:a:")) != -1) {
    switch (c)
      {
      case 'd':
        output_directory = optarg;
        break;
      case 'r':
        dummyrun = atoi(optarg);
        break;
      case 'p':
        current_ploidy = atoi(optarg);
        break;
      case 'w':
        TFBS = atoi(optarg);
        break;
      case 'a':
        array_size = atoi(optarg);
        break;        
     case 'v':
        verbose = 1;
        break;
      case 'h':
        fprintf(stderr, "%s [-d DIRECTORY] [-r DUMMYRUN] [-h] [-p PLOIDY] [-w WINDOWSIZE] [-a ARRAYSIZE]\n", argv[0]);
        exit(0);
        break;
      default:
        abort();
      }
  }

/* create output directory if needed */
#ifdef __unix__
  directory_success = mkdir(output_directory, S_IRUSR|S_IWUSR|S_IXUSR);
#else 
#ifdef __WIN32__
  directory_success = mkdir(output_directory);
#endif
#endif

  if (directory_success==-1) {
    if (errno == EEXIST) {
      fprintf(stderr, "directory '%s' already exists\n", output_directory);
    } else {
      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
      exit(-1);
    }
  }

  sprintf(fperrors_name, "%s/netsimerrors.txt", output_directory);
  fperrors = fopen(fperrors_name, "w");

  /* get the kdis.txt values */
  fpkdis = fopen("kdis.txt","r");
  for (j = 0; j < NUM_K_DISASSEMBLY; j++) {
    fscanf(fpkdis,"%f", &kdis[j]);
  }
  fclose(fpkdis); 

  /* change random number */
  for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);


  /**********************************************************************
   *  modifiable part starts here
   **********************************************************************/

  /* initialize protein and mRNA concentrations */
  /* Jasmin: you can use these initial concentrations for computing your kons */
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    //printf("%f\n", initProteinConc[i]);
  }
 // printf("\n");
 
  
    
    //populate Koff
    float Koff[5];
    float Gibbs;
    float RTlnKr;
    float koffCheck;
    float temperature = 293.0;
    RTlnKr = GasConstant * temperature * log(Kr);
    int jo;
    for( jo=0; jo<3; jo++){
         Gibbs = (((float) jo)/3.0 - 1.0) * RTlnKr;
         //printf("\n gibbs = %f\n", Gibbs);
         koffCheck = NumSitesInGenome*kon*0.25/exp(-Gibbs/(GasConstant*temperature));
         //printf("\n koffCheck = %f\n", koffCheck);
         Koff[jo] = koffCheck;
    }
    int jojo;
    for(jojo=0; jojo<3; jojo++){
       printf("Koff[%d] = %f\n", jojo, Koff[jojo]);
    }
    //*koff = NumSitesInGenome*kon*0.25/exp(-Gibbs/(GasConstant*state->temperature));
    //Gibbs = (((float) allBindingSites[k].hammingDist)/3.0 - 1.0) * state->RTlnKr;
    //Koff = [0 mismatch, 1 mismatch, 2 mismatch, coop on 1 side, coop on 2 sides]
    

  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, kdis);
  
  /* print binding sites */
 /* print_all_binding_sites(indiv.copies, indiv.allBindingSites, indiv.bindSiteCount, 
			  indiv.transcriptionFactorSeq, indiv.cisRegSeq, indiv.tfsStart); 
			  printf("tfsPerGene = %d", indiv.tfsPerGene); */


  /* Jasmin:
     pseudo-code for looping (see also print_all_binding_sites in netsim.c):*/
     //system("PAUSE");
    //int sitePos[10];
    //int transFactor[10];

    int *leftEdgePos;
    int *startPos;
    int *hammDist;
    float *diag;
    int *TFon;
    float *final;
    final = malloc(array_size*sizeof(float));
    
      //populate Kon
    float Kon[TFBS];
    
    leftEdgePos = malloc(indiv.tfsPerGene[0]*sizeof(int));
    startPos=malloc(TFBS*sizeof(int));
    //startPos=malloc(indiv.tfsPerGene[0]*sizeof(int));
    hammDist = malloc(TFBS *sizeof(int));
    diag = malloc(array_size*sizeof(float));
    TFon = malloc(TFBS*sizeof(int));
    //Kon = malloc(TFBS*sizeof(float));
      
    int *bits = calloc(TFBS, sizeof(int));
    unsigned long *viableStates;
    struct Ttype *arrayT;
    
    viableStates = malloc((array_size)*sizeof(unsigned long));
    arrayT = malloc(array_size*sizeof(struct Ttype));
    // arrayT = malloc((pow(2,TFBS))*sizeof(struct Ttype));
     int array = 0; 
     
     int finalPos =0;
     
    float *previous;
   previous = malloc(10*sizeof(float));
   float prev[10] = {1,0,0,0,0,0,0,0,0,0};
    int pi;
    for(pi=0; pi<10; pi++){
       previous[pi]=prev[pi];
    }
    
    int fey;
    for( fey=0; fey<TFBS; fey++){
         printf("binding site %3d:  ", fey);
         printf("%d\n", indiv.allBindingSites[fey].leftEdgePos);
         }
         
    leftEdgePositions = fopen("leftEdgePositions.txt", "w");
     if ((leftEdgePositions = fopen("leftEdgePositions.txt", "w"))) {
         fprintf(leftEdgePositions , "\n");
         for (i=0; i <indiv.tfsPerGene[0] ; i++) {  
             fprintf(leftEdgePositions, "binding site %3d:  ", i);
             fprintf(leftEdgePositions, "%d\n", indiv.allBindingSites[i].leftEdgePos);
         }
         fprintf(leftEdgePositions, "\n");
    
         qsort((void *) &(indiv.allBindingSites[0]), indiv.tfsPerGene[0],                                 
            sizeof(struct AllTFBindingSites),(compfn)compare );
         printf("\n");
  
         for( fey=0; fey<TFBS; fey++){
             printf("binding site %3d:  ", fey);
             printf("%d\n", indiv.allBindingSites[fey].leftEdgePos);
         }
         printf("tfsPerGene = %d", indiv.tfsPerGene[0]);
         system("PAUSE");
               //leftEdgePositions = fopen("leftEdgePositions.txt", "w");
                // if ((leftEdgePositions = fopen("leftEdgePositions.txt", "w"))) {
         for (i=0; i <indiv.tfsPerGene[0] ; i++) {
              /*if(indiv.allBindingSites[i].cisregID ==1){
                  printf("One:%d  LeftEdge:%d\n", i, indiv.allBindingSites[i].leftEdgePos);
               }*/
    
              fprintf(leftEdgePositions, "binding site %3d:  ", i);
              leftEdgePos[i] =  indiv.allBindingSites[i].leftEdgePos;
            /*fprintf( leftEdgePositions, "binding site %3d:  ", i);
            fprintf(leftEdgePositions, "       cis-reg region: %3d",indiv.allBindingSites[i].cisregID);
            fprintf(leftEdgePositions, "         cis-reg copy: %3d", indiv.allBindingSites[i].geneCopy);
            fprintf(leftEdgePositions, " (sequence %.*s)\n", CISREG_LEN, indiv.cisRegSeq[indiv.allBindingSites[i].cisregID][indiv.allBindingSites[i].geneCopy]);
            fprintf(leftEdgePositions, " transcription-factor: %3d", indiv.allBindingSites[i].tfID);
            fprintf(leftEdgePositions, " (sequence: %.*s)\n", TF_ELEMENT_LEN, indiv.transcriptionFactorSeq[indiv.allBindingSites[i].tfID][indiv.allBindingSites[i].geneCopy]); 
            fprintf(leftEdgePositions, "  L-edge of %2dbp hind: %3d\n", HIND_LENGTH, indiv.allBindingSites[i].leftEdgePos);        
            
            */
              fprintf(leftEdgePositions, "%d\n", indiv.allBindingSites[i].leftEdgePos);
    
            /*//fprintf(leftEdgePositions,  "%d\n", indiv.allBindingSites[i].leftEdgePos);
            fprintf(leftEdgePositions, "  Hind offset position: %3d\n", indiv.allBindingSites[i].hindPos); 
            fprintf(leftEdgePositions, "               strand: %3d\n", indiv.allBindingSites[i].strand);
            fprintf(leftEdgePositions, "         Hamming dist: %3d\n\n", indiv.allBindingSites[i].hammingDist); 
           */}
    }
    fclose(leftEdgePositions);
    printf("leftEdgePos[1] = %d", leftEdgePos[1]);
    system("PAUSE");
   
     printf("\n");
     int startSite;
     startSite=0;
     //while loop for sliding window
     while(startSite<indiv.tfsPerGene[0]){
     int lem;
     
     int bob;
      
     for(lem =0; lem<TFBS; lem++){
             startPos[lem] = indiv.allBindingSites[lem+startSite].leftEdgePos;
             hammDist[lem] = indiv.allBindingSites[lem+startSite].hammingDist;
             TFon[lem] = indiv.allBindingSites[lem+startSite].tfID;
             bob = indiv.allBindingSites[lem+startSite].tfID;
             //printf("bob= %d\n", bob);
             Kon[lem] = initProteinConc[bob]*kon;
             
             printf("%d", indiv.allBindingSites[lem+startSite].leftEdgePos);
             printf(" Hd = %d   tf = %d", hammDist[lem], TFon[lem]);
             printf(" Kon[lem] = %f\n", Kon[lem]);    
        }
     system("PAUSE");
     
    /* printf("startPos[0] = %d\n", startPos[0]);
      int hello;
      hello = nextStartSite(0,startPos);
      printf("hello = %d\n",hello);
      system("PAUSE");*/
    
     printf("\n");
     configure(0,bits,&array,viableStates,TFBS,startPos);

    int sh;
    for(sh=0; sh<array; sh++){
              printf("%lu\n", viableStates[sh]);
     }
 
     int y;
     int ones, twos, none;
     unsigned long *onesArray, *twosArray, *noneArray;
     noneArray=malloc(500*sizeof(unsigned long));
     onesArray= malloc(200*sizeof(unsigned long));
     twosArray = malloc(200*sizeof(unsigned long));
     ones =0;
     twos=0;
     none =0;
     for(y=0; y<array; y++){
         printf("%lu\n",viableStates[y]);    
           
        if(((viableStates[y])>>1)%2==1 ){
           twosArray[twos] = viableStates[y];
           twos++;
        }else if((viableStates[y]-1)%2==0){
              onesArray[ones] = viableStates[y];
              ones++;
        } else{
               noneArray[none] = viableStates[y];
               none++;
        }
     }        
      printf("\n none=%d, ones=%d, twos=%d\n\nonesArray\n",none, ones, twos); 
      
      for(y=0; y<ones; y++){
         printf("%lu\n", onesArray[y]);
         }
         printf("twosArray\n");
      for(y=0; y<twos; y++){
         printf("%lu\n", twosArray[y]);
         }
         printf("noneArray\n");
      for(y=0; y<none; y++){
         printf("%lu\n", noneArray[y]);
         } 
 
     printf("%d\n", array);
     
     print_vector_MATLAB(array);
     
     system("PAUSE");
     
     transitions(startSite,array,viableStates,TFBS,arrayT, Kon, Koff, hammDist, diag, TFon);
  
     system("PAUSE");

     print_arrayT_MATLAB(arrayT,array,viableStates);
     
     if(system("matlab -nodisplay -nojvm -nodesktop -nosplash -r \"pVect; exit;\"")!=0){
                       printf("Problem running Matlab.");
                       exit(1);   
     }else{
        printf("\n");
     }
     
   unsigned long *statesS;
   statesS = malloc(array*sizeof(unsigned long));
   float *outcome;
/*   float *previous;
   previous = malloc(10*sizeof(float));
   float prev[10] = {1,0,0,0,0,0,0,0,0,0};
    int pi;
    for(pi=0; pi<10; pi++){
       previous[pi]=prev[pi];
    }*/
    
    system("PAUSE");
   
    outcome = malloc( array*sizeof(float));
    
    float *vector;
    vector = malloc(array*sizeof(float));
    
    convertFile1("statesV1.txt", statesS, array);
    
    convertFile("b.txt", vector, array);
    
    int n;
    printf("\nvector\n");
    for(n=0; n<array; n++){
            printf( "%f\n", (vector[n]));
    }
    
    probSlide(statesS, vector, outcome, array, previous);
   populateFinal(outcome[0], outcome[1], outcome[2], final,pow(2,startSite));
   printf("FINAL\n");
   for(i=0; i<(pow(2,startSite)+2); i++){
       printf("%.4f\n", final[i]);
    }
    system("PAUSE");
     int nextS = nextStartSite(startSite, leftEdgePos);
     startSite = nextS;
     printf("startSiteHERE = %d\n", startSite);
     system("PAUSE");
   }
     //printf("long=%ud\n", sizeof(long));
    // printf("hindlength=%d\n", HIND_LENGTH);
      //printf("tfsPerGene = %d", indiv.tfsPerGene[0]);
     //  system("PAUSE");
  /* free dynamically allocated all binding sites list */
  free(indiv.allBindingSites);
   int d;
     for (d=0; d<array; d++) {
       free(arrayT[d].row);
     }   
  
  free(arrayT);
  free(startPos);
  free(viableStates);
  free(bits);
  free(hammDist);
  free(TFon);
  free(diag);
  /* close error file */
  fclose(fperrors);
  //fclose(sparseMatrixV1);
}
