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

FILE *sparseMatrixV1;

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

 struct Coltype {
    int colnum;
    float *kval;
};
    
struct Ttype {
  int row;
  struct Coltype *col;
  int colCount;
};  

int convertToDecimal(int *bits, int TFBS){
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

int isHindered(int bindSite, int *bits, int *startPos){
    int s = bindSite - 1;
    int count=0;
    //printf("%d\n", s);
    //printf("%d\n", startPos[bindSite]);
    int check = startPos[bindSite]-HIND_LENGTH;
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
  
void configure(int bindSite, int *bits, int *numStates, int *statesArray, int TFBS, int *startPos){

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
          // printf("%d", bits[i]);
       
        }
        //printf("\n");
        if(!isHindered(bindSite, bits,startPos)){
           bits[TFBS-1]=1;
           convertToDecimal(bits, TFBS);
           statesArray[(*numStates)] = convertToDecimal(bits, TFBS);
           (*numStates)++;
           for(i=0; i<TFBS; i++){
             // printf("%d", bits[i]);
           }
          // printf("\n");
        }
      } 
          //system("PAUSE");
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
  
 void transitions(int size, int *viableStates, int TFBSites, struct Ttype *arrayT, float kon[], float koff[5],int *hammDist, float *diag, int *TFon){
       int i, p,j,m, tf;
       int n=0;
       for( i=0;i<size; i++){
          arrayT[n].row = i;
          arrayT[n].col = malloc((TFBSites+2)*sizeof(struct Coltype));
         //printf("viableStates:%d, row num:%d\n",viableStates[i], i);
        m=0;
        for(p=0;p<TFBSites;p++){
          int col = viableStates[i];
           // system("PAUSE");
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
                 arrayT[n].col[m].kval = &kon[p];
                 m++;
                 //printf("    col=%d, j=%d, p=%d\n", col, j, p);
                // printf(" %d  %d  \n", i, j);
               }
             }
           }//printf("\n");
          }
       }
      
      // printf("kval00: %f   kon0: %f %m= %d  n= %d  i= %d\n", *arrayT[0].col[0].kval, kon[0], m, n, i);
       diagonal(i,diag, arrayT, m, n);
      
      // printf("CHECK HERE!!!!!!! kval00: %f   kon0: %f\n", *arrayT[0].col[0].kval, kon[0]);
       m++;
       arrayT[n].colCount = m;
       n++;
     }
  }
  
  void print_arrayT(struct Ttype *arrayT, int size, int *viableStates){
     int p, q;  
     printf("\n"); 
    p=0;  
    while (p < size) {
       q=0;
       while (q < arrayT[p].colCount) {
          printf( "%d  %d | %d  %d  %.2f\n",p,arrayT[p].col[q].colnum, viableStates[p], viableStates[arrayT[p].col[q].colnum],  *arrayT[p].col[q].kval); 
    	  //printf( "col%d: %d\n",q, arrayT[p].col[q].colnum);
	      //printf( "Value%d: %.2f\n",q, *arrayT[p].col[q].kval);     
    	  q++;
       }
       p++; 
       printf("\n");   
    }
    printf("\n");
} 

void print_arrayT_MATLAB(struct Ttype *arrayT, int size, int *viableStates){
      int p, q;  
     printf("\n Row    Col      kval\n"); 
    p=0;  
    while (p < size) {
       q=0;
       while (q < arrayT[p].colCount) {
           printf(" %d      %d      %f\n", p,arrayT[p].col[q].colnum,  *arrayT[p].col[q].kval); 
           q++;
       }
       p++;
    }
    
    sparseMatrixV1 = fopen("sparseMatrixV1.txt", "w");
     if (sparseMatrixV1 = fopen("sparseMatrixV1.txt", "w")){
        //ask alex about this line of code...
       //fprintf(sparseMatrixV1, "1Row  Column  Value\n\n");
        //printf( "1Row  Column  Value\n\n");
        p=0;  
      while (p < size) {
       q=0;
       while (q < arrayT[p].colCount) {
       fprintf(sparseMatrixV1, "%d,   %d,   %f\n", p+1,arrayT[p].col[q].colnum +1,  *arrayT[p].col[q].kval); 
       //printf( "%d,   %d,   %.2f\n\n", arrayT[6].row, arrayT[6].col[3].colnum, *arrayT[6].col[3].kval);
       q++;
       }
       p++;
      }
      }
    
}
  

int main(int argc, char *argv[])
{
  FILE *fpkdis;
  char fperrors_name[80];
  int i, j, k;
  Genotype indiv;
  float initProteinConc[NGENES], kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int hold_genotype_constant = 0;
  int curr_seed;

  verbose = 0;

  /* change to get a different genotype */
  dummyrun = 4;

  /* parse command-line options */
  while ((c = getopt (argc, argv, "hvgd:r:p:t:c:")) != -1) {
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
      case 't':
        tdevelopment = atof(optarg);
        break;
      case 'c':
        critical_size = atof(optarg);
        break;
      case 'g':
        hold_genotype_constant = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        fprintf(stderr, "%s [-d DIRECTORY] [-r DUMMYRUN] [-h] [-g] [-p PLOIDY] [-t DEVELOPMENTTIME] [-c CRITICALSIZE]\n", argv[0]);
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

  if (directory_success==-1) 
    if (errno == EEXIST) {
      fprintf(stderr, "directory '%s' already exists\n", output_directory);
    } else {
      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
      exit(-1);
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
  /*print_all_binding_sites(indiv.ploidy, indiv.allBindingSites, indiv.bindSiteCount, 
			  indiv.transcriptionFactorSeq, indiv.cisRegSeq); 
  printf("tfsPerGene = %d", indiv.tfsPerGene);

  /* Jasmin:
     pseudo-code for looping (see also print_all_binding_sites in netsim.c):*/
     //system("PAUSE");
    //int sitePos[10];
    //int transFactor[10];
    int TFBS;
    TFBS = 15;
    int *startPos;
    int *hammDist;
    float *diag;
    int *TFon;
    
      //populate Kon
    float Kon[TFBS];
    
    
    
    startPos=malloc(TFBS*sizeof(int));
    //startPos=malloc(indiv.tfsPerGene[0]*sizeof(int));
    hammDist = malloc(TFBS *sizeof(int));
    diag = malloc(TFBS*sizeof(float));
    TFon = malloc(TFBS*sizeof(int));
    //Kon = malloc(TFBS*sizeof(float));
      
    int *bits = calloc(TFBS, sizeof(int));
    int *viableStates;
    struct Ttype *arrayT;
    
    viableStates = malloc((100)*sizeof(int));
     arrayT = malloc((pow(2,TFBS)+1)*sizeof(struct Ttype));
     int array = 0;
    
    
    qsort((void *) &(indiv.allBindingSites[0]), indiv.tfsPerGene[0],                                 
           sizeof(struct AllTFBindingSites),(compfn)compare );
   printf("tfsPerGene = %d", indiv.tfsPerGene[0]);
   
   
  /* for (i=0; i <indiv.tfsPerGene[0] ; i++) {
    printf("binding site %3d:\n", i);
    printf("       cis-reg region: %3d",indiv.allBindingSites[i].cisregID);
    printf("         cis-reg copy: %3d", indiv.allBindingSites[i].geneCopy);
    //printf(" (sequence %.*s)\n", CISREG_LEN, cisRegSeq[allBindingSites[i].cisregID][indiv.allBindingSites[i].geneCopy]);
    printf(" transcription-factor: %3d", indiv.allBindingSites[i].tfID);
    //printf(" (sequence: %.*s)\n", TF_ELEMENT_LEN, transcriptionFactorSeq[indiv.allBindingSites[i].tfID][indiv.allBindingSites[i].geneCopy]); 
    printf("  L-edge of %2dbp hind: %3d\n", HIND_LENGTH, indiv.allBindingSites[i].leftEdgePos);        
    //printf("  Hind offset position: %3d\n", indiv.allBindingSites[i].hindPos); 
    printf("               strand: %3d\n", indiv.allBindingSites[i].strand);
    //printf("         Hamming dist: %3d\n", indiv.allBindingSites[i].hammingDist); 
  }*/
   
  // system("PAUSE");
   
     printf("\n");
     int lem;
      int bob;
     for(lem =0; lem<TFBS; lem++){
             startPos[lem] = indiv.allBindingSites[lem].leftEdgePos;
             hammDist[lem] = indiv.allBindingSites[lem].hammingDist;
             TFon[lem] = indiv.allBindingSites[lem].tfID;
             bob = indiv.allBindingSites[lem].tfID;
             //printf("bob= %d\n", bob);
             Kon[lem] = initProteinConc[bob]*kon;
             
             printf("%d", indiv.allBindingSites[lem].leftEdgePos);
             printf(" Hd = %d   tf = %d", hammDist[lem], TFon[lem]);
             printf(" Kon[lem] = %f\n", Kon[lem]);
            
     }
     /*int mat;
     for(mat=0; mat<TFBS; mat++){
        printf("tf[%d] = %d\n", mat, TFon[mat]);
     }*/
     printf("\n");
     configure(0,bits,&array,viableStates,TFBS,startPos);
     printf("%d\n", array);
      /* for(mat=0; mat<TFBS; mat++){
        printf("tf[%d] = %d\n", mat, TFon[mat]);
     }*/
     system("PAUSE");
     transitions(array,viableStates,TFBS,arrayT, Kon, Koff, hammDist, diag, TFon);
     //printf("%f\n", *arrayT[0].col[0].kval);
      // system("PAUSE");
     print_arrayT(arrayT,array,viableStates);
     print_arrayT_MATLAB(arrayT,array,viableStates);
       system("PAUSE");
  /* free dynamically allocated all binding sites list */
  free(indiv.allBindingSites);
   int d;
     for (d=0; d<array; d++) {
       free(arrayT[d].col);
  }   
  
  free(arrayT);
  free(viableStates);
  free(bits);
  free(startPos);
  free(hammDist);
  free(TFon);
  /* close error file */
  fclose(fperrors);
}
