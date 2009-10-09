/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Computes probability cloud of occupancy of transcription factors in cis-regulatory regions
 * Authors: Jasmin Uribe, Alex Lancaster
 * Copyright (c) 2008, 2009 Arizona Board of Regents (University of Arizona)
 */

#ifdef __WIN32__
#include <windows.h>
#endif

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

#define BUFSIZE 250
#define NITER 1000

FILE *recordFile;
FILE *outputFile;

struct Dtype{
       float active;
       float repress;
       float ratio;
       int count;
};

struct Wtype{
       int tfbsNum;
       int startPos;
       int hammDist;
       int tfIDon;
       float conc;
       float weight;
};    

typedef int (*compfn)(const void*, const void*);
//typedef float (*compfnf)(const void*, const void*);

int intcmp(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

/*Reorders binding sites from smallest left_edge_position to largest left_edge_position */
int compare(struct AllTFBindingSites  *elem1, struct AllTFBindingSites *elem2){
    if(elem1->left_edge_pos < elem2->left_edge_pos)
        return -1;
    else if(elem1->left_edge_pos > elem2->left_edge_pos)
         return 1;
    else
         return 0;
}

int compareWeights(struct Wtype *elem1, struct Wtype *elem2){
      if(elem1->weight < elem2->weight)
          return -1;
      else if(elem1->weight > elem2->weight)
          return 1;
      else 
          return 0;
}

//given current bound position get next avaliable position
int nextPos( int leftEdge, int *leftPos){
 
    int right;
    right = leftEdge + HIND_LENGTH;
    int num;
    num =0;
    while(leftPos[num] < right) {num++;}
   return num;
}

float active_to_repress(Genotype indiv, float initProteinConc[NGENES], int start, int gene){
       int n;
       
       int  TFBS, startNum, lem, bob, posNext, activeCount;
  int  size, val, m, b, k, kTF, add, p;
  int *left_edge_pos; 
  float *partition;
  float check, checkP, percent;
  float weightSum;
  int structStart;
       
       qsort((void *) &(indiv.all_binding_sites[start]), indiv.sites_per_gene[gene],                                 
            sizeof(struct AllTFBindingSites),(compfn)compare );
       
      recordFile = fopen("recordFile.txt", "w");
      if ((recordFile = fopen("recordFile.txt", "w"))) {
     
      fprintf(recordFile, "BS num  = %d\n", indiv.binding_sites_num);
       fprintf(recordFile, "BS num gene 0 = %d\n", indiv.sites_per_gene[0]);
      for(n=0; n<indiv.sites_per_gene[gene]; n++){
               fprintf(recordFile, "%d   %d\n",n,  indiv.all_binding_sites[start+n].left_edge_pos);
      }
      
      
      struct Dtype *arrayD;
 arrayD = malloc(NITER*sizeof(struct Dtype));
 
 /*each struct holds ID, BS numer, start positon, hamming distance, concentration and weight*/
 struct Wtype *arrayWT;
  int A, R;
          
  float Koff[5]; //rethink 5?
  float Gibbs, RTlnKr;
  float temperature = 293.0;
  RTlnKr = GASCONSTANT * temperature * log(KR);
  
  //Koff prob according to Gibbs
  for( n=0; n<3; n++){
    Gibbs = (((float) n)/3.0 - 1.0) * RTlnKr;
    Koff[n] = -Gibbs;
  }
  
  int unboundCount;
  unboundCount=0;
  activeCount = 0;
  percent = 0.000;
  left_edge_pos = malloc(indiv.sites_per_gene[gene]*sizeof(int));
  
  //populate left_edge_pos
  for(n=0; n<indiv.sites_per_gene[gene]; n++){
           left_edge_pos[n] = indiv.all_binding_sites[start+n].left_edge_pos;
  }
  
  startNum = indiv.all_binding_sites[start].left_edge_pos;
  int endNum =indiv.sites_per_gene[gene]+start;
 
  //count num BS in window
  TFBS=0;
  while( indiv.all_binding_sites[start+TFBS].left_edge_pos < startNum + HIND_LENGTH){
     TFBS++;
  }
  srand(time(NULL));  
    
    posNext=0;
    structStart = start;
    A=R=0.; 
    val =0;  
    size=0;
    while(val < NITER){
      weightSum = 0.;
      A=R=0.;
      posNext=0;
      structStart = start;
      //fprintf(recordFile, "\n\n -----ITER = %d ----- \n\n", val);
      while(left_edge_pos[posNext]<= (CISREG_LEN - HIND_LENGTH)){//loop through all sites in window
          structStart = start + posNext;              
         TFBS=0;
       
         while( indiv.all_binding_sites[structStart+TFBS].left_edge_pos < indiv.all_binding_sites[structStart].left_edge_pos + HIND_LENGTH && 
                       structStart+TFBS < endNum){
            TFBS++;
           
         }
         //fprintf(recordFile, "TFBS=%d\n", TFBS);
        arrayWT = malloc((TFBS+1) *sizeof(struct Wtype));
   
         for (lem =0; lem<TFBS; lem++) {
           arrayWT[lem].tfbsNum = lem+structStart;
           arrayWT[lem].startPos = indiv.all_binding_sites[lem+structStart].left_edge_pos;
           arrayWT[lem].hammDist = indiv.all_binding_sites[lem+structStart].hamming_dist;
           arrayWT[lem].tfIDon = indiv.all_binding_sites[lem+structStart].tf_id;
           bob = indiv.all_binding_sites[lem+structStart].tf_id;
           arrayWT[lem].conc = initProteinConc[bob];
           arrayWT[lem].weight = (float)(initProteinConc[bob] * (float)Koff[(indiv.all_binding_sites[lem+structStart].hamming_dist)]);
 
        // fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d Kon = %f weight = %.2f\n",arrayWT[lem].tfbsNum, arrayWT[lem].startPos, arrayWT[lem].hammDist, arrayWT[lem].tfIDon, arrayWT[lem].conc, arrayWT[lem].weight);
      }//closes array for-loop
      
      qsort((void *) &(arrayWT[0]), TFBS, sizeof(struct Wtype), (compfn)compareWeights);
      
      arrayWT[TFBS].conc = 0;
      arrayWT[TFBS].hammDist = 0;
      arrayWT[TFBS].startPos = 299;
      arrayWT[TFBS].tfbsNum = 111;
      arrayWT[TFBS].tfIDon = 11;
      arrayWT[TFBS].weight = arrayWT[0].weight;
      
      weightSum =0;
      for(n=0; n<(TFBS); n++){
         weightSum = (float)weightSum +  (float)arrayWT[n].weight;
        // printf("       weightSum = %f\n", weightSum);
               /*FIX THIS! SUM IS NOT CORRECT! ROUNDING ERRORS!*/    
      }
      weightSum += arrayWT[0].weight;
      
      float *prob;
      prob = malloc((TFBS+1)*sizeof(float));
    
      for(n =0; n<TFBS+1; n++){
         if(n == TFBS) { prob[n] = (arrayWT[0].weight) / weightSum;}
         else {prob[n] = (arrayWT[n].weight) / weightSum;}
      }
      m = rand()%1000;
      check =0.;
      check = m/1000.;
      partition = malloc((TFBS+1)*sizeof(float));
      for(n=0; n<TFBS+1; n++){
         if(n==0){partition[n] = prob[n];}
         else if(n!=TFBS){ partition[n] = (partition[(n-1)]+prob[n]);}
         else {partition[n] = 1;}
      }
      b = 0;
      n =0;
     checkP =check;
      
      kTF =0;
      if(0<=checkP && checkP<=partition[0]){b=0;}
      else{ for(k=0; k<TFBS; k++){
              if(partition[k]<checkP && checkP<=partition[k+1]){
                 b=k+1;
                 kTF=1;
              }
            }
            if(kTF ==0){b=TFBS;}
      } 
      if(arrayWT[b].tfIDon == 11){/*do nothing, there is nothing bound */
         unboundCount++;}
      else{
        if(indiv.activating[arrayWT[b].tfIDon][indiv.all_binding_sites[(arrayWT[b].tfbsNum)].gene_copy] ==1){A++;}
        else{R++;}
      }
      if(A+R > 9){
             fprintf(recordFile, "PROBLEM!!!!!!!");
             system("PAUSE");
             //TO DO: fix this problem!! Too many things get bound, not sure what is wrong!
      }
      fprintf(recordFile, "A=%d, R=%d\n", A, R);  
      
      if(arrayWT[b].startPos != 299){
        startNum = arrayWT[b].startPos + HIND_LENGTH;
        
        posNext = nextPos(arrayWT[b].startPos, left_edge_pos);
       
      }else{
            posNext++;
      }
      }//closes while(left_edge_pos)
     
      if(val!=0){
      p=0;
      add =1;
      while(p<val){
         if(A == arrayD[p].active){ 
            if(R == arrayD[p].repress){
               arrayD[p].count = arrayD[p].count +1;
               add =0;
               break;
            }
         }
         p++;
      }
    }
    if(val ==0 || add ==1){ 
       arrayD[size].active = A;
       arrayD[size].repress = R;
       arrayD[size].ratio = arrayD[size].active * .33 + .31;
       arrayD[size].count = 1;
       size++;
    }
    val++;
      }//closes while(val-NITER)
      
      activeCount=0;
      printf("%f\n", arrayD[0].active);
      
      for(lem=0; lem<size; lem++){
                
      if(arrayD[lem].repress < arrayD[lem].ratio){
         activeCount += arrayD[lem].count;} 
      }
       printf("active= %d, repress= %f, ratio=%f, count = %d\n", arrayD[0].active, arrayD[0].repress, arrayD[0].ratio, arrayD[0].count);
      percent = activeCount/(float)NITER;
      fprintf(recordFile, "unboundCount = %d\n", unboundCount);
      fprintf(recordFile, "activeCount = %d\n", activeCount);
      fprintf(recordFile, "percent = %f\n", percent);
      
      }//closes recordFile loop
      fclose(recordFile);
      return percent;
 }
 
 void active_vect(Genotype indiv, float initProteinConc[NGENES], float *gene_active){
       int geneSum = 0;
       int i;
       for(i=0; i<NGENES; i++){ 
           printf("geneSum = %d,  i=%d\n", geneSum, i);
           gene_active[i] = active_to_repress(indiv, initProteinConc, geneSum, i);
           fprintf(outputFile, "geneSum = %d, i = %d, WHAT = %f\n\n", geneSum, i, gene_active[i] );
           geneSum += indiv.sites_per_gene[i];
        
  }
 }
     

int main(int argc, char *argv[])
{
  FILE *fpkdis;
  char fperrors_name[80];
  int i, j;
  Genotype indiv;
  float initProteinConc[NGENES], kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int curr_seed;
  
#ifdef __unix__
  const char *startcmd = "/usr/local/MathSoft/bin/matlab -nojvm -nosplash -nodisplay -nodesktop";
#else 
#ifdef __WIN32__
    const char *startcmd = NULL;
#endif
#endif
  
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
  outputFile = fopen("outputFile.txt", "w");
      if ((outputFile = fopen("outputFile.txt", "w"))) {


  /**********************************************************************
   *  modifiable part starts here
   **********************************************************************/

  /* Initialize protein and mRNA concentrations */
  /* Jasmin: you can use these initial concentrations for computing your kons */
  
  /*Creates array of inial protein concentrations for computing kons*/
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    printf("%f\n", initProteinConc[i]);
    fprintf(outputFile, "%f\n", initProteinConc[i]);
  }
  printf("\n");
  
  /* new API requires a Genotype clone */
  Genotype *UNUSED_clone = NULL;
 
  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, UNUSED_clone, kdis, 0);
  
  /*Sort binding sites from smallest left_edge_position to largest left_edge_position*/
  
  float *gene_active;
  gene_active = malloc(NGENES*sizeof(float));
  printf("OUTSIDE WHAT = %.2f\n\n", active_to_repress(indiv, initProteinConc, 0, 0));
  active_vect(indiv, initProteinConc, gene_active);
  for(i=0; i<NGENES; i++){
           fprintf(outputFile, " WHAT = %f\n", gene_active[i] );
          /* printf("geneSum = %d,  i=%d\n", geneSum, i);
           gene_active[i] = active_to_repress(indiv, initProteinConc, geneSum, i);
           fprintf(outputFile, "geneSum = %d, i = %d, WHAT = %f\n\n", geneSum, i, gene_active[i] );
           geneSum += indiv.sites_per_gene[i];*/
        
  }
 }//closes file loop
  fclose(outputFile);

  system("PAUSE");	

}
