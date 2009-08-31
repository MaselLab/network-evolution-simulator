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

void active_to_repress(Genotype indiv, float initProteinConc[NGENES]){
      /*each struct holds # activators, # repressors, ratio A/R, and how many times ratio has occurred */          
 struct Dtype *arrayD;
 arrayD = malloc(NITER*sizeof(struct Dtype));
 
 /*each struct holds ID, BS numer, start positon, hamming distance, concentration and weight*/
 struct Wtype *arrayWT;
  int A, R;
          
  float Koff[5]; //rethink 5?
  float Gibbs, RTlnKr;
  float temperature = 293.0;
  RTlnKr = GASCONSTANT * temperature * log(KR);
  
  int startSite, TFBS, startNum, lem, bob, posNext, count, activeCount;
  int tfCount, siteCount, size, val, n, b, k, kTF, add, p;
  int *left_edge_pos; 
  float *partition, *active_prob;
  float check, checkP, percent;
  
  int unboundCount;
  unboundCount=0;
  activeCount = 0;
  percent = 0.000;
  active_prob = malloc(TFGENES*sizeof(float));
 // int bsNum;
  
  int gene_num;
  //int number =0;
  
  gene_num =0;
  recordFile = fopen("recodeFile.txt", "w");
  if ((recordFile = fopen("recordFile.txt", "w"))) {
  /*for(bsNum = 0; bsNum < indiv.binding_sites_num; bsNum++){
    for(gene_num =0; gene_num<TFGENES; gene_num++){
       int bsBegin = indiv.sites_per*/
  percent = 0.000;
  unboundCount=0;
  activeCount = 0;
  percent = 0.000;
  left_edge_pos = malloc(indiv.sites_per_gene[0]*sizeof(int));
  
  printf("NUM BS = %d\n", indiv.binding_sites_num);
  
  int *gene_start_pos;
  gene_start_pos = malloc(TFGENES*sizeof(int));
  int position =0;
  gene_start_pos[0] = 0;
  for(lem=0; lem < TFGENES; lem++){
             //printf("%d  %d  \n", lem, indiv.sites_per_gene[lem]);
            // position = indiv.sites_per_gene[lem];
             if(lem==0){ gene_start_pos[lem] = 0;
                         position = indiv.sites_per_gene[lem];}
             else{
             gene_start_pos[lem] = position;
            // printf("%d  %d  ", lem, indiv.sites_per_gene[lem]);
             position+=indiv.sites_per_gene[lem];
             
             
             }
             printf("%d %d\n",lem, gene_start_pos[lem]);
  }
  printf("sites_per_gene[] = %d\n TFGENES = %d\n", indiv.sites_per_gene[TFGENES], TFGENES);
  system("PAUSE");
  //TO DO: all_binding_sites holds ALL BS. need to divide up into genes. where gene[0] ends, gene[1] begins.
  //need to loop through all genes and store probabilities for each gene.
  /* if(gene_num == 0){number =0;}
   else{ number += indiv.sites_per_gene[gene_num-1];}
   
  printf("number = %d\n", number);
  int otherNum = indiv.all_binding_sites[number].left_edge_pos;
  printf("otherNum = %d\n", otherNum);*/
  for(gene_num =0; gene_num < NGENES; gene_num++){
  int start_pos = gene_start_pos[gene_num];
  startNum = indiv.all_binding_sites[start_pos].left_edge_pos;
  printf("startNum = %d\n", startNum);
  TFBS=0;
  printf("TFBS = %d\n", TFBS);
  while( indiv.all_binding_sites[TFBS+start_pos].left_edge_pos < startNum + HIND_LENGTH){
     TFBS++;
  }
  printf("TFBS = %d\n", TFBS);
  system("PAUSE");
  posNext=0;
  for( count =0; count<indiv.sites_per_gene[gene_num]; count++){
     left_edge_pos[count] = indiv.all_binding_sites[count+start_pos].left_edge_pos;
  }
  printf("count = %d\n", count);
  for( lem=0; lem<3; lem++){
    Gibbs = (((float) lem)/3.0 - 1.0) * RTlnKr;
    Koff[lem] = -Gibbs;
  }
 //system("PAUSE"); 
  siteCount = TFBS;
   
  A=R=0.;        
  startSite = 0;
  size =0;
  val =0;
 
 
  srand(time(NULL));
   //recordFile = fopen("recodeFile.txt", "w");
  //if ((recordFile = fopen("recordFile.txt", "w"))) {
  startNum = indiv.all_binding_sites[start_pos].left_edge_pos;
  printf("val=%d, NITER = %d\n", val, NITER);
  while(val < NITER){
    A=R=0.;
    
    fprintf(recordFile, "startNum = %d\n", startNum);
    //printf("TFBS!!! = %d\n", TFBS);
    //system("PAUSE");
    //TFBS=0;
    //while( indiv.all_binding_sites[TFBS].left_edge_pos < startNum + HIND_LENGTH){TFBS++;}
    posNext=0;
    A=R=0.;        
    startSite = 0;
    while(left_edge_pos[posNext]<= (CISREG_LEN - HIND_LENGTH)){
       fprintf(recordFile,  "left_edge_pos[posNext] = %d\n", left_edge_pos[posNext]);
     // fprintf(recordFile, "left_edge_pos = %d\n", left_edge_pos[posNext]);
      
      tfCount = posNext;
      TFBS =0;
      while(indiv.all_binding_sites[start_pos + tfCount].left_edge_pos < (indiv.all_binding_sites[posNext].left_edge_pos) + HIND_LENGTH && tfCount < indiv.sites_per_gene[0]){                                         
         TFBS++;
         tfCount++;
      }
  
      arrayWT = malloc((TFBS+1) *sizeof(struct Wtype));
   
      for (lem =0; lem<TFBS; lem++) {
         arrayWT[lem].tfbsNum = lem+start_pos;
         arrayWT[lem].startPos = indiv.all_binding_sites[start_pos+lem+posNext].left_edge_pos;
         arrayWT[lem].hammDist = indiv.all_binding_sites[start_pos+lem+posNext].hamming_dist;
         arrayWT[lem].tfIDon = indiv.all_binding_sites[start_pos+lem+posNext].tf_id;
         bob = indiv.all_binding_sites[start_pos+lem+posNext].tf_id;
         arrayWT[lem].conc = initProteinConc[bob];
         arrayWT[lem].weight = (float)(initProteinConc[bob] * (float)Koff[(indiv.all_binding_sites[start_pos+lem+posNext].hamming_dist)]);
 
         //fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d Kon = %f weight = %.2f\n",arrayWT[lem].tfbsNum, arrayWT[lem].startPos, arrayWT[lem].hammDist, arrayWT[lem].tfIDon, arrayWT[lem].conc, arrayWT[lem].weight);
      }
      
      qsort((void *) &(arrayWT[0]), TFBS, sizeof(struct Wtype), (compfn)compareWeights);
      
      arrayWT[TFBS].conc = 0;
      arrayWT[TFBS].hammDist = 0;
      arrayWT[TFBS].startPos = 299;
      arrayWT[TFBS].tfbsNum = 111;
      arrayWT[TFBS].tfIDon = 11;
      arrayWT[TFBS].weight = arrayWT[0].weight;
 
      //fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d conc = %f weight = %.2f\n",arrayWT[TFBS].tfbsNum, arrayWT[TFBS].startPos, arrayWT[TFBS].hammDist, arrayWT[TFBS].tfIDon, arrayWT[TFBS].conc, arrayWT[TFBS].weight);
 
      float weightSum = 0.;
      for(lem=0; lem<(TFBS); lem++){
         weightSum = (float)weightSum +  (float)arrayWT[lem].weight;
        // printf("       weightSum = %f\n", weightSum);
               /*FIX THIS! SUM IS NOT CORRECT! ROUNDING ERRORS!*/
         
      }
      weightSum += arrayWT[0].weight;
    
      float *prob;
      prob = malloc((TFBS+1)*sizeof(float));
    
      for(lem =0; lem<TFBS+1; lem++){
         if(lem == TFBS) { prob[lem] = (arrayWT[0].weight) / weightSum;}
         else {prob[lem] = (arrayWT[lem].weight) / weightSum;}
      }
     
     // fprintf(recordFile, "BEFORE n=%d\n",n);
      
      //fprintf(recordFile, "gasdev(&seed) = %f\n", gasdev(&seed));
      n = rand()%1000;
      //fprintf(recordFile, "n1=%d\n",n);
     // fprintf(recordFile, "n2=%d\n",n);
      check =0.;
      check = n/1000.;
     // fprintf(recordFile, "\ncheck = %f\n", check);
      partition = malloc((TFBS+1)*sizeof(float));
      for(lem=0; lem<TFBS+1; lem++){
         if(lem==0){partition[lem] = prob[lem];}
         else if(lem!=TFBS){ partition[lem] = (partition[(lem-1)]+prob[lem]);}
         else {partition[lem] = 1;}
        // fprintf(recordFile, "partition[%d] = %f\n", lem, partition[lem]);
      }
      b = 0;
      lem =0;
     checkP =check;
     // fprintf(recordFile, "checkP=%f\n", checkP);
      
      kTF =0;
      if(0<=checkP && checkP<=partition[0]){b=0;}
      else{ for(k=0; k<TFBS; k++){
             // fprintf(recordFile, "part[%d] = %.2f, part[%d] = %.2f\n",k,  partition[k],k+1,partition[k+1]);
              if(partition[k]<checkP && checkP<=partition[k+1]){
                 b=k+1;
                 kTF=1;
              }
            }
            if(kTF ==0){b=TFBS;}
      }     
      
     // fprintf(recordFile, "b=%d\n", b);
     // fprintf(recordFile, "TF = %d \n", arrayWT[b].tfIDon);
     // fprintf(recordFile, "lem = %d\n", arrayWT[b].tfbsNum);
     // fprintf(recordFile, "activating[][] = %d\n", indiv.activating[arrayWT[b].tfIDon][0]);
      
      if(arrayWT[b].tfIDon == 11){/*do nothing, there is nothing bound */
         unboundCount++;}
      else{
        if(indiv.activating[arrayWT[b].tfIDon][indiv.all_binding_sites[(arrayWT[b].tfbsNum)].gene_copy] ==1){A++;}
        else{R++;}
      }
      //fprintf(recordFile, "A = %d   R = %d\n", A, R);
      if(A+R > 9){
             fprintf(recordFile, "PROBLEM!!!!!!!");
             system("PAUSE");
             //TO DO: fix this problem!! Too many things get bound, not sure what is wrong!
      }
 
      startSite++;
      //fprintf(recordFile, "startSite = %d\n", startSite);
     // fprintf(recordFile, "ARGH!! = %d\n", arrayWT[b].startPos);
      
      if(arrayWT[b].startPos != 299){
        startNum = arrayWT[b].startPos + HIND_LENGTH;
        
         fprintf(recordFile, "startNum = %d\n", startNum);

        posNext = nextPos(arrayWT[b].startPos, left_edge_pos);
        fprintf("%d\n", posNext);
      }else{
            posNext++;
      }
      //fprintf(recordFile, "posNext = %d\n\n", posNext);
    }
    
    //fprintf(recordFile, "\n\nEND OF GENE!!!!!!!!\n\n");
   // fprintf(recordFile, "A = %d   R = %d\n", A, R);
   
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
  } 
  activeCount=0;
  for(lem=0; lem<size; lem++){
     printf("%d  count = %d  active = %.2f, repress = %.2f, ratio = %.3f\n", lem, arrayD[lem].count, arrayD[lem].active, arrayD[lem].repress, arrayD[lem].ratio);
     fprintf(recordFile, "%d  count = %d  active = %.2f, repress = %.2f, ratio = %.3f\n", lem, arrayD[lem].count, arrayD[lem].active, arrayD[lem].repress, arrayD[lem].ratio);
     if(arrayD[lem].repress < arrayD[lem].ratio){
         activeCount += arrayD[lem].count;} 
  }
  percent = activeCount/(float)NITER;
  fprintf(recordFile, "unboundCount = %d\n", unboundCount);
  fprintf(recordFile, "activeCount = %d\n", activeCount);
  fprintf(recordFile, "percent = %f\n", percent);
  //active_prob[gene_num] = percent;
  //}
    /*for(lem =0; lem<TFGENES; lem ++){
          printf("%.4f ", active_prob[lem]);
          fprintf(recordFile, "%.4f ", active_prob[lem]);
    }*/
    printf("here");
}
  }//file for-loop
  
  fclose(recordFile);
  
  free(left_edge_pos); 
  free(arrayD);
  free(arrayWT);
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


  /**********************************************************************
   *  modifiable part starts here
   **********************************************************************/

  /* Initialize protein and mRNA concentrations */
  /* Jasmin: you can use these initial concentrations for computing your kons */
  
  /*Creates array of inial protein concentrations for computing kons*/
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    printf("%f\n", initProteinConc[i]);
  }
  printf("\n");
  
  /* new API requires a Genotype clone */
  Genotype *UNUSED_clone = NULL;
 
  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, UNUSED_clone, kdis, 0);
  
  /*Sort binding sites from smallest left_edge_position to largest left_edge_position*/
  qsort((void *) &(indiv.all_binding_sites[0]), indiv.sites_per_gene[0],                                 
            sizeof(struct AllTFBindingSites),(compfn)compare );
 
 active_to_repress(indiv, initProteinConc);
 
  printf("Here\n");  
  system("PAUSE");	

}
