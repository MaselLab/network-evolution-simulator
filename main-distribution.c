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


#include "engine.h"

#include "random.h"
#include "lib.h"
#include "netsim.h"

#define BUFSIZE 250
#define NITER 3000

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
   // if(leftEdge < CISREG_LEN - HIND_LENGTH){
    right = leftEdge + HIND_LENGTH;
    int num;
    num =0;
    while( leftPos[num] < right){
           num++;
           }
   //printf("num = %d, leftEdge = %d \n\n", num, leftPos[num]);
   return num;
  //}else{ return 0;}
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
     /* case 'w':
        numBp = atoi(optarg);
        break;
      case 'a':
        array_size = atoi(optarg);
        break;  */      
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
  /*printf("\n");
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = initProteinConc[i]* (rand()%1000)/1000;
    printf("%f\n", initProteinConc[i]);
  }*/
  //system("PAUSE");
  printf("\n");
  
  /* new API requires a Genotype clone */
  Genotype *UNUSED_clone = NULL;
 
  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, UNUSED_clone, kdis, 0);
  
  /*Sort binding sites from smallest left_edge_position to largest left_edge_position*/
  qsort((void *) &(indiv.all_binding_sites[0]), indiv.sites_per_gene[0],                                 
            sizeof(struct AllTFBindingSites),(compfn)compare );
 
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
  
  int startSite, TFBS, startNum, lem, bob, posNext, count, jo, jojo;
  int tfCount, siteCount, size, val, n, b, k, kTF, add, p, chek;
  int *left_edge_pos, *startPos, *hammDist, *TFon; 
  float *weight, *partition;
  float check, checkP;
  
  int unboundCount;
  unboundCount=0;
 
  left_edge_pos = malloc(indiv.sites_per_gene[0]*sizeof(int));
  
  startNum = indiv.all_binding_sites[0].left_edge_pos;
 
  while( indiv.all_binding_sites[TFBS].left_edge_pos < startNum + HIND_LENGTH){
     TFBS++;
  }
  
  posNext=0;
  for( count =0; count<indiv.sites_per_gene[0]; count++){
     left_edge_pos[count] = indiv.all_binding_sites[count].left_edge_pos;
  }

  for( jo=0; jo<3; jo++){
    Gibbs = (((float) jo)/3.0 - 1.0) * RTlnKr;
    Koff[jo] = -Gibbs;
  }
  
  for(jojo=0; jojo<3; jojo++){
    //printf("Koff[%d] = %f\n", jojo, Koff[jojo]);
  }
  
  siteCount = TFBS;
  //printf("siteCount= %d\n", siteCount);
   
  A=R=0.;        
  startSite = 0;
  size =0;
  val =0;
 
  srand(time(NULL));
   recordFile = fopen("recodeFile.txt", "w");
  if ((recordFile = fopen("recordFile.txt", "w"))) {
  while(val < NITER){
    A=R=0.;
    startNum = indiv.all_binding_sites[0].left_edge_pos;
    while( indiv.all_binding_sites[TFBS].left_edge_pos < startNum + HIND_LENGTH){TFBS++;}
    posNext=0;
    A=R=0.;        
    startSite = 0;
    while(left_edge_pos[posNext]<= (CISREG_LEN - HIND_LENGTH)){
      fprintf(recordFile, "left_edge_pos = %d\n", left_edge_pos[posNext]);
      
      tfCount = posNext;
      TFBS =0;
      //printf("startNum = %d\n", tfCount);
      //printf("LEP+HIND_LENGTH = %d\n", (indiv.all_binding_sites[posNext].left_edge_pos) + HIND_LENGTH);
      while(indiv.all_binding_sites[tfCount].left_edge_pos < (indiv.all_binding_sites[posNext].left_edge_pos) + HIND_LENGTH && tfCount < indiv.sites_per_gene[0]){
         //printf("tfcount = %d\n", tfCount);                                             
         TFBS++;
         tfCount++;
      }
      //printf("TFBS=%d\n", TFBS);
  
      arrayWT = malloc((TFBS+1) *sizeof(struct Wtype));
      // system("PAUSE");
   
      for (lem =0; lem<TFBS; lem++) {
         arrayWT[lem].tfbsNum = lem;
         arrayWT[lem].startPos = indiv.all_binding_sites[lem+posNext].left_edge_pos;
         arrayWT[lem].hammDist = indiv.all_binding_sites[lem+posNext].hamming_dist;
         arrayWT[lem].tfIDon = indiv.all_binding_sites[lem+posNext].tf_id;
         bob = indiv.all_binding_sites[lem+posNext].tf_id;
         arrayWT[lem].conc = initProteinConc[bob];
         arrayWT[lem].weight = (float)(initProteinConc[bob] * (float)Koff[(indiv.all_binding_sites[lem+posNext].hamming_dist)]);
 
         fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d Kon = %f weight = %.2f\n",arrayWT[lem].tfbsNum, arrayWT[lem].startPos, arrayWT[lem].hammDist, arrayWT[lem].tfIDon, arrayWT[lem].conc, arrayWT[lem].weight);
      }
      
      qsort((void *) &(arrayWT[0]), TFBS, sizeof(struct Wtype), (compfn)compareWeights);
      
      arrayWT[TFBS].conc = 0;
      arrayWT[TFBS].hammDist = 0;
      arrayWT[TFBS].startPos = 299;
      arrayWT[TFBS].tfbsNum = 111;
      arrayWT[TFBS].tfIDon = 11;
      arrayWT[TFBS].weight = arrayWT[0].weight;
 
      fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d conc = %f weight = %.2f\n",arrayWT[TFBS].tfbsNum, arrayWT[TFBS].startPos, arrayWT[TFBS].hammDist, arrayWT[TFBS].tfIDon, arrayWT[TFBS].conc, arrayWT[TFBS].weight);
 
      float weightSum = 0.;
      for(lem=0; lem<(TFBS); lem++){
         //printf(" weight = %.2f\n", arrayWT[lem].weight);  
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
         //printf("prob[%d] = %f\n", lem, prob[lem]);
      }
     
      fprintf(recordFile, "BEFORE n=%d\n",n);
      
      fprintf(recordFile, "gasdev(&seed) = %f\n", gasdev(&seed));
      n = rand()%1000;
      fprintf(recordFile, "n1=%d\n",n);
      fprintf(recordFile, "n2=%d\n",n);
      check =0.;
      check = n/1000.;
      fprintf(recordFile, "\ncheck = %f\n", check);
      partition = malloc((TFBS+1)*sizeof(float));
      for(lem=0; lem<TFBS+1; lem++){
         if(lem==0){partition[lem] = prob[lem];}
         else if(lem!=TFBS){ partition[lem] = (partition[(lem-1)]+prob[lem]);}
         else {partition[lem] = 1;}
         fprintf(recordFile, "partition[%d] = %f\n", lem, partition[lem]);
      }
      b = 0;
      lem =0;
     checkP =check;
      fprintf(recordFile, "checkP=%f\n", checkP);
      
      kTF =0;
      if(0<=checkP && checkP<=partition[0]){b=0;}
      else{ for(k=0; k<TFBS; k++){
              fprintf(recordFile, "part[%d] = %.2f, part[%d] = %.2f\n",k,  partition[k],k+1,partition[k+1]);
              if(partition[k]<checkP && checkP<=partition[k+1]){
                 b=k+1;
                 kTF=1;
              }
            }
            if(kTF ==0){b=TFBS;}
      }     
      
      fprintf(recordFile, "b=%d\n", b);
      fprintf(recordFile, "TF = %d \n", arrayWT[b].tfIDon);
      fprintf(recordFile, "lem = %d\n", arrayWT[b].tfbsNum);
      fprintf(recordFile, "activating[][] = %d\n", indiv.activating[arrayWT[b].tfIDon][0]);
      
      if(arrayWT[b].tfIDon == 11){/*do nothing, there is nothing bound */
         unboundCount++;}
      else{
        if(indiv.activating[arrayWT[b].tfIDon][indiv.all_binding_sites[(arrayWT[b].tfbsNum)].gene_copy] ==1){A++;}
        else{R++;}
      }
      fprintf(recordFile, "A = %d   R = %d\n", A, R);
      if(A+R > 9){
             fprintf(recordFile, "PROBLEM!!!!!!!");
             system("PAUSE");
             //TO DO: fix this problem!! Too many things get bound, not sure what is wrong!
      }
 
      startSite++;
      fprintf(recordFile, "startSite = %d\n", startSite);
      fprintf(recordFile, "ARGH!! = %d\n", arrayWT[b].startPos);
      
      if(arrayWT[b].startPos != 299){
        startNum = arrayWT[b].startPos + HIND_LENGTH;
        
         fprintf(recordFile, "startNum = %d\n", startNum);

        posNext = nextPos(arrayWT[b].startPos, left_edge_pos);
      }else{
            posNext++;
      }
      fprintf(recordFile, "posNext = %d\n\n", posNext);
    }
    
    fprintf(recordFile, "\n\nEND OF GENE!!!!!!!!\n\n");
    fprintf(recordFile, "A = %d   R = %d\n", A, R);
   
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
  chek =0;
  for(chek=0; chek<size; chek++){
     printf("%d  count = %d  active = %.2f, repress = %.2f, ratio = %.3f\n", chek, arrayD[chek].count, arrayD[chek].active, arrayD[chek].repress, arrayD[chek].ratio);
     fprintf(recordFile, "%d  count = %d  active = %.2f, repress = %.2f, ratio = %.3f\n", chek, arrayD[chek].count, arrayD[chek].active, arrayD[chek].repress, arrayD[chek].ratio);
  }
  fprintf(recordFile, "unboundCount = %d\n", unboundCount);
  }//file for-loop
  fclose(recordFile);
  printf("Here\n");  
  system("PAUSE");	
    
  free(left_edge_pos);
  free(startPos);
  free(hammDist);
  free(TFon);
  free(weight);
  free(arrayD);
  free(arrayWT); 
}
