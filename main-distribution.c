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
#define NITER 200

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
  /*for (i=0; i<NGENES; i++) {
    initProteinConc[i] = 1000;
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
  //printf("startNum = %d\n", startNum);
 
  while( indiv.all_binding_sites[TFBS].left_edge_pos < startNum + HIND_LENGTH){
     TFBS++;
  }
  //printf("TFBS=%d\n", TFBS);
  
  posNext=0;
  for( count =0; count<indiv.sites_per_gene[0]; count++){
     left_edge_pos[count] = indiv.all_binding_sites[count].left_edge_pos;
     //printf("%d\n", left_edge_pos[count]);
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
  n= time(NULL);
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
 
         /*startPos[lem] = indiv.all_binding_sites[lem+posNext].left_edge_pos;
         hammDist[lem] = indiv.all_binding_sites[lem+posNext].hamming_dist;
         TFon[lem] = indiv.all_binding_sites[lem+posNext].tf_id;
         Kon[lem] = initProteinConc[bob];
         weight[lem] = (float)(initProteinConc[bob] * Koff[(indiv.all_binding_sites[lem+posNext].hamming_dist)]);*/
      
         //printf("%d  LEP = %d  Hd = %d   tf = %d Kon = %f weight = %.2f\n",lem, indiv.all_binding_sites[lem+posNext].left_edge_pos, hammDist[lem], TFon[lem], Kon[lem], weight[lem]);
         //printf("%d  LEP = %d  Hd = %d   tf = %d Kon = %f weight = %.2f\n",arrayWT[lem].tfbsNum, arrayWT[lem].startPos, arrayWT[lem].hammDist, arrayWT[lem].tfIDon, arrayWT[lem].conc, arrayWT[lem].weight);
         fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d Kon = %f weight = %.2f\n",arrayWT[lem].tfbsNum, arrayWT[lem].startPos, arrayWT[lem].hammDist, arrayWT[lem].tfIDon, arrayWT[lem].conc, arrayWT[lem].weight);
      }
      
      //system("PAUSE");
      //printf("\n");
      qsort((void *) &(arrayWT[0]), TFBS, sizeof(struct Wtype), (compfn)compareWeights);
      
      //system("PAUSE");
      //printf("last entry\n");
      arrayWT[TFBS].conc = 0;
      arrayWT[TFBS].hammDist = 0;
      arrayWT[TFBS].startPos = 299;
      arrayWT[TFBS].tfbsNum = 111;
      arrayWT[TFBS].tfIDon = 11;
      arrayWT[TFBS].weight = arrayWT[0].weight;
     // printf("%d  LEP = %d  Hd = %d   tf = %d conc = %f weight = %.2f\n",arrayWT[TFBS].tfbsNum, arrayWT[TFBS].startPos, arrayWT[TFBS].hammDist, arrayWT[TFBS].tfIDon, arrayWT[TFBS].conc, arrayWT[TFBS].weight);
      fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d conc = %f weight = %.2f\n",arrayWT[TFBS].tfbsNum, arrayWT[TFBS].startPos, arrayWT[TFBS].hammDist, arrayWT[TFBS].tfIDon, arrayWT[TFBS].conc, arrayWT[TFBS].weight);
      //printf("END last entry\n");
      
      
      float weightSum = 0.;
      for(lem=0; lem<(TFBS); lem++){
         //printf(" weight = %.2f\n", arrayWT[lem].weight);  
         weightSum = (float)weightSum +  (float)arrayWT[lem].weight;
        // printf("       weightSum = %f\n", weightSum);
               /*FIX THIS! SUM IS NOT CORRECT! ROUNDING ERRORS!
         if(lem==0) {weightSum += 2 * arrayWT[lem].weight;}
         else  {weightSum += arrayWT[lem].weight;}*/
      }
      weightSum += arrayWT[0].weight;
      //printf("weightSum = %.2f\n", weightSum);
      //fprintf(recordFile, "weightSum = %.2f\n", weightSum);
    
      float *prob;
      prob = malloc((TFBS+1)*sizeof(float));
    
      for(lem =0; lem<TFBS+1; lem++){
         if(lem == TFBS) { prob[lem] = (arrayWT[0].weight) / weightSum;}
         else {prob[lem] = (arrayWT[lem].weight) / weightSum;}
         //printf("prob[%d] = %f\n", lem, prob[lem]);
      }
      //system("PAUSE");
      fprintf(recordFile, "BEFORE n=%d\n",n);
      
      srand(n);  
      //NOT RANDOM!! FIX THIS!!
      
      n = rand()%100;
      //printf("n=%d\n",n);
      check =0.;
      check = n/1.;
      //printf("\ncheck = %f\n", check);
      fprintf(recordFile, "\ncheck = %f\n", check);
      partition = malloc((TFBS+1)*sizeof(float));
      for(lem=0; lem<TFBS+1; lem++){
         if(lem==0){partition[lem] = prob[lem];}
         else if(lem!=TFBS){ partition[lem] = (partition[(lem-1)]+prob[lem]);}
         else {partition[lem] = 1;}
         //printf("partition[%d] = %f\n", lem, partition[lem]);
         fprintf(recordFile, "partition[%d] = %f\n", lem, partition[lem]);
      }
      b = 0;
      lem =0;
      checkP =0.;
      checkP = check /100.;
      //printf("checkP=%f\n", checkP);
      fprintf(recordFile, "checkP=%f\n", checkP);
     
     /* if(0<=checkP && checkP<=partition[0]){b=0;}
      else if(partition[0]<checkP && checkP<=partition[1]){b=1;}
      else if(partition[1]<checkP && checkP<=partition[2]){b=2;}
      else if(partition[2]<checkP && checkP<=partition[3]){b=3;}
      else if(partition[3]<checkP && checkP<=partition[4]){b=4;}
      else if(partition[4]<checkP && checkP<=partition[5]){b=5;}
      else if(partition[5]<checkP && checkP<=partition[6]){b=6;}
      else {b=7;}*/
      //system("PAUSE");
      
      kTF =0;
      if(0<=checkP && checkP<=partition[0]){b=0;}
      else{ for(k=0; k<TFBS; k++){
              //printf("part[%d] = %.2f, part[%d] = %.2f\n",k,  partition[k],k+1,partition[k+1]);
              fprintf(recordFile, "part[%d] = %.2f, part[%d] = %.2f\n",k,  partition[k],k+1,partition[k+1]);
              if(partition[k]<checkP && checkP<=partition[k+1]){
                 b=k+1;
                 kTF=1;
              }/*else{
                 b=TFBS;
              }*/
            }
            if(kTF ==0){b=TFBS;}
            //else{b=TFBS;}
      }     
      
 
      //printf("b=%d\n", b);
      //printf("TF = %d \n", arrayWT[b].tfIDon);
     // printf("lem = %d\n", arrayWT[b].tfbsNum);
     // printf("activating[][] = %d\n", indiv.activating[arrayWT[b].tfIDon][0]);
      
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
      //printf("A = %d   R = %d\n", A, R);
      fprintf(recordFile, "A = %d   R = %d\n", A, R);
      if(A+R > 9){
            // printf("PROBLEM!!!!!!!");
             fprintf(recordFile, "PROBLEM!!!!!!!");
             system("PAUSE");
             //TO DO: fix this problem!! Too many things get bound, not sure what is wrong!
      }
 
      startSite++;
     // printf("startSite = %d\n", startSite);
     // printf("ARGH!! = %d\n", arrayWT[b].startPos);
      
      fprintf(recordFile, "startSite = %d\n", startSite);
      fprintf(recordFile, "ARGH!! = %d\n", arrayWT[b].startPos);
      
      if(arrayWT[b].startPos != 299){
        startNum = arrayWT[b].startPos + HIND_LENGTH;
       // printf("startNum = %d\n", startNum);
        
         fprintf(recordFile, "startNum = %d\n", startNum);

        posNext = nextPos(arrayWT[b].startPos, left_edge_pos);
      }else{
            posNext++;
      }
      //printf("posNext = %d\n\n", posNext);
      fprintf(recordFile, "posNext = %d\n\n", posNext);
      //system("PAUSE");
    }
    //printf("\n\nEND OF GENE!!!!!!!!\n\n");
   // printf("A = %d   R = %d\n", A, R);
    
    fprintf(recordFile, "\n\nEND OF GENE!!!!!!!!\n\n");
    fprintf(recordFile, "A = %d   R = %d\n", A, R);
    //system("PAUSE");
    
    if(val!=0){
      p=0;
      add =1;
      while(p<val){
         if(A == arrayD[p].active){ 
            //printf("ACITVE MATCH ");
            if(R == arrayD[p].repress){
             //  printf("REPRESS MATCH");
               arrayD[p].count = arrayD[p].count +1;
             //  printf("count = %d\n", arrayD[p].count);
             //  printf("\n");
               add =0;
               break;
            }
         }
         //printf("\n");
         p++;
      }
    }
    if(val ==0 || add ==1){ 
       arrayD[size].active = A;
       arrayD[size].repress = R;
       arrayD[size].ratio = arrayD[size].active * .33442 + .31303;
       /*if(R==0) {arrayD[size].ratio = A*10.;}
       else{arrayD[size].ratio = arrayD[size].active * .33442 + .31303;}*/
       arrayD[size].count = 1;
       size++;
    }

    val++;
    //system("PAUSE");
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
