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

typedef int (*compfn)(const void*, const void*);
//typedef float (*compfnf)(const void*, const void*);

int intcmp(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

float fltcmp(const void *a, const void *b)
{
      return *(float *)a - *(float *)b;
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

float compareWeights(float *elem1, float *elem2){
      if(elem1 < elem2)
          return -1;
      else if(elem1 > elem2)
          return 1;
      else 
          return 0;
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
  printf("\n");
  
  /* new API requires a Genotype clone */
  Genotype *UNUSED_clone = NULL;
  /* create sequences and binding site matrix */
  initialize_genotype(&indiv, UNUSED_clone, kdis, 0);
  
  /*Sort binding sites from smallest left_edge_position to largest left_edge_position*/
  qsort((void *) &(indiv.all_binding_sites[0]), indiv.sites_per_gene[0],                                 
            sizeof(struct AllTFBindingSites),(compfn)compare );
            
  
          
  float Koff[5];
  float Gibbs;
  float RTlnKr;
  float temperature = 293.0;
  RTlnKr = GASCONSTANT * temperature * log(KR);
  
  int startSite, TFBS, startNum;
  int lem, bob;
  int *left_edge_pos, *startPos, *hammDist, *TFon; 
  
  
  startNum = indiv.all_binding_sites[0].left_edge_pos;
  printf("startNum = %d\n", startNum);
  
  TFBS =0;
  while(indiv.all_binding_sites[TFBS].left_edge_pos < startNum + HIND_LENGTH){
     TFBS++;
  }
  printf("TFBS=%d\n", TFBS);
  
  
  float Kon[TFBS];
  float *weight;
  
  left_edge_pos = malloc(indiv.sites_per_gene[0]*sizeof(int));
  startPos=malloc(TFBS*sizeof(int));
  hammDist = malloc(TFBS *sizeof(int));
  TFon = malloc(TFBS*sizeof(int));
  weight = malloc(TFBS*sizeof(float));
  
  int jo, jojo;
  for( jo=0; jo<3; jo++){
    Gibbs = (((float) jo)/3.0 - 1.0) * RTlnKr;
    //koffCheck = NUMSITESINGENOME*kon*0.25/exp(-Gibbs/(GASCONSTANT*temperature));
    Koff[jo] = -Gibbs;
  }
  
  for(jojo=0; jojo<3; jojo++){
    printf("Koff[%d] = %f\n", jojo, Koff[jojo]);
  }
            
  startSite = 0;
  
  system("PAUSE");
  
  for (lem =0; lem<TFBS; lem++) {
      startPos[lem] = indiv.all_binding_sites[lem+startSite].left_edge_pos;
      hammDist[lem] = indiv.all_binding_sites[lem+startSite].hamming_dist;
      TFon[lem] = indiv.all_binding_sites[lem+startSite].tf_id;
      bob = indiv.all_binding_sites[lem+startSite].tf_id;
      Kon[lem] = initProteinConc[bob];
      weight[lem] = initProteinConc[bob] * Koff[(indiv.all_binding_sites[lem+startSite].hamming_dist)];
      
      printf("%d", indiv.all_binding_sites[lem+startSite].left_edge_pos);
      printf(" Hd = %d   tf = %d", hammDist[lem], TFon[lem]);
      printf(" Kon = %f", Kon[lem]);
      printf(" weight = %.2f\n", weight[lem]);    
    }
    system("PAUSE");
    printf("\n");
    qsort((void *) &(weight[0]), TFBS, sizeof(float*), (compfn)compareWeights);
    for(lem=0; lem<TFBS; lem++){
        printf(" weight = %.2f\n", weight[lem]);    
    }
     /* TO DO: Sort possible TFs according to weight. Decide on suitable 'unbound weight.' 
               Sum all weights and find prob of each TF. Pick random num to decide.*/          
    
  printf("Here\n");
  
  
  system("PAUSE");	
  return 0;
}
