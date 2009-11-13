/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
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

#include "randomJ.h"
#include "libJ.h"
#include "netsimJ.h"


//#define NITER 3000

//FILE *recordFile;




/*Reorders binding sites from smallest left_edge_position to largest left_edge_position 
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
      for(n=0; n<indiv.sites_per_gene[gene]; n++){
               fprintf(recordFile, "%d   %d\n",n,  indiv.all_binding_sites[start+n].left_edge_pos);
      }
      
      
      struct Dtype *arrayD;
 arrayD = malloc(NITER*sizeof(struct Dtype));
 
 //each struct holds ID, BS numer, start positon, hamming distance, concentration and weight
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
               //FIX THIS! SUM IS NOT CORRECT! ROUNDING ERRORS!    
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
      if(arrayWT[b].tfIDon == 11){//do nothing, there is nothing bound 
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
      for(lem=0; lem<size; lem++){
      if(arrayD[lem].repress < arrayD[lem].ratio){
         activeCount += arrayD[lem].count;} 
      }
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
           printf("geneSum = %d, i = %d, WHAT = %f\n\n", geneSum, i, gene_active[i] );
           //fprintf(outputFile, "geneSum = %d, i = %d, WHAT = %f\n\n", geneSum, i, gene_active[i] );
           geneSum += indiv.sites_per_gene[i];
        
  }
 }  */


int main(int argc, char *argv[])
{
  int i, j;
  CellState state[POP_SIZE];
  Genotype indivs[POP_SIZE];
  TimeCourse *timecoursestart[POP_SIZE][NPROTEINS]; /* array of pointers to list starts */
  TimeCourse *timecourselast[POP_SIZE][NPROTEINS];
  float kdis[NUM_K_DISASSEMBLY];

  int output_binding_sites = 0; /*verbose flag*/
  int no_fixed_dev_time = 0; /* on = fixed development time, off = divides when ready  */
  int max_divisions = 0; /* per population */
  int curr_seed; /*still needs to be fixed by Barry*/

  int c;  /* argv index */

  burn_in = 0;  /* don't do burn-in of kon by default */
  verbose = 0;  /* no verbose out by default */

  /* parse command-line options */
  while (1) {
    static struct option long_options[] =
      {
        /* These options set a flag. */
        /* These options don't set a flag.
           We distinguish them by their indices. */
        {"directory",  required_argument, 0, 'd'},
        {"randomseed", required_argument, 0, 'r'},
        {"ploidy",  required_argument, 0, 'p'},
        {"timedev",  required_argument, 0, 't'},
        {"timemax",  required_argument, 0, 0},
        {"criticalsize",  required_argument, 0, 'c'},
        {"divisions",  required_argument, 0, 's'},
        {"random-replication",  no_argument, 0, 0},
        {"no-recompute-koff",  no_argument, 0, 0},
        {"no-recompute-kon",  no_argument, 0, 0},
        {"nofixedtime",  no_argument, 0, 'n'},
        {"burnin",  no_argument, 0, 'b'},
        {"kon",  required_argument, 0, 0},
        {"konafter",  required_argument, 0, 0},
        {"timesphase",  required_argument, 0, 0},
        {"timeg2phase",  required_argument, 0, 0},
        {"growthscaling",  required_argument, 0, 0},
        {"verbose", no_argument,  0, 'v'},
        {"help",  no_argument, 0, 'h'},
        {"outputbindingsites",  no_argument, 0, 'o'},
        {0, 0, 0, 0}
      };
    
    /* `getopt_long' stores the option index here. */
    int option_index = 0;
    char *endptr;
    
    c = getopt_long (argc, argv, "d:r:p:t:c:s:vhonb",
                     long_options, &option_index);
    
    /* Detect the end of the options. */
    if (c == -1)
      break;
    
    switch (c)  {
    case 0:  /* long option without a short arg */
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      if (strcmp("kon", long_options[option_index].name) == 0) {
        kon  = strtof(optarg, &endptr);
        printf("setting kon=%g\n", kon);
      } else if (strcmp("konafter", long_options[option_index].name) == 0) {
        kon_after_burnin  = strtof(optarg, &endptr);
        printf("setting konafter=%g\n", kon_after_burnin);
        burn_in = 1;
      } else if (strcmp("timesphase", long_options[option_index].name) == 0) {
        time_s_phase  = strtof(optarg, &endptr);
        printf("setting time_s_phase=%g\n", time_s_phase);
      } else if (strcmp("timeg2phase", long_options[option_index].name) == 0) {
        time_g2_phase  = strtof(optarg, &endptr);
        printf("setting time_g2_phase=%g\n", time_g2_phase);
      } else if (strcmp("growthscaling", long_options[option_index].name) == 0) {
        growth_rate_scaling  = strtof(optarg, &endptr);
        printf("setting growthscaling=%g\n", growth_rate_scaling);
      } else if (strcmp("timemax", long_options[option_index].name) == 0) {
        timemax  = strtof(optarg, &endptr);
        printf("setting timemax=%g\n", timemax);
      } else if (strcmp("no-random-replication", long_options[option_index].name) == 0) {
        printf("don't make replication times random in S phase\n");
        random_replication_time = 0;
      } else if (strcmp("no-recompute-koff", long_options[option_index].name) == 0) {
        printf("don't recompute rates->koff\n");
        recompute_koff = 0;
      } else if (strcmp("no-recompute-kon", long_options[option_index].name) == 0) {
        printf("don't recompute rates->kon\n");
        recompute_kon = 0;
      } else { 
        printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
        printf ("\n");
      }
      break;
    case 'd':
      output_directory = optarg;
      break;
    case 'r':
      dummyrun = atoi(optarg);
      break;
    case 'p':
      current_ploidy = atoi(optarg);
      break;
    case 's':
      max_divisions = atoi(optarg);
      break;
    case 't':
      tdevelopment = atof(optarg);
      break;
    case 'n':
      no_fixed_dev_time = 1;
      break;
    case 'c':
      critical_size = atof(optarg);
      break;
    case 'b':
      burn_in = 1;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
      fprintf(stderr, "Usage: %s [OPTION]\n\
\n\
 -o,  --outputbindingsites  print out binding sites\n\
 -d,  --directory=DIRECTORY directory to store output\n\
 -r,  --randomseed=SEED     random seed\n\
 -p,  --ploidy=PLOIDY       ploidy (1=haploid, 2=diploid)\n\
 -t,  --timedev=TIME        length of time to run development\n\
 -t,  --timemax=TIME        set a maximum length of time to run development\n\
                             (no upper limit by default)\n\
 -n,  --nofixedtime         no fixed development time\n\
 -s,  --divisions=DIVISONS  maximum number of divisions\n\
 -c,  --criticalsize=SIZE   critical size for cell division\n\
 -b,  --burnin              whether to do burn-in (off by default)\n\
      --kon=KON             initial kon value\n\
      --konafter=KON        kon value post-burnin\n\
      --no-random-replication  don't make replication times in S phase random\n\
      --timesphase=TIME     length of S-phase (30 mins by default)\n\
      --timeg2phase=TIME    length of G2-phase (30 mins by default)\n\
      --growthscaling=GS    amount to accelerate the growth rate\n\
                              (2.0 by default)\n\
      --no-recompute-koff   don't recompute koff rates after a fixed number of operations\n\
      --no-recompute-kon    don't recompute kon rates after a fixed number of operations\n\
 -h,  --help                display this help and exit\n\
 -v,  --verbose             verbose output to error file\n\
\n", argv[0]);
      exit(0);
      break;
    case 'o':
      output_binding_sites = 1;
      break;
    default:
      abort();
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc) {
    printf ("non-option ARGV-elements: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
    putchar ('\n');
  }

  /* create output directory if needed */
  create_output_directory(output_directory);

  /* create error output file */
  create_output_file("netsimerrors.txt", output_directory, &(fperrors), -1);

  /* create output files for cell size, growth rate, TFs */
  for (j = 0; j < POP_SIZE; j++) {
    create_output_file("cellsize", output_directory, &(fp_cellsize[j]), j);
#if 0 /* currently disable these file outputs */
    //create_output_file("koff", output_directory, &(fp_koff[j]), j);
    create_output_file("growthrate", output_directory, &(fp_growthrate[j]), j);
    //create_output_file("tfsbound", output_directory, &(fp_tfsbound[j]), j);
    create_output_file("rounding", output_directory, &(fp_rounding[j]), j);
    /* print header */
    fprintf(fp_rounding[j], "t koff transport mRNAdecay picDisassembly salphc maxSalphc minSalphc \
          acetylationCount deacetylationCount picAssemblyCount \
          transcriptInitCount picDisassemblyCount\n");
#endif

  }
  
  /* slight hack to initialize seed  */
  for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);

  initialize_growth_rate_parameters();

  /* get the kdis.txt values */
  read_kdisassembly(kdis);

  /* now create and run the population of cells */
  init_run_pop(indivs, state, timecoursestart, timecourselast, (float) 293.0, 
               kdis, output_binding_sites, no_fixed_dev_time, max_divisions);

  print_all_protein_time_courses(timecoursestart, timecourselast);

  system("PAUSE");
  /* cleanup memory */
  for (j = 0; j < POP_SIZE; j++) {
    fprintf(fperrors,"cleanup cell %03d\n", j);
    for (i=0; i < NGENES; i++) {
      LOG_VERBOSE("deleting protein %02d timecourse\n", i);
      delete_time_course(timecoursestart[j][i]);
      timecoursestart[j][i] = timecourselast[j][i] = NULL;
    }
    //free_mem_CellState(&state[j]);
    free(indivs[j].all_binding_sites);
  }
  //system("PAUSE");
  /* close file descriptors */
  fclose(fperrors);
  for (j = 0; j < POP_SIZE; j++) {
    fclose(fp_cellsize[j]);
    
#if 0  /* currently disable */
   // fclose(fp_koff[j]);
    fclose(fp_growthrate[j]);
    //fclose(fp_tfsbound[j]);
    fclose(fp_rounding[j]);
#endif
  }
}
