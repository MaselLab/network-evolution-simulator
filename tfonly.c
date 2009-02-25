/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Standalone TF-only simulator
 * Authors: Alex Lancaster
 * Copyright (c) 2009 Arizona Board of Regents (University of Arizona)
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

#include "random.h"
#include "lib.h"
#include "netsim.h"

int main(int argc, char *argv[])
{
  int i, j;

  float x, t, dt;
  GillespieRates rates;
  Genotype genotype;
  CellState state;
  KonStates konStates;
  float *koffvalues;
  float konrate;

  // FIXME: ultimately won't be used, declare here for moment 
  float transport[NGENES];  /* transport rates of each mRNA */
  float mRNAdecay[NGENES];  /* mRNA decay rates */

  int maxbound2, maxbound3;

  TimeCourse *timecoursestart[NPROTEINS]; /* array of pointers to list starts */
  TimeCourse *timecourselast [NPROTEINS];
  float initmRNA[NGENES], initProteinConc[NPROTEINS], kdis[NUM_K_DISASSEMBLY];

  int hold_genotype_constant = 0;
  int output_binding_sites = 0;
  // TODO, FIXME: problem with commenting these out
  // simulation goes haywire, probably an uninitialized variable issue
  // will look at later - AKL
  int no_fixed_dev_time = 0; /* switch on/off fixed development time  */
  int max_divisions = 0;
  int curr_seed;

  int c;  /* for getopt_long() */

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
        {"timedev",  required_argument, 0, 't'},
        {"genotypeconst",  no_argument, 0, 'g'},
        {"recompute-koff",  no_argument, 0, 0},
        {"burnin",  no_argument, 0, 'b'},
        {"kon",  required_argument, 0, 0},
        {"konafter",  required_argument, 0, 0},
        {"verbose", no_argument,  0, 'v'},
        {"help",  no_argument, 0, 'h'},
        {"outputbindingsites",  no_argument, 0, 'o'},
        {0, 0, 0, 0}
      };
    
    /* `getopt_long' stores the option index here. */
    int option_index = 0;
    char *endptr;
    
    c = getopt_long (argc, argv, "d:r:t:s:gvhob",
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
      } else if (strcmp("recompute-koff", long_options[option_index].name) == 0) {
        printf("recompute rates->koff\n");
        recompute_koff = 1;
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
    case 't':
      tdevelopment = atof(optarg);
      break;
    case 'g':
      hold_genotype_constant = 1;
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
 -g,  --genotypeconst       hold genotype constant\n\
 -o,  --outputbindingsites  print out binding sites\n\
 -d,  --directory=DIRECTORY directory to store output\n\
 -r,  --randomseed=SEED     random seed\n\
 -t,  --timedev=TIME        length of time to run development\n\
 -b,  --burnin              whether to do burn-in (off by default)\n\
      --kon=KON             initial kon value\n\
      --konafter=KON        kon value post-burnin\n\
      --recompute-koff      recompute koff every time step\n\
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
  
  /* create output directory */
  create_output_directory(output_directory);
  /* create error output file */
  create_output_file("netsimerrors.txt", output_directory, &(fperrors), -1);
  /* create output files for koff values and TF occupancies */
  create_output_file("koff", output_directory, &(fp_koff[0]), 0);
  create_output_file("tfsbound", output_directory, &(fp_tfsbound[0]), 0);
  
  /* slight hack to initialize seed  */
  if (!hold_genotype_constant)
    for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);

  /* initialize protein and mRNA (not used) concentrations */
  for (i=0; i < NPROTEINS; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }

  /* get the kdis.txt values */
  read_kdisassembly(kdis);

  maxbound2 = maxbound;
  maxbound3 = 10*maxbound;

  /* set to be haploid */
  current_ploidy = 1;  

  /* initialize genotype (1 cis-reg, 10 TFs) */
  initialize_genotype(&genotype, kdis);

  //printf("cisregseq=%s\n", (char *)genotype.cisRegSeq);

  /* initialize cell state */
  initialize_cell(&state, 0, genotype.copies, genotype.mRNAdecay, initmRNA, initProteinConc, burn_in);

  // FIXME: hack
  /* reset mRNA decay rates */
  for (i=0; i < NGENES; i++) 
    genotype.mRNAdecay[i] = 0;

  if (output_binding_sites)  /* if command-line option was passed in */
    print_all_binding_sites(genotype.copies, genotype.allBindingSites, genotype.bindSiteCount, 
                            genotype.transcriptionFactorSeq, genotype.cisRegSeq, genotype.tfsStart); 

  /* set cell temperature and value of RTlnKr constant */
  state.temperature = 293.0;
  state.RTlnKr = GasConstant * state.temperature * log(Kr);

  /* initialize time courses */
  for (i=0; i < NPROTEINS; i++) {
    timecoursestart[i] = NULL;
    timecourselast[i] = NULL;
  } 

  add_time_points((float) 0.0, state.proteinConc,  timecoursestart,  timecourselast);

  /* initialize cached data structure: allocate memory */

  initialize_cell_cache(&state, genotype, &konStates, &koffvalues, maxbound2, maxbound3);  
  rates.total = 0.0; //LOCAL HACK FIXME

  printf("state.tfBoundCount=%d, state.tfHinderedCount=%d\n", state.tfBoundCount, state.tfHinderedCount);

  /* initialize transcriptional state of genotype */
  calc_from_state(&genotype, &state, &rates, &konStates, transport, mRNAdecay);


  printf("rates.koff=%g, rates.salphc=%g, rates.maxSalphc=%g, rates.minSalphc=%g, rates.total=%g\n", 
         rates.koff, rates.salphc, rates.maxSalphc, rates.minSalphc, rates.total);

  while (t < tdevelopment) {       /* run until tdevelopment reached */

    /* compute total konrate (which is constant over the Gillespie step) */
    if (konStates.nkon==0) konrate = (-rates.salphc);
    else calc_kon_rate(dt, &konStates, &konrate); 

    x = expdev(&seed);        /* draw random number */
    
    /* do first Gillespie step to chose next event */
    calc_dt(&x, &dt, &rates, &konStates, mRNAdecay, genotype.mRNAdecay,
            state.mRNACytoCount, state.mRNATranslCytoCount);


    printf("t=%g, dt=%g state.tfBoundCount=%d, state.tfHinderedCount=%d\n", 
           t, dt, state.tfBoundCount, state.tfHinderedCount);
    printf("rates.koff=%g, rates.salphc=%g, rates.maxSalphc=%g, rates.minSalphc=%g, rates.total=%g\n", 
           rates.koff, rates.salphc, rates.maxSalphc, rates.minSalphc, rates.total);

    print_tf_occupancy(&state, genotype.allBindingSites, t);

    log_snapshot(&rates,
                 &state,
                 &genotype,
                 &konStates,
                 &koffvalues,
                 mRNAdecay,
                 transport, 
                 konrate,
                 x,
                 t);
    

    for (i=0; i < TFGENES; i++)
      printf("protein[%d]=%g\n", i, state.proteinConc[i]);

    if (t+dt < tdevelopment) {

      /* 
       * STOCHASTIC EVENT: TF binding event
       */
      if (x < rates.salphc + (konrate)) {   /* add variable (salphc) and constant (konrate) */
        tf_binding_event(&rates, &state, &genotype, &konStates, koffvalues,
                         timecoursestart, timecourselast, (konrate), dt, t, 
                         maxbound2, maxbound3, 0);
      } else {
        x -= (rates.salphc + (konrate));        
        /* 
         * STOCHASTIC EVENT: a TF unbinds (koff) 
         */
        if (x < rates.koff) {  
          tf_unbinding_event(&rates, &state, &genotype, &konStates, koffvalues,
                             timecoursestart,  timecourselast, (konrate), dt, t, x, 0);
        } else {
          //x -= rates.koff;  
          /*
           * FALLBACK: shouldn't get here, previous
           * events should be exhaustive
           */
          printf("ERROR: shouldn't get here!\n");
        }
      }
      
      /* Gillespie step: advance time to next event at dt */
      t += dt;
      LOG_VERBOSE_NOCELLID("dt=%g t=%g\n", dt, t);
      printf("dt=%g t=%g\n", dt, t);
    } else {
      /* we will reach the end of development in dt */
      LOG_VERBOSE_NOCELLID("finish at t=%g dt=%g\n", t, dt);
      
      /* do remaining dt */
      dt = tdevelopment - t;
      
      /* advance to end of development (this exits the outer while loop) */
      t = tdevelopment;
    }
  }

  for (i=0; i < NPROTEINS; i++) {
    //if ((output) && j==POP_SIZE-1) print_time_course(timecoursestart[j][i], i, j);
    if ((output)) print_time_course(timecoursestart[i], i, 0);
    if (verbose) fprintf(fperrors, "deleting gene %d\n", i);
    delete_time_course(timecoursestart[i]);
    timecoursestart[i] = timecourselast[i] = NULL;
  }
    
  fprintf(fperrors,"indiv %d\n", j);
  free_mem_CellState(&state);
  free(genotype.allBindingSites);

  fclose(fperrors);
  fclose(fp_koff[0]);
  fclose(fp_tfsbound[0]);
}
