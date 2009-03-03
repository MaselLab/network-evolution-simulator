/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008 Arizona Board of Regents (University of Arizona)
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
  CellState state[POP_SIZE];
  Genotype indivs[POP_SIZE];
  TimeCourse *timecoursestart[POP_SIZE][NGENES]; /* array of pointers to list starts */
  TimeCourse *timecourselast[POP_SIZE][NGENES];
  float initmRNA[NGENES], initProteinConc[NGENES], kdis[NUM_K_DISASSEMBLY];

  int hold_genotype_constant = 0;
  int output_binding_sites = 0;
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
        {"ploidy",  required_argument, 0, 'p'},
        {"timedev",  required_argument, 0, 't'},
        {"timemax",  required_argument, 0, 0},
        {"criticalsize",  required_argument, 0, 'c'},
        {"divisions",  required_argument, 0, 's'},
        {"genotypeconst",  no_argument, 0, 'g'},
        {"random-replication",  no_argument, 0, 0},
        {"recompute-koff",  no_argument, 0, 0},
        {"recompute-kon",  no_argument, 0, 0},
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
    
    c = getopt_long (argc, argv, "d:r:p:t:c:s:gvhonb",
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
      } else if (strcmp("random-replication", long_options[option_index].name) == 0) {
        printf("making replication times random in S phase\n");
        random_replication_time = 1;
      } else if (strcmp("recompute-koff", long_options[option_index].name) == 0) {
        printf("recompute rates->koff\n");
        recompute_koff = 1;
      } else if (strcmp("recompute-kon", long_options[option_index].name) == 0) {
        printf("recompute rates->kon\n");
        recompute_kon = 1;
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
      --random-replication  make replication times in S phase random\n\
      --timesphase=TIME     length of S-phase (30 mins by default)\n\
      --timeg2phase=TIME    length of G2-phase (30 mins by default)\n\
      --growthscaling=GS    amount to accelerate the growth rate\n\
                              (2.0 by default)\n\
      --recompute-koff      recompute koff rates after a fixed number of operations\n\
      --recompute-kon       recompute kon rates after a fixed number of operations\n\
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
    create_output_file("koff", output_directory, &(fp_koff[j]), j);
#if 0 // TODO: currently disable
    create_output_file("growthrate", output_directory, &(fp_growthrate[j]), j);
#endif
    create_output_file("tfsbound", output_directory, &(fp_tfsbound[j]), j);
    create_output_file("rounding", output_directory, &(fp_rounding[j]), j);
    // print header
    fprintf(fp_rounding[j], "t koff transport mRNAdecay picDisassembly salphc maxSalphc minSalphc \
          acetylationCount deacetylationCount picAssemblyCount \
          transcriptInitCount picDisassemblyCount\n");
  }
  
  /* slight hack to initialize seed  */
  if (!hold_genotype_constant)
    for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);

  initialize_growth_rate_parameters();

  /* initialize protein concentrations */
  for (i=0; i < NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }

  /* get the kdis.txt values */
  read_kdisassembly(kdis);

  /* now create and run the population of cells */
  develop(indivs, state, timecoursestart, timecourselast, (float) 293.0, initmRNA, initProteinConc, 
          kdis, hold_genotype_constant, output_binding_sites, no_fixed_dev_time, max_divisions);

  for (j = 0; j < POP_SIZE; j++) {
    fprintf(fperrors,"indiv %d\n", j);
    for (i=0; i < NGENES; i++) {
      //if ((output) && j==POP_SIZE-1) print_time_course(timecoursestart[j][i], i, j);
      if ((output)) print_time_course(timecoursestart[j][i], i, j);
      if (verbose) fprintf(fperrors, "deleting gene %d\n", i);
      delete_time_course(timecoursestart[j][i]);
      timecoursestart[j][i] = timecourselast[j][i] = NULL;
    }
    free_mem_CellState(&state[j]);
    free(indivs[j].allBindingSites);
  }
  fclose(fperrors);
  for (j = 0; j < POP_SIZE; j++) {
    fclose(fp_cellsize[j]);
    fclose(fp_koff[j]);
    // TODO: currently disable
#if 0
    fclose(fp_growthrate[j]);
#endif
    fclose(fp_tfsbound[j]);
    fclose(fp_rounding[j]);
  }
}
