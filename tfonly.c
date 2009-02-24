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
  FILE *fpkdis;
  char fperrors_name[80];
  char fp_cellsize_name[80];
  char fp_koff_name[80];
#if 0
  char fp_growthrate_name[80];
#endif
  char fp_tfsbound_name[80];
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

  int directory_success;
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

  /* create error output file */
  sprintf(fperrors_name, "%s/netsimerrors.txt", output_directory);
  fperrors = fopen(fperrors_name, "w");

  /* create output files for cell size, growth rate, TFs */
  for (j = 0; j < POP_SIZE; j++) {
    sprintf(fp_cellsize_name, "%s/cellsize-%03d.dat", output_directory, j);
    if ((fp_cellsize[j] = fopen(fp_cellsize_name,"w"))==NULL)
      fprintf(fperrors,"error: Can't open %s file\n", fp_cellsize_name);

    sprintf(fp_koff_name, "%s/koff-%03d.dat", output_directory, j);
    if ((fp_koff[j] = fopen(fp_koff_name,"w"))==NULL)
      fprintf(fperrors,"error: Can't open %s file\n", fp_koff_name);


    // TODO: currently disable
#if 0
    sprintf(fp_growthrate_name, "%s/growthrate-%03d.dat", output_directory, j);
    if ((fp_growthrate[j] = fopen(fp_growthrate_name,"w"))==NULL)
      fprintf(fperrors,"error: Can't open %s file\n", fp_growthrate_name);
#endif
    sprintf(fp_tfsbound_name, "%s/tfsbound-%03d.dat", output_directory, j);
    if ((fp_tfsbound[j] = fopen(fp_tfsbound_name,"w"))==NULL)
      fprintf(fperrors,"error: Can't open %s file\n", fp_tfsbound_name);
  }
  
  /* slight hack to initialize seed  */
  if (!hold_genotype_constant)
    for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);

  initialize_growth_rate_parameters();

  /* initialize protein and mRNA (not used) concentrations */
  for (i=0; i < NPROTEINS; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }

  /* get the kdis.txt values */
  if ((fpkdis = fopen("kdis.txt","r"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", "kdis.txt");
  for (j = 0; j < NUM_K_DISASSEMBLY; j++) {
    fscanf(fpkdis, "%f", &kdis[j]);
  }
  fclose(fpkdis);

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

  /* output binding sites */
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

  printf("state.tfBoundCount=%d, state.tfHinderedCount=%d\n", state.tfBoundCount, state.tfHinderedCount);

  /* initialize transcriptional state of genotype */
  calc_from_state(&genotype, &state, &rates, &konStates, transport, mRNAdecay);

  printf("rates.koff=%g, rates.salphc=%g, rates.maxSalphc=%g, rates.minSalphc=%g, rates.total=%g\n", 
         rates.koff, rates.salphc, rates.maxSalphc, rates.minSalphc, rates.total);

  while (t < tdevelopment) {       /* run until tdevelopment reached */

    /* compute total konrate (which is constant over the Gillespie step) */
    if (konStates.nkon==0) konrate = (-rates.salphc);
    else calc_kon_rate(dt, &konStates, &konrate); 
    
    /* do first Gillespie step to chose next event */
    calc_dt(&x, &dt, &rates, &konStates, mRNAdecay, genotype.mRNAdecay,
            state.mRNACytoCount, state.mRNATranslCytoCount);

    x = expdev(&seed);        /* draw random number */

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
  fclose(fp_cellsize[0]);
  fclose(fp_koff[0]);
  // TODO: currently disable
#if 0
  fclose(fp_growthrate[0]);
#endif
  fclose(fp_tfsbound[0]);
}
