/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- 
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
  FILE *fpout, *fpkdis;
  char fperrors_name[80];
  char fp_cellsize_name[80];
  char fp_growthrate_name[80];
  char fp_tfsbound_name[80];
  int i, j, k, gen;
  CellState state[POP_SIZE];
  Genotype indivs[POP_SIZE];
  TimeCourse *timecoursestart[POP_SIZE][NGENES]; /* array of pointers to list starts */
  TimeCourse *timecourselast[POP_SIZE][NGENES];
  float initmRNA[NGENES], initProteinConc[NGENES], x, kdis[NUM_K_DISASSEMBLY];

  int directory_success;
  int hold_genotype_constant = 0;
  int output_binding_sites = 0;
  int curr_seed;

  int c;  /* for getopt_long() */

  verbose = 0;
  initialize_growth_rate_parameters();

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
        {"criticalsize",  required_argument, 0, 'c'},
        {"genotypeconst",  no_argument, 0, 'g'},
        {"kon",  required_argument, 0, 0},
        {"verbose", no_argument,  0, 'v'},
        {"help",  no_argument, 0, 'h'},
        {"outputbindingsites",  no_argument, 0, 'o'},
        {0, 0, 0, 0}
      };
    
    /* `getopt_long' stores the option index here. */
    int option_index = 0;
    
    c = getopt_long (argc, argv, "d:r:p:t:c:gvho",
                     long_options, &option_index);
    
    /* Detect the end of the options. */
    if (c == -1)
      break;
    
    switch (c)  {
    case 0:  /* long option without a short arg */
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      if( strcmp("kon", long_options[option_index].name) == 0) {
        char *endptr;
        kon  = strtof(optarg, &endptr);
        printf("setting kon=%f\n", kon);
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
      fprintf(stderr, "Usage: %s [OPTION]\n\
\n\
 -g,  --genotypeconst       hold genotype constant\n\
 -o,  --outputbindingsites  print out binding sites\n\
 -d,  --directory=DIRECTORY directory to store output\n\
 -r,  --randomseed=SEED     random seed\n\
 -p,  --ploidy=PLOIDY       ploidy (1=haploid, 2=diploid)\n\
 -t,  --timedev=TIME        length of time to run development\n\
 -c,  --criticalsize=SIZE   critical size for cell division\n\
      --kon=KON             kon value\n\
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

  if (directory_success==-1) 
    if (errno == EEXIST) {
      fprintf(stderr, "directory '%s' already exists\n", output_directory);
    } else {
      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
      exit(-1);
    }

  /* create error output file */
  sprintf(fperrors_name, "%s/netsimerrors.txt", output_directory);
  fperrors = fopen(fperrors_name, "w");

  /* create output files for cell size */
  sprintf(fp_cellsize_name, "%s/cellsize.dat", output_directory);
  if ((fp_cellsize = fopen(fp_cellsize_name,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", fp_cellsize_name);

  /* create output files for growth rate */
  sprintf(fp_growthrate_name, "%s/growthrate.dat", output_directory);
  if ((fp_growthrate = fopen(fp_growthrate_name,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", fp_growthrate_name);

  sprintf(fp_tfsbound_name, "%s/tfsbound.dat", output_directory);
  if ((fp_tfsbound = fopen(fp_tfsbound_name,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", fp_tfsbound_name);

  /* slight hack to initialize seed  */
  if (!hold_genotype_constant)
    for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);

  /* initialize protein concentrations */
  for (i=0; i<NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }

  /* get the kdis.txt values */
  fpkdis = fopen("kdis.txt","r");
  for (j = 0; j < NUM_K_DISASSEMBLY; j++) {
    fscanf(fpkdis,"%f", &kdis[j]);
  }
  fclose(fpkdis);

  /* now create and run the population of cells */
  develop(indivs, state, timecoursestart, timecourselast, (float) 293.0, initmRNA, initProteinConc, kdis, hold_genotype_constant, output_binding_sites);

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
  fclose(fp_cellsize);
  fclose(fp_growthrate);
  fclose(fp_tfsbound);
}
