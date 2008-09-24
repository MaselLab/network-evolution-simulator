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
  CellState state;
  Genotype indivs[PopSize];
  TimeCourse *timecoursestart[NGENES]; /* array of pointers to list starts */
  TimeCourse *timecourselast[NGENES];
  TimeCourse *start;
  float initmRNA[NGENES], initProteinConc[NGENES], x, kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int hold_genotype_constant = 0;
  int curr_seed;

  verbose = 0;
  initialize_growth_rate_parameters();

  /* parse command-line options */
  while ((c = getopt (argc, argv, "hvgd:r:p:t:c:")) != -1) {
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
        fprintf(stderr, "%s [-d DIRECTORY] [-r DUMMYRUN] [-h] [-g] [-p PLOIDY] [-t DEVELOPMENTTIME] [-c CRITICALSIZE]\n", argv[0]);
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

  /* now create the population of cells */
  for (j = 0; j < PopSize; j++) {
    if (j==PopSize-1) output=1;
    initialize_genotype(&indivs[j], kdis);
    /* if genotype is held constant, start varying the seed *after*
       initialize_genotype, so we can run the same genotype with
       slightly different amounts of noise  */
    if (hold_genotype_constant)
      for (curr_seed=0; curr_seed<dummyrun; curr_seed++) 
         ran1(&seed);
   
    initialize_cell(&state, indivs[j].copies, indivs[j].mRNAdecay, initmRNA, initProteinConc);

    /* print binding sites */
    print_all_binding_sites(indivs[j].copies, indivs[j].allBindingSites, indivs[j].bindSiteCount, 
                            indivs[j].transcriptionFactorSeq, indivs[j].cisRegSeq); 

    develop(&indivs[j], &state, (float) 293.0, timecoursestart, timecourselast);
    fprintf(fperrors,"indiv %d\n",j);
    for (i=0; i < NGENES; i++) {
      if ((output) && j==PopSize-1) print_time_course(timecoursestart[i], i);
      if (verbose) fprintf(fperrors, "deleting gene %d\n", i);
      delete_time_course(timecoursestart[i]);
      timecoursestart[i] = timecourselast[i] = NULL;
    }
    free_mem_CellState(&state);
  }

  for (j = 0; j < PopSize; j++) {
    free(indivs[j].allBindingSites);
  }
  fclose(fperrors);
  fclose(fp_cellsize);
  fclose(fp_growthrate);
  fclose(fp_tfsbound);
}
