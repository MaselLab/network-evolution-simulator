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

#include "random.h"
#include "lib.h"
#include "netsim.h"



int main(int argc, char *argv[])
{
  int i, j;
  CellState state[N_para_threads];
  Genotype indivs[N_para_threads];
//  TimeCourse *timecoursestart[2][NPROTEINS]; /* array of pointers to list starts */
//  TimeCourse *timecourselast[2][NPROTEINS];
  float kdis[NUM_K_DISASSEMBLY];

  int output_binding_sites = 0; /*verbose flag*/
  int no_fixed_dev_time = 0; /* 0 = fixed development time, 1 = divides when ready  */
  int curr_seed; /*still needs to be fixed by Barry*/

  /* create output directory if needed */
  create_output_directory(output_directory);

  /* create error output file */
  create_output_file("netsimerrors.txt", output_directory, &(fperrors), -1);

  /* create output files for cell size, growth rate, TFs */
  for (j = 0; j < 2; j++) {

    create_output_file("cellsize", output_directory, &(fp_cellsize[j]), j);
#if 0 /* currently disable these file outputs */
    create_output_file("growthrate", output_directory, &(fp_growthrate[j]), j);
    create_output_file("tfsbound", output_directory, &(fp_tfsbound[j]), j);
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

  init_run_pop(indivs, state, 
//  timecoursestart, timecourselast, 
				(float) 293.0,
              kdis, output_binding_sites,no_fixed_dev_time);
//			  , no_fixed_dev_time, max_divisions);

//  print_all_protein_time_courses(timecoursestart, timecourselast);
  //system("PAUSE");
  /* cleanup memory */
  
  for (j = 0; j < 1; j++) {
    fprintf(fperrors,"cleanup cell %03d\n", j);
    
//    for (i=0; i < NPROTEINS; i++) {
//		
//	
//      //TODO: FIX LOGGING
//      //LOG_VERBOSE("deleting protein %02d timecourse\n", i);
//      delete_time_course(timecoursestart[j][i]); 
//           
//      timecoursestart[j][i] = NULL;
//      
//	  timecourselast[j][i] = NULL;
//	  
//    }
   
    free_mem_CellState(&state[j]);
    
    free(indivs[j].all_binding_sites); // memo leakage here. but should not matter when the code is complete
  }

  /* close file descriptors */
  fclose(fperrors);
  
  for (j = 0; j < 2; j++) {
    fclose(fp_cellsize[j]);



#if 0  /* currently disable */
    fclose(fp_growthrate[j]);
    fclose(fp_tfsbound[j]);
    fclose(fp_rounding[j]);
#endif
  }

  return 0;
}
