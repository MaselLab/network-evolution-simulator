/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Standalone TF-only simulator
 * Authors: Alex Lancaster, Barry Rountree
 * Copyright (c) 2009 Arizona Board of Regents (University of Arizona)
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

#include "random.h"
#include "lib.h"
#include "netsim.h"


// no_argument, required_argument and optional_argument are 
// #defined in <getopt.h>.
static struct option long_options[] =
{
	{"directory",  		required_argument, 	0, 'd'},
	{"randomseed", 		required_argument, 	0, 'r'},
	{"timedev",  		required_argument, 	0, 't'},
	{"genotypeconst",  	no_argument, 		0, 'g'},
	{"recompute-koff",  	no_argument, 		0,  3 },
	{"burnin",  		no_argument, 		0, 'b'},
	{"kon",  		required_argument, 	0,  1 },
	{"konafter",  		required_argument, 	0,  2 },
	{"verbose", 		no_argument,  		0, 'v'},
	{"help",  		no_argument, 		0, 'h'},
	{"outputbindingsites",  no_argument, 		0, 'o'},
	{0, 0, 0, 0}
};

static void print_help(char *arg){
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
	\n", arg);
}
static int hold_genotype_constant = 0;
static int output_binding_sites = 0;
static void parse_options(int argc, char **argv){
	int c;  
	int option_index = 0;
	char *endptr;

	// Declared in netsim.h
	burn_in = 0;  // don't do burn-in of kon by default
	verbose = 0;  // no verbose out by default

	while (1){  
		c = getopt_long (argc, argv, "d:r:t:s:gvhob", long_options, &option_index);
		if(c==-1){
			break;
		}
    
		switch (c)  {
			case 1: //global variable
				kon  = strtof(optarg, &endptr);	
				printf("setting kon=%g\n", kon);
				break;
			case 2: //global variable
				kon_after_burnin  = strtof(optarg, &endptr);
				printf("setting konafter=%g\n", kon_after_burnin);
				burn_in = 1;
				break;
			case 3: //global variable
				printf("recompute rates->koff\n");
				recompute_koff = 1;
				break;
			case 'd': //global variable 
				output_directory = optarg;
				break;
			case 'r': //global variable
				dummyrun = atoi(optarg);
				break;
			case 't': //global variable
				tdevelopment = atof(optarg);
				break;
			case 'g': //static
				hold_genotype_constant = 1;
				break;
			case 'b': //global variable
				burn_in = 1;
				break;
			case 'v': //global variable
				verbose = 1;
				break;
			case 'h':
				print_help(argv[0]);
				exit(0);
				break;
			case 'o': //static
				output_binding_sites = 1;
				break;
			default:
				fprintf(stderr, "%s::%d  Unknown option %c.  Aborting.\n", 
					__FILE__, __LINE__, c);
				abort();
		}
	}

	// Print any remaining command line arguments (not options).
	if (optind < argc) {
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
		printf ("%s ", argv[optind++]);
		putchar ('\n');
  	}
} 

static void system_init(){
	int curr_seed;  // This doesn't do what you think it does.
	create_output_directory(output_directory);
	create_output_file("netsimerrors.txt", output_directory, &(fperrors), -1);
	//create_output_file("koff", output_directory, &(fp_koff[0]), 0);
	//create_output_file("tfsbound", output_directory, &(fp_tfsbound[0]), 0);

	if (!hold_genotype_constant){
		// For now, the random number generator is initialized to
		// the first seed value and any subsequent negative seed.
		// This will be fixed shortly.  (Also need to put in a way
		// of specifying the seed on the command line.
		for (curr_seed=0; curr_seed<dummyrun; curr_seed++){
			ran1(&seed);
		}
	}
}

static float initmRNA[NGENES]; 
static float initProteinConc[NPROTEINS]; 
static float kdis[NUM_K_DISASSEMBLY];
static int maxbound2; 
static int maxbound3;
static Genotype genotype;
static CellState state;
static TimeCourse *timecoursestart[NPROTEINS]; /* array of pointers to list starts */
static TimeCourse *timecourselast [NPROTEINS];
static GillespieRates rates;
static KonStates konStates;
static float *koffvalues;
static float mRNAdecay[NGENES];  /* mRNA decay rates */
static int UNUSED_genotypeID=0;    	// required by initialize_genotype.
static Genotype *UNUSED_clone=NULL; 	// required by initialize_genotype.
static int UNUSED_cellID=0;	      	// required by calc_dt.
static float transport[NGENES];  /* transport rates of each mRNA */

static void model_init(){
	int i, j;

	// initialize mRNA concentrations  (not used) 
  for (i=0; i < NGENES; i++) {
      initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }

	// initialize protein concentrations 
	for (i=0; i < NPROTEINS; i++) {
		initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
	}

	// get the kdis.txt values 
	read_kdisassembly(kdis);

	maxbound2 = MAXBOUND;
	maxbound3 = 10*MAXBOUND;

	// set to be haploid (global)
	current_ploidy = 1;  

	// initialize genotype (1 cis-reg, 10 TFs) 
	initialize_genotype(&genotype, UNUSED_clone, kdis, UNUSED_genotypeID);

	// initialize cell state 
	initialize_cell(&state, 0, genotype.copies, genotype.mRNAdecay, initmRNA, initProteinConc, burn_in);

  //printf("cisregseq=%s\n", (char *)genotype.cisreg_seq);


	// FIXME: hack
	// reset mRNA decay rates 
	for (i=0; i < NGENES; i++) { 
		genotype.mRNAdecay[i] = 0;
	}

	if (output_binding_sites){  // if command-line option was passed in 
		print_all_binding_sites(genotype.copies, genotype.all_binding_sites, genotype.binding_sites_num, 
				    genotype.tf_seq, genotype.cisreg_seq, genotype.site_id_pos); 
	}

	/* set cell temperature and value of RTlnKr constant */
	state.temperature = 293.0;
	state.RTlnKr = GASCONSTANT * state.temperature * log(KR);

	// initialize time courses 
	for (i=0; i < NPROTEINS; i++) {
		timecoursestart[i] = NULL;
		timecourselast[i] = NULL;
	} 

	add_time_points((float) 0.0, state.protein_conc,  timecoursestart,  timecourselast);

	// initialize cached data structure: allocate memory 

	initialize_cell_cache(&state, genotype, &konStates, &koffvalues, maxbound2, maxbound3);  

	printf("state.tf_bound_num=%d, state.tf_hindered_num=%d\n", state.tf_bound_num, state.tf_hindered_num);

  for (i=0; i < NGENES; i++) {
    // FIXME: here we override the initial cell state and ensure there is at least one
    // mRNA in the cytoplasm before we initialize the cell state, otherwise there is no initial
    // kon rates, also ensure that there are no mRNAs in the nucleus so that no transport occurs
    state.mRNA_cyto_num[i] = 1;
    state.mRNA_nuclear_num[i] = 0;
    printf("mRNA_cyto_num[%d]=%d\n", i, state.mRNA_cyto_num[i]);
    printf("mRNA_nuclear_num[%d]=%d\n", i, state.mRNA_nuclear_num[i]);
  }


	// initialize transcriptional state of genotype 
	calc_from_state(&genotype, &state, &rates, &konStates, transport, mRNAdecay);

  // FIXME: hack to work around the fact that calc_from_state() sets
  // the acetylationCount to a non-zero factor, meaning it has a
  // non-zero rate, but we don't model the acetylation process
  // TODO: probably should copy calc_from_state here, and modify it locally
  // and refactor back into netsim.c later
  for (i=0; i < NGENES; i++) {
    for (j=0; j < genotype.copies[i]; j++) {
      rates.acetylation_num[j] = 0;
    }
  }

  for (i=0; i < NGENES; i++) {
    printf("after mRNA_cyto_num[%d]=%d\n", i, state.mRNA_cyto_num[i]);
    printf("after mRNA_nuclear_num[%d]=%d\n", i, state.mRNA_nuclear_num[i]);
  }

  printf("rates->transport=%g\n", rates.transport);
  
  // printf("rates.koff=%g, rates.salphc=%g, rates.max_salphc=%g, rates.min_salphc=%g, rates.subtotal=%g\n", 
  //       rates.koff, rates.salphc, rates.max_salphc, rates.min_salphc, rates.subtotal);
	//printf("in model_init(): rates.koff=%g, rates.salphc=%g, rates.max_salphc=%g, rates.min_salphc=%g, rates.subtotal=%g\n", 
	// rates.koff, rates.salphc, rates.max_salphc, rates.min_salphc, rates.subtotal);
}

void run(){
	int i;
	float x, t, dt;
	float konrate = 0.0;
	// FIXME: ultimately won't be used, declare here for moment 
	// TODO, FIXME: problem with commenting these out
	// simulation goes haywire, probably an uninitialized variable issue
	// will look at later - AKL
	// int no_fixed_dev_time = 0; // switch on/off fixed development time  
	// int max_divisions = 0;

	while (t < tdevelopment) {       // run until tdevelopment reached 

		x = expdev(&seed);        // draw random number

    // log a snapshot to netsimerrors.txt if verbose flag is passed
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

		// do first Gillespie step to chose next event 
		calc_dt(&x, &dt, &rates, &konStates, mRNAdecay, genotype.mRNAdecay,
			state.mRNA_cyto_num, state.mRNA_transl_cyto_num, UNUSED_cellID);

		// compute total konrate (which is constant over the Gillespie step) 
		if (konStates.nkon==0){
			konrate = (-rates.salphc);
		}else{ 
			calc_kon_rate(dt, &konStates, &konrate); 
		}

		printf("t=%g, dt=%g state.tf_bound_num=%d, state.tf_hindered_num=%d\n", t, dt, state.tf_bound_num, state.tf_hindered_num);
		printf("OFF [rates.koff=%g] ON [konrate=%g, rates.salphc=%g, (rates.max_salphc=%g, rates.min_salphc=%g)] rates.subtotal=%g, TOTAL [%g]\n", 
           rates.koff, konrate, rates.salphc, rates.max_salphc, rates.min_salphc, rates.subtotal, rates.subtotal+konrate);
    printf("TRANSPORT rates.transport=%g, rates.mRNAdecay=%g, rates.pic_disassembly=%g, rates.acetylation_num=%d, rates.deacetylation_num=%d, rates.transcript_init_num=%d, rates.pic_assembly_num=%d, rates.pic_disassembly_num=%d\n",
               rates.transport, rates.mRNAdecay, rates.pic_disassembly, rates.acetylation_num[0], rates.deacetylation_num[0], rates.transcript_init_num[0], 
               rates.pic_assembly_num[0], rates.pic_disassembly_num[0]);


    float diff = (rates.subtotal + konrate) - ((rates.koff + rates.salphc) + konrate);
    if (abs(diff > 1e-7)) {
      printf("\nDIFF=TOTAL(%g)- SUM(%g)=%g comprised of OFF [rates.koff=%g] and ON [konrate=%g, rates.salphc=%g]\n\n", 
             rates.subtotal + konrate, (rates.koff + rates.salphc) + konrate, diff, rates.koff, konrate, rates.salphc);
    }
    system("PAUSE");
		//print_tf_occupancy(&state, genotype.all_binding_sites, t);

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


		x = ran1(&seed)*(rates.subtotal + konrate); // Note that rates.subtotal = rates.salphc + rates.koff
    
		if (t+dt < tdevelopment) {

			if (x < rates.salphc + (konrate)) {   // add variable (salphc) and constant (konrate) 
				// STOCHASTIC EVENT: TF binding event
				tf_binding_event(&rates, &state, &genotype, &konStates, koffvalues,
					timecoursestart, timecourselast, (konrate), dt, t, 
					maxbound2, maxbound3, 0);
			} else {
				x -= (rates.salphc + (konrate));        
				// STOCHASTIC EVENT: a TF unbinds (koff) 
				if (x < rates.koff) {  
          // gets set to true (1) if unbinding failed
          int unbinding_failed = 0;
					tf_unbinding_event(&rates, &state, &genotype, &konStates, koffvalues,
                             timecoursestart,  timecourselast, (konrate), dt, t, x, 0, &unbinding_failed);
				} else {
					//x -= rates.koff;  
					//
					// FALLBACK: shouldn't get here, previous
					// events should be exhaustive
					// 
					printf("ERROR: shouldn't get here!\n");
				}
			}

			// Gillespie step: advance time to next event at dt 
			t += dt;
			LOG_VERBOSE_NOCELLID("dt=%g t=%g\n", dt, t);
			printf("dt=%g t=%g\n", dt, t);
		} else {
			// we will reach the end of development in dt 
			LOG_VERBOSE_NOCELLID("finish at t=%g dt=%g\n", t, dt);

			// do remaining dt 
			dt = tdevelopment - t;

			// advance to end of development (this exits the outer while loop) 
			t = tdevelopment;
		}
	}

	for (i=0; i < NPROTEINS; i++) {
		//if ((output) && j==POP_SIZE-1) print_time_course(timecoursestart[j][i], i, j);
		if ((output)){ 
			print_time_course(timecoursestart[i], i, 0);
		}
		if (verbose) {
			fprintf(fperrors, "deleting gene %d\n", i);
		}
		delete_time_course(timecoursestart[i]);
		timecoursestart[i] = timecourselast[i] = NULL;
	}

	free_mem_CellState(&state);
	free(genotype.all_binding_sites);

	fclose(fperrors);
	//fclose(fp_koff[0]);
	//fclose(fp_tfsbound[0]);
}

int main(int argc, char *argv[])
{
	parse_options(argc, argv);
	system_init();
	model_init();
	run();
}

