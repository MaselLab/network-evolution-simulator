/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
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
#include <omp.h>

/* local includes */
#include "random.h"
#include "lib.h"
#include "netsim.h"

const int MAXELEMENTS=500*MAX_COPIES; 
/* start by allocating maxelements when initializing a genotype, double as needed, reduce at end */
const int MAXBOUND=500*MAX_COPIES;
const int NMIN=4;
const float KRNA=618.0;
const float TTRANSLATION=1.0;
const float TTRANSCRIPTION=1.0;
const float PROB_ACTIVATING=0.62;
const float TRANSCRIPTINIT=8.5; /* replace betaon and betaoff */
const float DEACETYLATE=0.462;
const float ACETYLATE=0.1155;
const float PICASSEMBLY=0.0277;
const float STARTNUCLEUS=0.1;
const float KR=10.0;                  /* don't put this less than 1, weird things happen to koff calculation */
const float GASCONSTANT=8.31447;
const float COOPERATIVITY=1.0;        /* dGibbs, relative to 1 additional specific nt */
const float COOPERATIVE_DISTANCE=11;  /* distance co-operativity operates (changed from 20) */ 
const float NUMSITESINGENOME = 1.3e+6; /* updated from 1.8e+6 */

const float mN = 3.3e-10 * 5e6;    /* Lynch et al. (2008), Tsai et al. (2008)  */

const int avg_protein_conc = 12000;
const float penalty = 0.000000;

/* below are default options, can be changed on command line*/

float kon=1e-4;              /* lower value is so things run faster */
                             /* actual value should be kon=0.2225 is
                                based on 1 molecule taking
                                240seconds=4 minutes and 89% of the
                                proteins being in the nucleus*/
float kon_after_burnin=1e-4; /* lower value value after burn is so things run faster */

int burn_in = 1;             /* disable burn-in by default */

float tdevelopment = 60.0;/* default  development time: can be changed at runtime */
float timemax = -1.0;      /* set an upper limit to development time (default to -1.0=no limit) */
int current_ploidy = 1;    /* ploidy can be changed at run-time: 1 = haploid, 2 = diploid */
int output = 0;
long seed = -58;        
int dummyrun = 0;          /* used to change seed */
int recompute_koff = 1;    /* toggle whether to recompute certain features at each time to avoid
                              compounding rounding error (off by default) */
int recompute_kon = 1;     /* likewise for kon */
float critical_size = 1.0; /* critical size at which cell divides, 
                              set to negative to prevent division  */
float growth_rate_scaling = 2.0; /* set default growth rate scaling factor */
float duration_env0 = 60.1; // in minutes
float duration_env1 = 0.1;
int N_replicates=10;
/* end default options */

/* protein aging term: used when c=c'+g=0, set to 1e-4 < mean-3*sd of
   Belle et al. (2006) and small with respect to most growth rates  */
float protein_aging = 1e-4;

/* file output parameters */
char *output_directory = "output";   /* default output directory */
int verbose = 0;                     /* don't log verbosely by default */ 

/* initialize the growth rate parameters: 
 * do computations here so that we can easily change the scaling factor and Pp */
void initialize_growth_rate_parameters() {
  float hc, gpeak, Ltf;
  gpeak = 0.005776*growth_rate_scaling;  /* in min^-1 based on doubling time of 120 min: ln(2)/(120 min)=0.005776 */
  Pp_a = 12000;
  Pp_b = 12000;              /* mean gene expression of all proteins is 12064.28 */
  Ltf= 1418;               /* mean gene expression of only TFs is 1418 */
  hc = (gpeak/avg_protein_conc)*(1-(log(2-2*0.2)/log(2)));      /* in min^-1 cost of doubling gene expression, based on Wagner (2005) 
                                                   * using {s=0.2, N=500} matches {s=10^-5, N=10^7} combination (both Ns=100) */
  h = hc/0.023;            /* using c=0.023/min from mean of distribution from Belle et al (2006)*/
  gmax_a = gpeak + hc*(Pp_a+(TFGENES*Ltf));    /* compute the gmax coefficient based on gpeak and other parameters */
  gmax_b = gpeak + hc*(Pp_b+(TFGENES*Ltf));
}

char set_base_pair(float x) {
  char base;

  if (x<0.25)
    base = 'a';
  else if (x<0.5)
    base = 'c';
  else if (x<0.75)
    base = 'g';
  else 
    base = 't';
  
  return base;
}

void initialize_sequence(char *Seq, 
                         int len,
                         int ploidy,
                         int num_elements)
{
  float x;
  int i;
  int current_element = len/(num_elements*ploidy);
  int first, second, third, fourth;

  LOG_VERBOSE_NOCELLID("len=%d, NGENES=%d, ploidy=%d, current_element=%d\n", len, NGENES, ploidy, current_element); 
  for (i=0; i<len/ploidy; i++) {
    first = (i / current_element)*ploidy*current_element + i % current_element;
    second = first + current_element;
    third = second + current_element;
    fourth = third + current_element;
    LOG_VERBOSE_NOCELLID("first=%d, second=%d, third=%d, fourth=%d\n", first, second, third, fourth); 
    x = ran1(&seed);
    
    Seq[first] = set_base_pair(x);
    /* clone the randomly chosen sequence for all other sites */
    Seq[second] = Seq[first];
    Seq[third] = Seq[first];
    Seq[fourth] = Seq[first];
  }
  LOG_VERBOSE_NOCELLID("length: %d, sequence is %s\n", strlen(Seq), Seq);
}

void print_genotype(Genotype *genotype, int genotype_id) {
  int i, p;

  printf("[genotype %03d] hind_pos: ", genotype_id);
  for (i=0; i < TFGENES; i++) {
    printf("%d ", genotype->hindrance_positions[i]);
  }
  printf("\n");

  for (i=0; i < NGENES; i++) {
    printf("[genotype %03d gene %02d] ", genotype_id, i);
//    printf("repl=%g, ", genotype->replication_time[i]);
    printf("mRNAdy=%g, ", genotype->mRNAdecay[i]);
    printf("protdy=%g, ", genotype->proteindecay[i]);
    printf("transl=%g, ", genotype->translation[i]);

    for (p=0; p < MAX_COPIES; p++) {
      printf("act[%d]=%d, ", p, genotype->activating[i][p]);
      printf("PICd[%d]=%g ", p, genotype->pic_disassembly[i][p]);
    }
    printf("\n");
  }
}


void print_all_binding_sites(int copies[NGENES],
                             AllTFBindingSites *all_binding_sites, 
                             int numElements,
                             char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                             char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN],
                             int site_id_pos[NGENES][MAX_COPIES][2])
{
  int i, j;

  // TODO: tidy up logic, a bit messy now, but works
  /* loop through the maximum of either the number of genes or number
     of TFs, as they can now be different */
  int max_elements = NGENES >= TFGENES ? NGENES : TFGENES;

  for (i=0; i < max_elements; i++) {
    j=0; 
    /* if copies for this element isn't set, still print out the TFs */
    while (i >= NGENES || j < copies[i]) {
      if (i < TFGENES) 
        printf("TF sequence gene %2d (copy %d): %.*s\n", i, j, TF_ELEMENT_LEN, tf_seq[i][j]);
      else
        printf("            gene %2d (copy %d) does not encode a TF\n", i, j);
      if (i < NGENES) {
        printf("cis-reg     gene %2d (copy %d): %.*s\n", i, j, CISREG_LEN, cisreg_seq[i][j]);
        printf("ID range    gene %2d (copy %d): [%3d, %3d]\n", i, j, site_id_pos[i][j][0], site_id_pos[i][j][1]);
        printf("\n");
      } else {
        printf("            gene %2d (copy %d): no cis-reg gene here\n", i, j);
        printf("\n");
        break;       /* run out of cisreg genes skip to next element */
      }
      j++;
    }
  } 

  printf("numElements: %3d\n", numElements);
  
  for (i=0; i < numElements; i++) {
    printf("binding site %3d:\n", i);
    printf("       cis-reg region: %3d", all_binding_sites[i].cisreg_id);
    printf("         cis-reg copy: %3d", all_binding_sites[i].gene_copy);
    printf(" (sequence %.*s)\n", CISREG_LEN, cisreg_seq[all_binding_sites[i].cisreg_id][all_binding_sites[i].gene_copy]);
    printf(" transcription-factor: %3d", all_binding_sites[i].tf_id);
    printf(" (sequence: %.*s)\n", TF_ELEMENT_LEN, tf_seq[all_binding_sites[i].tf_id][all_binding_sites[i].gene_copy]); 
    printf("  L-edge of %2dbp hind: %3d\n", HIND_LENGTH, all_binding_sites[i].left_edge_pos);        
    printf("  Hind offset position: %3d\n", all_binding_sites[i].hind_pos); 
    printf("               strand: %3d\n", all_binding_sites[i].strand);
    printf("         Hamming dist: %3d\n", all_binding_sites[i].hamming_dist); 
  }
}

#if 0
void print_tf_occupancy(CellState *state,
                        AllTFBindingSites *all_binding_sites,
                        float t)
{
  int bound_count[NGENES][MAX_COPIES];
  int i, j, gene_id, gene_copy;

  for (i = 0; i < NGENES; i++) 
    for (j = 0; j < MAX_COPIES; j++) 
      bound_count[i][j] = 0;
  
  /* for all the currently bound sites */
  for (j = 0; j < state->tf_bound_num; j++) {
    gene_id = all_binding_sites[state->tf_bound_indexes[j]].cisreg_id;
    gene_copy = all_binding_sites[state->tf_bound_indexes[j]].gene_copy;
    bound_count[gene_id][gene_copy]++;
  }
  fprintf(fp_tfsbound[state->cell_id], "%g %d ", t, state->tf_bound_num);
  for (i = 0; i < NGENES; i++) 
    for (j = 0; j < MAX_COPIES; j++) 
      fprintf(fp_tfsbound[state->cell_id], "%d ", bound_count[i][j]);

  fprintf(fp_tfsbound[state->cell_id], "\n");

}

void print_rounding(CellState *state, GillespieRates *rates, float t)
{
  fprintf(fp_rounding[state->cell_id], "%g %d %d %d %d %d %d %d\n", 
          t, rates->koff_operations, rates->transport_operations, rates->mRNAdecay_operations, 
          rates->pic_disassembly_operations, rates->salphc_operations, rates->max_salphc_operations, rates->min_salphc_operations);
}
#endif

void initialize_genotype_fixed(Genotype *genotype, 
                               float kdis[],
                               int genotype_id)
{
  int i, j, p;

  LOG_NOCELLID("[genotype %03d] activators vs repressors ", genotype_id);

  /* initialize hindrance for all TFGENES */
  for (p=0; p < TFGENES; p++) {
    if (HIND_LENGTH == TF_ELEMENT_LEN) {
      genotype->hindrance_positions[p]=0;
    } else  {
      genotype->hindrance_positions[p]=rint(ran1(&seed)*(HIND_LENGTH - TF_ELEMENT_LEN));
    }
  } 
  
  for (i=0; i < NGENES; i++) {
    genotype->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    while (genotype->mRNAdecay[i]<0.0)
      genotype->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    genotype->proteindecay[i]=-1.0;
    while (genotype->proteindecay[i] < 0.0) {
      if (ran1(&seed) < 0.08421)
        genotype->proteindecay[i] = 0.0;
      else genotype->proteindecay[i] = exp(0.7874*gasdev(&seed)-3.7665);
    }
    /* dilution no longer done here, because it is now variable (function of instantaneous growth rate) */
    genotype->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    while (genotype->translation[i] < 0.0)
      genotype->translation[i] = exp(0.7406*gasdev(&seed)+4.56);

    /* make the activations the same in each copy */
    if (ran1(&seed)<PROB_ACTIVATING) {
      for (p=0; p < MAX_COPIES; p++) 
        genotype->activating[i][p] = 1;
    } else {
      for (p=0; p < MAX_COPIES; p++) 
        genotype->activating[i][p] = 0;
    }

    for (p=0; p < MAX_COPIES; p++) 
      LOG_NOFUNC("%d ", genotype->activating[i][p]);

    j = trunc(NUM_K_DISASSEMBLY * ran1(&seed));
    
    for (p=0; p < MAX_COPIES; p++) 
      genotype->pic_disassembly[i][p] = kdis[j];
  }
  LOG_NOFUNC("\n");
 
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, but the same random TF sequence,
 * hindrance positions, replication times, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype, 
                         Genotype *clone,
                         float kdis[],
                         int genotype_id)
{
  int p;
  
  initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*MAX_COPIES*NGENES, MAX_COPIES, NGENES);

  /* within a replicate:
   * initialize all individuals with the same:
   *
   * TF sequence
   * hindrance positions
   * per-gene replication times
   * mRNA and protein decay rates
   * translation rates
   * activators or inhibitors
   * PIC disassembly rates
   */
  // if initializing the first cell, then create from scratch
//  if (genotype_id == 0) {
    initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*MAX_COPIES*TFGENES, MAX_COPIES, TFGENES);
    initialize_genotype_fixed(genotype, kdis, genotype_id);
  // otherwise clone it from the first instance
//  } else {
//    /* copy the TF sequence data */
//    memcpy(genotype->tf_seq, clone->tf_seq, sizeof(char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN]));
//    initialize_new_cell_genotype(genotype, clone);
//  }

  /* start number of copies of gene at current_ploidy */
  for (p=0; p < NGENES; p++) {
    genotype->copies[p] = current_ploidy;
  }

  calc_all_binding_sites(genotype->copies, genotype->cisreg_seq, genotype->tf_seq, 
                         &(genotype->binding_sites_num), &(genotype->all_binding_sites), genotype->hindrance_positions,
                         genotype->sites_per_gene, genotype->site_id_pos);
}

/*
 * mutate a base pair on the cis-regulatory sequence of gene_id and
 * gene_copy with probability m
 */
void mutate(Genotype *gene,
            int gene_id,
            int gene_copy,
            float m)
{
  int k, p;
  char x;
  float q;
  
  for (k=0; k<CISREG_LEN; k++) {
    p = gene_copy;
    if (m > ran1(&seed)) {
      x = gene->cisreg_seq[gene_id][p][k]; /* because sometimes gene and new are the same */
      LOG_VERBOSE_NOCELLID("in gene_id=%2d, mutating pos=%3d, copy=%1d:", gene_id, k, p);
      while (gene->cisreg_seq[gene_id][p][k] == x) {
        q = ran1(&seed);
        gene->cisreg_seq[gene_id][p][k] = set_base_pair(q);
      }
      if (verbose)
        LOG_NOFUNC(" '%c' -> '%c'\n", x, gene->cisreg_seq[gene_id][p][k]);
    }
  }
}

/*
 * compute the list binding sites for specified gene and gene copy
 */
int calc_all_binding_sites_copy(char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN],
                                char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                int binding_sites_num,
                                AllTFBindingSites **all_binding_sites,
                                int *max_alloc,
                                int gene_id,
                                int gene_copy,
                                int hind_pos[TFGENES])
{
  int i, j, tfind, match, max_binding_site_alloc;

  max_binding_site_alloc = *max_alloc;

  for (i=0; i < CISREG_LEN-TF_ELEMENT_LEN; i++) {  /* scan forwards */
    for (tfind=0; tfind < TFGENES; tfind++) {      /* only loop through TF genes */
      match=0;
      for (j=i; j < i+TF_ELEMENT_LEN; j++) {
        if (cisreg_seq[gene_id][gene_copy][j] == tf_seq[tfind][gene_copy][j-i])
          match++;
      }
      if (match >= NMIN){
        if (binding_sites_num + 1 >= max_binding_site_alloc) {
          max_binding_site_alloc = 2*max_binding_site_alloc;
          *all_binding_sites = realloc(*all_binding_sites, max_binding_site_alloc*sizeof(AllTFBindingSites));
          if (!(*all_binding_sites)) {
            LOG_ERROR_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d failed.\n", max_binding_site_alloc);
            exit(1);
          }
          else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d succeeded\n", max_binding_site_alloc);
        }
        if(((i - hind_pos[tfind]) >=0) && (((i - hind_pos[tfind]) + (HIND_LENGTH -1)) < CISREG_LEN)) {
          (*all_binding_sites)[binding_sites_num].cisreg_id = gene_id;
          (*all_binding_sites)[binding_sites_num].gene_copy = gene_copy; 
          (*all_binding_sites)[binding_sites_num].tf_id = tfind;
          (*all_binding_sites)[binding_sites_num].strand = 0;
          (*all_binding_sites)[binding_sites_num].hamming_dist = TF_ELEMENT_LEN-match;
          (*all_binding_sites)[binding_sites_num].hind_pos = hind_pos[tfind];
          (*all_binding_sites)[binding_sites_num].left_edge_pos = i - hind_pos[tfind];
          binding_sites_num++;
        }
      }
    }
  }
  for (i=CISREG_LEN-1; i>=TF_ELEMENT_LEN-1; i--) {  /* scan backwards */
    for (tfind=0; tfind < TFGENES; tfind++) {       /* only loop through TF genes */
      match=0;
      for (j=i; j>i-TF_ELEMENT_LEN; j--)
        if (
            (cisreg_seq[gene_id][gene_copy][j]=='a' && tf_seq[tfind][gene_copy][i-j]=='t')
            || (cisreg_seq[gene_id][gene_copy][j]=='t' && tf_seq[tfind][gene_copy][i-j]=='a')
             || (cisreg_seq[gene_id][gene_copy][j]=='c' && tf_seq[tfind][gene_copy][i-j]=='g')
            || (cisreg_seq[gene_id][gene_copy][j]=='g' && tf_seq[tfind][gene_copy][i-j]=='c')            
            ) match++;
      if (match >= NMIN){
        if (binding_sites_num + 1 >= max_binding_site_alloc){
          max_binding_site_alloc = 2*max_binding_site_alloc;
          *all_binding_sites = realloc(*all_binding_sites, max_binding_site_alloc*sizeof(AllTFBindingSites));
          if (!(*all_binding_sites)){
            LOG_ERROR_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d failed.\n", max_binding_site_alloc);
            exit(1);
          }
          else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d succeeded\n", max_binding_site_alloc);
        }
        if (((i-TF_ELEMENT_LEN+1 - hind_pos[tfind]) >=0) && (((i-TF_ELEMENT_LEN+1 - hind_pos[tfind])+(HIND_LENGTH-1))< CISREG_LEN)) {
          (*all_binding_sites)[binding_sites_num].cisreg_id = gene_id;
          (*all_binding_sites)[binding_sites_num].gene_copy = gene_copy; 
          (*all_binding_sites)[binding_sites_num].tf_id = tfind;
          (*all_binding_sites)[binding_sites_num].strand = 1;
          (*all_binding_sites)[binding_sites_num].hamming_dist = TF_ELEMENT_LEN-match;
          (*all_binding_sites)[binding_sites_num].hind_pos = hind_pos[tfind];
          (*all_binding_sites)[binding_sites_num].left_edge_pos = i-TF_ELEMENT_LEN+1 - hind_pos[tfind];
          binding_sites_num++;
        }
      }
    }
  }
  *max_alloc = max_binding_site_alloc;
  return binding_sites_num;
}

/*
 * compute the list of binding sites for the specified number of gene
 * copies
 */
void calc_all_binding_sites(int copies[NGENES],
                            char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN],
                            char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                            int *new_binding_sites_num,
                            AllTFBindingSites **all_binding_sites,
                            int hind_pos[TFGENES],
                            int sites_per_gene[NGENES],
                            int site_id_pos[NGENES][MAX_COPIES][2])
{
  int p, max_binding_site_alloc, binding_sites_num;
  int gene_id;

  max_binding_site_alloc = MAXELEMENTS;
  *all_binding_sites = malloc(max_binding_site_alloc*sizeof(AllTFBindingSites));
  if (!(*all_binding_sites)) {
    LOG_ERROR_NOCELLID("initial setting of all_binding_sites failed.\n");
    exit(1);
  }
  binding_sites_num = 0;

  /* initialize per gene # of sites */
  for (gene_id=0; gene_id < NGENES; gene_id++) {
    sites_per_gene[gene_id] = 0;
    for (p=0; p < MAX_COPIES; p++) {
      site_id_pos[gene_id][p][0] = -1;
      site_id_pos[gene_id][p][1] = -1;
    }
  }

  for (gene_id=0; gene_id < NGENES; gene_id++) {  /* now which cis-reg region */
    for (p=0; p < MAX_COPIES; p++) {     /* loop through the maximum copies possible */
      if (p < copies[gene_id]) {  
        int before = binding_sites_num;

        /* record initial site_id for this copy */
        site_id_pos[gene_id][p][0] = binding_sites_num;

        /* if this particular gene has this copy number then generate binding sites
           for all the relevant gene copies (assume no gene divergence) */
        binding_sites_num = calc_all_binding_sites_copy(cisreg_seq, 
                                                    tf_seq, 
                                                    binding_sites_num,
                                                    all_binding_sites,
                                                    &max_binding_site_alloc,
                                                    gene_id,
                                                    p, 
                                                    hind_pos);

        /* record end site_id for this copy */
        site_id_pos[gene_id][p][1] = binding_sites_num - 1;

        /* add the new number of sites */
        sites_per_gene[gene_id] += (binding_sites_num-before);
      }
    }
  }

  *all_binding_sites = realloc(*all_binding_sites, binding_sites_num*sizeof(AllTFBindingSites));
  if (!(*all_binding_sites)) {
    LOG_ERROR_NOCELLID("realloc of all_binding_sites down to binding_sites_num = %d failed.\n", binding_sites_num);
    exit(1);
  }
  *new_binding_sites_num = binding_sites_num;
}


int add_fixed_event(int i,
                    int p,
                    float t,
                    FixedEvent **start,
                    FixedEvent **last)
{
  FixedEvent *newtime;
  int pos;

  newtime = (FixedEvent *)malloc(sizeof(FixedEvent));
  if (!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->gene_id = i;
  newtime->copy = p;
  newtime->time = t;
  LOG_VERBOSE_NOCELLID("adding event at time=%f for gene=%d (copy=%d)\n", t, i, p);
  pos = sls_store(newtime, start, last);
  return pos;
}

void add_time_point(float time,
                    float conc,
                    TimeCourse **start,
                    TimeCourse **last)
{
  TimeCourse *newtime;
  
  newtime = (TimeCourse *)malloc(sizeof(TimeCourse));
  if (!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->time = time;
  newtime->concentration = conc;
  sls_store_end2(newtime, start, last);
}

void add_fixed_event_end(int gene_id,
                         int gene_copy,
                         float t,
                         FixedEvent **start,
                         FixedEvent **last)
{
  FixedEvent *newtime;
  
  newtime = (FixedEvent *)malloc(sizeof(FixedEvent));
  if (!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->gene_id = gene_id;
  newtime->copy = gene_copy;
  newtime->time = t;
  LOG_VERBOSE_NOCELLID("adding event end at time=%f for gene=%d (copy=%d)\n", t, gene_id, gene_copy);
  sls_store_end(newtime, start, last);
}

void delete_fixed_event(int gene_id,
                        int gene_copy,
                        int i,
                        FixedEvent **start,
                        FixedEvent **last)
{
  FixedEvent *info, *lastinfo;
  int j, done;
  
  j = -1;
  done = 0;
  info = *start;
  while (info) {
    if ((info->gene_id==gene_id && info->copy==gene_copy)) {
      j++;
      if (j == i) {
        if (info == *start) {
          *start = info->next;
          if (info == *last) *last = NULL;
        } else {
          lastinfo->next = info->next;
          if (info == *last) *last = lastinfo;
        }
        done = 1;
      } else {
        lastinfo = info;
        info = info->next;
      }
    } else {
      lastinfo = info;
      info = info->next;
    }
  }
  if (done == 0) {
    LOG_ERROR_NOCELLID("In %d elements, couldn't find element %d to delete in gene %d (copy %d)\n",
                       j+1, i, gene_id, gene_copy);
  }
  free(info);
}

void delete_fixed_event_start(FixedEvent **start,
                              FixedEvent **last)
{
  FixedEvent *info;
  
  info = *start;
  *start = info->next;
  if (*last == info) *last = NULL;
  free(info);
}

/*
 * initialize the cell state with the specified initial protein
 * concentration, mean mRNA number and mRNA decay and whether to do
 * burn-in of high kon rate or not
 */
void initialize_cell(CellState *state,
                     int id,
                     int copies[NGENES],
                     float mRNAdecay[NGENES],
                     float meanmRNA[NGENES],
                     float init_protein_conc[NPROTEINS],
                     int burn_in)
{
  int i, j, k, totalmRNA;
  float t;

  /* initialize the ID of the cell */
  state->cell_id = id;

  /* initialize whether to do kon burn-in or not */
  state->burn_in = burn_in;

  /* start cell size at 0.5 */
  state->cell_size = 1.0; //04/20 KUN

  /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
  state->growth_rate = 0.0;

  state->mRNA_transcr_time_end = NULL;
  state->mRNA_transcr_time_end_last = NULL;
  state->mRNA_transl_time_end = NULL;
  state->mRNA_transl_time_end_last = NULL;
  state->env0_time_end = NULL;
  state->env0_time_end_last = NULL;
  state->env1_time_end = NULL;
  state->env1_time_end_last = NULL;

  state->tf_bound_num = 0;  /* initialize with nothing bound */
  state->tf_hindered_num = 0;
   
//  state->tf_bound_indexes = NULL; // moved to the beginning of calc_avg_growth_rate
//  state->tf_hindered_indexes = NULL;

  for (i=0; i < NGENES; i++) {

    for (j=0; j < MAX_COPIES; j++) {
      state->active[i][j] = ON_WITH_NUCLEOSOME;
    }

    totalmRNA = (int) poidev(meanmRNA[i],&seed);
    state->mRNA_nuclear_num[i] = (int) bnldev(STARTNUCLEUS, totalmRNA, &seed);
    state->mRNA_cyto_num[i] = totalmRNA - state->mRNA_nuclear_num[i];
    state->mRNA_transl_cyto_num[i] = 0;

    for (k=0; k<state->mRNA_cyto_num[i]; k++) {
      t = expdev(&seed) / mRNAdecay[i];
      if (t < TTRANSLATION) {
        (state->mRNA_cyto_num[i])--;
        (state->mRNA_transl_cyto_num[i])++;
        LOG_VERBOSE("add translation event time=%g for gene=%d\n", (TTRANSLATION-t), i);
        add_fixed_event(i, -1, TTRANSLATION-t, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
      }
    } 

    int total_mRNA_transcribing = (int) poidev(meanmRNA[i]*TTRANSCRIPTION*mRNAdecay[i], &seed);
    
    /* split it up evenly between the copies */
    int mRNA_copy1 = trunc(total_mRNA_transcribing/current_ploidy);
    int mRNA_copy2 = total_mRNA_transcribing - mRNA_copy1;

    for (j=0; j < MAX_COPIES; j++) {
      if (j < current_ploidy)  {
        state->mRNA_transcr_num[i][j] = (j==0) ? mRNA_copy1 : mRNA_copy2;
        LOG_VERBOSE_NOCELLID("initializing state->mRNA_transcr_num[%2d][%2d]=%d\n", i, j, state->mRNA_transcr_num[i][j]);
        for (k=0; k < state->mRNA_transcr_num[i][j]; k++)
          add_fixed_event(i, j, ran1(&seed)*TTRANSCRIPTION, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));
      } else {
        state->mRNA_transcr_num[i][j] = 0;
      }
    }
  }
  for (i=0; i < NPROTEINS; i++) {
    state->protein_conc[i] = init_protein_conc[i];
  }
  // TODO: add more env changing 
  add_fixed_event(0,0,duration_env0,&(state->env0_time_end),&(state->env0_time_end_last));
}

/*
 * initialize the cell state that are "cached", i.e. that are
 * maintained for code efficiency and are not part of the underlying
 * molecular biology.  All of these data structures can be regenerated
 * from the cell state at any point in the simulation
 */
void initialize_cell_cache(CellState *state,
                           Genotype genes,
                           KonStates *kon_states,
                           float **koffvalues,
                           int maxbound2,
                           int maxbound3)
{
  int i,j;
  /* number of possible binding sites */
  /* currently create for all proteins, not just TFs, as we need info for protein decay */
  for (i=0; i < NPROTEINS; i++){
    kon_states->kon_list[i] = malloc(sizeof(KonList));
    kon_states->kon_list[i]->available_sites = malloc(genes.binding_sites_num*sizeof(int));
    for(j=0;j<genes.binding_sites_num;j++){   
        kon_states->kon_list[i]->available_sites[j]=-1;
    }
   }
  
 state->tf_bound_indexes = realloc(state->tf_bound_indexes,maxbound2*sizeof(int));
 
 *koffvalues = malloc(maxbound2*sizeof(float)); 
   for(i=0;i<maxbound2;i++){
      state->tf_bound_indexes[i]=-1;
      (*koffvalues)[i]=0.0;
  } 
  
  state->tf_hindered_indexes = realloc(state->tf_hindered_indexes, 2*maxbound3*sizeof(int));
 
  for(i=0;i<2*maxbound3;i++){
      state->tf_hindered_indexes[0][i]=-1;
      state->tf_hindered_indexes[1][i]=-1;
  }  
  if (!kon_states->konvalues || !state->tf_bound_indexes || !koffvalues ||
      !state->tf_hindered_indexes || !kon_states->kon_list) {
    LOG_ERROR("memory allocation error at start of develop\n");
    exit(1);
  }
}

/* could perhaps be a little faster with option to skip *df calculation for first 2 calls */
void calc_time (float t, 
                float x, 
                GillespieRates *rates,
                KonStates *kon_states,
                float *f, 
                float *df)
{
  float r, denom, numer, ct, ect;
  int i;
  //printf("x=%g, t=%g\n", x, t);
  r = numer = 0.0;

  /* loop over all only TFs that bind */
  for (i=0; i < TFGENES; i++) {

    LOG_VERBOSE_NOCELLID("t=%g, kon_states->kon_list[%d]->site_count=%d, kon_states->konvalues[%d][KON_DIFF_INDEX]=%g\n", 
                         t, i, kon_states->kon_list[i]->site_count, i, kon_states->konvalues[i][KON_DIFF_INDEX]);

    /* if this particular TF is bound somewhere */
    if (kon_states->kon_list[i]->site_count > 0) {
      ct = kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX] * t;
      if (fabs(ct)<EPSILON) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) kon_states->kon_list[i]->site_count) * kon_states->konvalues[i][KON_DIFF_INDEX] * ect;
      numer += ((float) kon_states->kon_list[i]->site_count) * kon_states->konvalues[i][KON_DIFF_INDEX] * (ect-ct*exp(-ct));
    }
  }
  numer *= kon;
  r *= kon;
  LOG_VERBOSE_NOCELLID("r=%g \n", r);
  denom = r + t*(rates->subtotal + rates->salphc);
  denom = denom * denom;
  r /= t;

  /* compute delta-t */
  *f = x/(r + rates->subtotal + rates->salphc) - t;

  LOG_VERBOSE_NOCELLID("x=%g, r=%g, rates->subtotal=%g, rates->salphc=%g, f=%g\n", x, r, rates->subtotal, rates->salphc, *f);

  /* compute derivative of equation */
  *df = x*numer/denom - 1.0;

  LOG_VERBOSE_NOCELLID("t=%g f=%g df=%g\n", t, *f, *df);
  //printf("x=%g, t=%g f=%g df=%g\n",x,  t, *f, *df);
  //system("PAUSE");
}

void calc_kon_rate(float t,
                   KonStates *kon_states,
                   float *konrate)
{
  float r,ct=0.0,ect=0.0;
  int i;
  
  r = 0.0;
  /* loop through all TFs */
  for (i=0; i < TFGENES; i++) {  
    if (kon_states->kon_list[i]->site_count > 0) {
      ct = kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX]*t;
      if (fabs(ct)<EPSILON) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) kon_states->kon_list[i]->site_count)*kon_states->konvalues[i][KON_DIFF_INDEX]*ect;
    }
  }
  *konrate = kon*r/t;
  LOG_VERBOSE_NOCELLID("r=%g t=%g konrate=%g\n", r, t, *konrate);
}

/* 
 * change in the number of mRNAs in the cytoplasm affects the kon
 * rates, note we must have already updated protein_conc first
 */
void change_mRNA_cytoplasm(int i,
                           Genotype *genotype,
                           CellState *state,
                           GillespieRates *rates,
                           KonStates *kon_states)
{
  float salphc; 
  
  /* number of mRNAs in cytoplasm affects kon rates */
  salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  
  LOG_VERBOSE("change_mRNA_cytoplasm[%d]: mRNA=%d, transl rate=%g, protein decay=%g, salphc=%g\n", 
              i, state->mRNA_cyto_num[i], genotype->translation[i], kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX], salphc);
  
  rates->salphc += kon_states->kon_list[i]->site_count*kon*(salphc - kon_states->konvalues[i][KON_SALPHC_INDEX]);
  rates->salphc_operations++;

  rates->max_salphc += kon_states->kon_list[i]->site_count*kon*(fmaxf(state->protein_conc[i], salphc) - fmaxf(state->protein_conc[i], kon_states->konvalues[i][KON_SALPHC_INDEX]));
  rates->max_salphc_operations++;
  rates->min_salphc += kon_states->kon_list[i]->site_count*kon*(fminf(state->protein_conc[i], salphc) - fminf(state->protein_conc[i], kon_states->konvalues[i][KON_SALPHC_INDEX]));    
  rates->min_salphc_operations++;

  kon_states->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  kon_states->konvalues[i][KON_SALPHC_INDEX] = salphc;
}

/*
 * calculate the koff value for the specified binding site, k
 */
void calc_koff(int k,
               AllTFBindingSites *all_binding_sites,
               CellState *state,
               float *koff,
               float t)
{
  float Gibbs;  /*free energy in kJ/mol*/
  int posdiff, front, back, j;
  
  front = back = 0;
  Gibbs = (((float) all_binding_sites[k].hamming_dist)/3.0 - 1.0) * state->RTlnKr; /* subject to revision of TF_ELEMENT_LEN */
  
  for (j=0; j < state->tf_bound_num; j++) { //search for any cooperation or hindrance
      
    if (all_binding_sites[k].cisreg_id==all_binding_sites[state->tf_bound_indexes[j]].cisreg_id && //which gene this bs is on 
        all_binding_sites[k].gene_copy==all_binding_sites[state->tf_bound_indexes[j]].gene_copy && // and which copy
        !(k==state->tf_bound_indexes[j])) 
    {
      posdiff = all_binding_sites[k].left_edge_pos - all_binding_sites[state->tf_bound_indexes[j]].left_edge_pos;
      
      if (abs(posdiff) < HIND_LENGTH) 
      {
        LOG_ERROR("t=%g steric hindrance breached with site %d (at pos %d, strand %d copy %d gene %d), %d away from site %d \
                   (at pos %d, strand %d copy %d of gene %d)\n",
                  t, k, all_binding_sites[k].left_edge_pos, all_binding_sites[k].strand, all_binding_sites[k].gene_copy, all_binding_sites[k].cisreg_id,
                  posdiff, state->tf_bound_indexes[j], all_binding_sites[state->tf_bound_indexes[j]].left_edge_pos, 
                  all_binding_sites[state->tf_bound_indexes[j]].strand, all_binding_sites[state->tf_bound_indexes[j]].gene_copy, 
                  all_binding_sites[state->tf_bound_indexes[j]].cisreg_id);
        exit(-1);
      }
      
      if (abs(posdiff) < COOPERATIVE_DISTANCE) 
      { 
            if (posdiff>0) front++; else back++;
      }
    }
    
    if ((front) && (back)) // if there's cooperation
      j=state->tf_bound_num;
  }
  if (front>0) Gibbs -= COOPERATIVITY*state->RTlnKr/3; 
  if (back>0) Gibbs -= COOPERATIVITY*state->RTlnKr/3;
  *koff = NUMSITESINGENOME*kon*0.25/exp(-Gibbs/(GASCONSTANT*state->temperature));
  LOG_VERBOSE("state->RTlnKr=%g front=%d back=%d H=%d Gibbs=%g koff=%g\n",state->RTlnKr,front,back,all_binding_sites[k].hamming_dist,Gibbs,*koff);
  /* 25% protein in nucleus is implicit in formula above */
}

/* 
 * when TF binding changes, adjust cooperativity at neighbouring sites 
 */
void scan_nearby_sites(int indexChanged,
                       AllTFBindingSites *all_binding_sites,
                       CellState *state,
                       GillespieRates *rates,
                       float *koffvalues,
                       float t)
{
  int posdiff, j;
  float diff;
  
  /* for all the currently bound sites */
  for (j = 0; j < state->tf_bound_num; j++) {
    /* we are on the same cisreg sequence and we aren't the same TF binding  */
    if (all_binding_sites[indexChanged].cisreg_id == all_binding_sites[state->tf_bound_indexes[j]].cisreg_id && 
        all_binding_sites[indexChanged].gene_copy == all_binding_sites[state->tf_bound_indexes[j]].gene_copy && 
        !(indexChanged==state->tf_bound_indexes[j])) {

      /* how close are we on the sequence */
      posdiff = all_binding_sites[indexChanged].left_edge_pos - all_binding_sites[state->tf_bound_indexes[j]].left_edge_pos;
       if (abs(posdiff) < HIND_LENGTH) { /* within HIND_LENGTH: bad: shouldn't happen */
        LOG_ERROR("t=%g steric hindrance 2 has been breached with site %d %d away from site %d\n",
                  t, indexChanged, posdiff, state->tf_bound_indexes[j]);
        exit(-1);
      }
      if (abs(posdiff) < COOPERATIVE_DISTANCE) {  /* within cooperative distance adjust koff */

        diff = -koffvalues[j];                /* save old value */
        calc_koff(state->tf_bound_indexes[j], /* recompute koffvalues */
                  all_binding_sites, state, 
                  &(koffvalues[j]), t);          
        diff += koffvalues[j];                /* calculating how koff changes  */
        rates->koff += diff;                  /* adjust rates by difference */
        rates->koff_operations++;
      }
    }
  }
}

/*
 * remove binding site site_id from the pool of available binding sies
 */
void remove_kon(int site_id,
                int TFID,
                GillespieRates *rates,
                float salphc,
                KonStates *kon_states,
                float protein_concTFID)
{
  int k;
  k = 0;
  
  while (!(kon_states->kon_list[TFID]->available_sites[k] == site_id) && k < kon_states->kon_list[TFID]->site_count) {
    k++;
  }

  LOG_VERBOSE_NOCELLID(">>> remove site %d kon_list (k=%d of %d total sites for TF %d, grandtotal=%d)\n", 
                       site_id, k, kon_states->kon_list[TFID]->site_count, TFID, kon_states->nkon);

  /* make sure that we have enough unoccupied sites left */
  if (k < kon_states->kon_list[TFID]->site_count && k < kon_states->nkon) { 
    /* adjust rates */
    rates->salphc -= kon*salphc;
    rates->salphc_operations++;

    rates->max_salphc -= kon*fmaxf(protein_concTFID, salphc);
    rates->max_salphc_operations++;

    rates->min_salphc -= kon*fminf(protein_concTFID, salphc);
    rates->min_salphc_operations++;

    /* one less site available for binding of total */
    (kon_states->nkon)--;

    /* also per gene */
    (kon_states->kon_list[TFID]->site_count)--;

    /* move the last element end of array into space vacated by site k */
    kon_states->kon_list[TFID]->available_sites[k] = kon_states->kon_list[TFID]->available_sites[kon_states->kon_list[TFID]->site_count];
    kon_states->kon_list[TFID]->available_sites[kon_states->kon_list[TFID]->site_count]=-1;
  } else {
    LOG_VERBOSE_NOCELLID("||| couldn't remove site %d from TF %d in kon_list (k=%d)\n", site_id, TFID, k);
  }
  /* else do nothing: there is likely a redundancy in steric
     hindrance, hence no site to remove */
}

/*
 * add binding site site_id back to the pool of available binding sies
 */
void add_kon(float protein_concTFID,
             float salphc,
             int TFID,
             int site_id,
             GillespieRates *rates,
             KonStates *kon_states)
{

  /* update rates because new site is now available */
  rates->salphc += kon*salphc;
  rates->salphc_operations++;

  rates->max_salphc += fmaxf(protein_concTFID, salphc);
  rates->max_salphc_operations++;

  rates->min_salphc += fminf(protein_concTFID, salphc);
  rates->min_salphc_operations++;

  /* add back site_id to pool of available sites */
  kon_states->kon_list[TFID]->available_sites[kon_states->kon_list[TFID]->site_count] = site_id;

  /* one more site available */
  (kon_states->kon_list[TFID]->site_count)++;
  (kon_states->nkon)++;
}

/* 
 * for specified gene_id and gene_copy, tests whether criterion for
 * transcription is met
 */
int ready_to_transcribe(int gene_id,
                        int gene_copy,
                        int *tf_bound_indexes,
                        int tf_bound_num,
                        AllTFBindingSites *all_binding_sites,
                        int activating[NGENES][MAX_COPIES],
                        int *on)
{
  int i, off;
  
  *on=off=0;
  for (i=0; i < tf_bound_num; i++) {
    if (gene_id==all_binding_sites[tf_bound_indexes[i]].cisreg_id &&
        gene_copy==all_binding_sites[tf_bound_indexes[i]].gene_copy)
      {
        if (activating[all_binding_sites[tf_bound_indexes[i]].tf_id][gene_copy]) (*on)++;
        else off++;
    }
  }
  if ((float)off <= 0.33442*(float)(*on) + 0.31303) 
    return (1);
  else 
    return (0);
}

/* 
 * returns true if at least one activator is bound
 */
int is_one_activator(int gene_id,
                     int gene_copy,
                     int *tf_bound_indexes,
                     int tf_bound_num,
                     AllTFBindingSites *all_binding_sites,
                     int activating[NGENES][MAX_COPIES])
{
  int i;
  
  for (i=0; i < tf_bound_num; i++)
    if (gene_id==all_binding_sites[tf_bound_indexes[i]].cisreg_id && 
        gene_copy==all_binding_sites[tf_bound_indexes[i]].gene_copy &&
        (activating[all_binding_sites[tf_bound_indexes[i]].tf_id][gene_copy])) 
      return (1);
  return (0);
}

/* 
 * calculate the rates based upon the initialization of the genotype
 * and the cell state.  Note this is only appropriate if nothing is
 * bound
 */
void calc_from_state(Genotype *genotype,
                     CellState *state,
                     GillespieRates *rates,
                     KonStates *kon_states,
                     float transport_rate[NGENES],
                     float mRNAdecay[NGENES]) 
{
  int i, j, k;
  float salphc, protein_conc_for_tfid; 
  float protein_decay;

  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to fall below cut-off, use aging term */
    
    protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;
    
    salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);
    
    kon_states->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / (protein_decay);
    
    kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
    
    kon_states->konvalues[i][KON_SALPHC_INDEX] = salphc;
   
    kon_states->kon_list[i]->site_count = 0;
    
    LOG_VERBOSE("protein decay[%d]=%g\n", i, kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX]);
  }
  
  state->tf_bound_num=0;
  rates->koff=0.0;
  rates->koff_operations = 0;
  rates->transport=0.0;
  rates->transport_operations = 0;
  rates->mRNAdecay=0.0; // these rates are calced in calc_dt  
  rates->mRNAdecay_operations = 0;
  rates->pic_disassembly=0.0;
  rates->pic_disassembly_operations = 0;
  rates->salphc=0.0;
  rates->salphc_operations=0;
  rates->max_salphc=0.0;
  rates->max_salphc_operations=0;
  rates->min_salphc=0.0;
  rates->min_salphc_operations=0;

  for (k=0; k < genotype->binding_sites_num; k++) {
    i = genotype->all_binding_sites[k].tf_id;
    protein_conc_for_tfid = state->protein_conc[i];
    salphc = kon_states->konvalues[i][KON_SALPHC_INDEX];

    rates->salphc += salphc;
    rates->salphc_operations++;
    rates->max_salphc += fmaxf(protein_conc_for_tfid, salphc);
    rates->max_salphc_operations++;
    rates->min_salphc += fminf(protein_conc_for_tfid, salphc);
    rates->min_salphc_operations++;

    /* update the list of sites that bind for a particular TF, i */
    kon_states->kon_list[i]->available_sites[kon_states->kon_list[i]->site_count] = k;
    (kon_states->kon_list[i]->site_count)++;
  }

  /* initialize kon_states->nkon as the total number of binding sites */
  kon_states->nkon = genotype->binding_sites_num;

  for (i=0; i < NGENES; i++) {
  	// TODO: fix this!
    // FIXME: this updates kon_states for non-TFs, even though this isn't used
//    LOG("after initializing kon_states for gene=%d, site_count=%d, sites_per_gene=%d, nkon=%d\n", 
//        i, kon_states->kon_list[i]->site_count, genotype->sites_per_gene[i], kon_states->nkon);

    transport_rate[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
//    printf("in calc-from-state: N_mRNA %d = %d, rate= %f, the address is %i\n", i, state->mRNA_nuclear_num[i], transport_rate[i],&transport_rate[i]);
    rates->transport += transport_rate[i];
    rates->transport_operations++;
//    LOG_VERBOSE("mRNA_nuclear_num=%d, initializing transport[%d]=%g\n", state->mRNA_nuclear_num[i], i, transport[i]);
  }
//  LOG_VERBOSE("initializing rates->transport=%g\n", rates->transport);

  /* start all genes in acteylated state */
  for (j=0; j < MAX_COPIES; j++) {  
    int pos = 0;
    for (i=0; i < NGENES; i++) {
      if (genotype->copies[i] > j) {  
        state->state_change_ids[ACETYLATION_STATE][j][pos] = i;
        LOG_VERBOSE("Initializing statechange gene=%d, ploidy=%d state_change_ids[%d][%d]=%d\n", i, j, j, 
                    pos, state->state_change_ids[ACETYLATION_STATE][j][pos]);
        pos++;
      }
    }
  }

  /* scale up all salphc rates by the global kon value */
  rates->salphc *= kon;
  rates->salphc_operations++;
  rates->max_salphc *= kon;
  rates->max_salphc_operations++;
  rates->min_salphc *= kon;
  rates->min_salphc_operations++;

  /* first initialize everything at zero */
  for (j=0; j < MAX_COPIES; j++) {
    rates->acetylation_num[j]=0;
    rates->deacetylation_num[j]=0;
    rates->pic_assembly_num[j]=0;
    rates->transcript_init_num[j]=0;
    rates->pic_disassembly_num[j]=0;
  }

  /* now set the per copy acetylation rate  */
  for (i=0; i < NGENES; i++) {
    for (j=0; j < genotype->copies[i]; j++) {
      rates->acetylation_num[j]++;
    }
  }

  if (verbose) 
    for (j=0; j < MAX_COPIES; j++) {
      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
    }
}

/* 
 * check to see if a fixed event ends within the time t
 *
 * returns:
 *  0 if there is no fixed event occuring before time t
 *  1 if a transcription event happens before time t
 *  2 if a translation event happens before time t
 *  3 if a gene replication event happens before time t
 */
 
int does_fixed_event_end(FixedEvent *mRNA_transl_time_end,
                         FixedEvent *mRNA_transcr_time_end,
                         FixedEvent *env0_time_end,
                         FixedEvent *env1_time_end,
			 float t) 
{
	int retval;
	float t1;
	float t2;
	float t3;
	float t4;
	
	if(mRNA_transcr_time_end == NULL && mRNA_transl_time_end==NULL && env0_time_end == NULL && env1_time_end == NULL){
		retval =0;
	}
	else{
		t1 = mRNA_transcr_time_end ? mRNA_transcr_time_end->time : TIME_INFINITY;
		t2 = mRNA_transl_time_end ? mRNA_transl_time_end->time : TIME_INFINITY;
		t3 = env0_time_end ? env0_time_end->time : TIME_INFINITY;
		t4 = env1_time_end ? env1_time_end->time : TIME_INFINITY;
		
		if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4)){
//			if (mRNA_transcr_time_end == NULL) retval = 0;
//      		else 
			  retval = 1;	
		}
		else{
			if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4)) {
//		        if (mRNA_transl_time_end == NULL)  retval = 0;
//		        else 
				retval = 2;
			}
			else{
				if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4)) {
//			        if (env0_time_end == NULL)  retval = 0;
//			        else 
					retval = 3;
			    }
			    else{
			    	if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3)) {
//			        if (env1_time_end == NULL)  retval = 0;
//			        else 
					retval = 4;
			    	}
			    	else{
						retval = 0;
					}
			    }
			}
		}	
	 							
	}
	return retval;
}


/*
 * calculate the length of the next (Gillespie) timestep, dt
 */
void calc_dt(float *x,
             float *dt,
             GillespieRates *rates,
             KonStates *kon_states,
             float mRNAdecay[],
             float mRNAdecayrates[],
             int mRNA_cyto_num[],
             int mRNA_transl_cyto_num[],
             int cell_id)
{
  float tbound1, tbound2;
  int i, j;

  /* reset the subtotal rate (excludes konrate) for current step */
  rates->subtotal=0.0;
  /* reset mRNA decay rate */
  rates->mRNAdecay=0.0;
  /* hence reset the number of rounding operations for this particular rate */
  rates->mRNAdecay_operations=0;

  /* update mRNAdecay rate based on the total number of mRNAs in both
     cytoplasm (mRNA_cyto_num) and ones that have only just recently arrived
     (mRNA_transl_cyto_num) */
  for (i=0; i < NGENES; i++) {
    mRNAdecay[i] = mRNAdecayrates[i] * ((float) mRNA_cyto_num[i] + (float) mRNA_transl_cyto_num[i]);
    rates->mRNAdecay += mRNAdecay[i];
    rates->mRNAdecay_operations++;
  }

  /* recompute and cache the total rate in data structure */
  rates->subtotal += rates->koff;
  rates->subtotal += rates->transport;
  rates->subtotal += rates->mRNAdecay;
  rates->subtotal += rates->pic_disassembly;
  rates->subtotal += rates->salphc;
  //LOG_ERROR_NOCELLID("rates subtotal = %f, salphc = %f\n", rates->subtotal, rates->salphc);
  //printf("rates subtotal = %f\n rates koff = %f\n rates transport = %f\n rates mrna = %f\n rates pic = %f\n rates salphc = %f \n",rates->subtotal, rates->koff, rates->transport, rates->mRNAdecay, rates->pic_disassembly, rates->salphc);
  /* 
   * convert the counts back into rates using the constants 
   */
  for (j=0; j < MAX_COPIES; j++) {
    rates->subtotal += (float) rates->acetylation_num[j] * ACETYLATE;
    rates->subtotal += (float) rates->deacetylation_num[j] * DEACETYLATE;
    rates->subtotal += (float) rates->pic_assembly_num[j] * PICASSEMBLY;
    rates->subtotal += (float) rates->transcript_init_num[j] * TRANSCRIPTINIT;    
  } 
  //printf("rates subtotal = %f\n",rates->subtotal);
  //system("PAUSE");
  tbound1 = *x/(rates->subtotal + rates->max_salphc);
  tbound2 = *x/(rates->subtotal + rates->min_salphc);
  LOG_VERBOSE_NOCELLID("[cell %03d] bounds %g %g\n", cell_id, tbound1, tbound2);

  /* if bounds are the same, simply choose tbound1 */
  if (tbound1==tbound2){
    if (kon_states->nkon!=0) {
      LOG_ERROR_NOCELLID("[cell %03d] nkon=%d when it should be zero x=%f rates->max_salphc=%g rates->min_salphc=%g rates->subtotal=%g\n",
                         cell_id, kon_states->nkon, *x, rates->max_salphc, rates->min_salphc, rates->subtotal);
    }
    *dt = tbound1;
  } else {
    /* otherwise get delta t by solving the equation using Newton-Raphson method */
    *dt = rtsafe(&calc_time, *x, rates, kon_states, tbound1, tbound2, (float) RT_SAFE_EPSILON); 
  }
}

/*
 * end transcription: update the mRNAs in the nucleus, cytoplasm
 * etc. accordingly and delete the event from the queue
 */
void end_transcription(float *dt,
                       float t,
                       CellState *state,
                       float transport[NGENES],
                       GillespieRates *rates)
{
  int i, j, total;
  
  /* recompute the delta-t based on difference between now and the
     time of transcription end */
  *dt = state->mRNA_transcr_time_end->time - t;

  if (verbose) {
    total = 0;
    for (i=0; i < NGENES; i++) 
      for (j=0; j < MAX_COPIES; j++) 
        total += state->mRNA_transcr_num[i][j];
    
    LOG_VERBOSE("\ntranscription event finishes out of %d possible t=%g dt=%g\n", total, t, *dt);
  }

  /* get the gene which is ending transcription */
  i = state->mRNA_transcr_time_end->gene_id;
  j = state->mRNA_transcr_time_end->copy;

  /* increase number of mRNAs in nucleus */
  (state->mRNA_nuclear_num[i])++;

  /* decrease the number of mRNAs undergoing transcription */
  (state->mRNA_transcr_num[i][j])--;

  /* delete the fixed even which has just occurred */
  delete_fixed_event_start(&(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));

  /* add rate KRNA to transport and Gillespie rates */
  transport[i] += KRNA;
  rates->transport += KRNA;
  rates->transport_operations++;

  LOG_VERBOSE("add one new mRNA in nucleus, updating transport[%d]=%g, rates->transport=%g\n", i, transport[i], rates->transport);

}


/*
 * do the actual disassembling of the pre-initiation complex
 */
void disassemble_PIC(CellState *state,
                     Genotype *genotype,
                     int gene_id,
                     int gene_copy,
                     GillespieRates *rates)
{
  float disassembly = genotype->pic_disassembly[gene_id][gene_copy];
  LOG_ERROR("PIC DIS2\n");
  remove_from_array(gene_id, TRANSCRIPTINIT_STATE, state->state_change_ids[TRANSCRIPTINIT_STATE][gene_copy], 
                    &(rates->transcript_init_num[gene_copy]), (int) 1);
  remove_from_array(gene_id, PICDISASSEMBLY_STATE, state->state_change_ids[PICDISASSEMBLY_STATE][gene_copy], 
                    &(rates->pic_disassembly_num[gene_copy]), (int) 1);
  rates->pic_disassembly -= disassembly;
  rates->pic_disassembly_operations++;
  
  /* disassemble PIC in OFF state */
  if (state->active[gene_id][gene_copy] == OFF_PIC) {
    (state->active[gene_id][gene_copy]) = OFF_NO_PIC;
    state->state_change_ids[DEACETYLATION_STATE][gene_copy][rates->deacetylation_num[gene_copy]] = gene_id;
    (rates->deacetylation_num[gene_copy])++;
  }
  /* disassemble PIC in ON state */
  if (state->active[gene_id][gene_copy] == ON_FULL) {
    (state->active[gene_id][gene_copy]) = ON_NO_PIC;
  }
}

/*
 * update the activity state of 'gene_id' based on the new balance of
 * bound transcription factors and the current state of the gene
 */
void revise_activity_state(int gene_id,
                           int gene_copy,
                           Genotype *genotype,
                           CellState *state,
                           GillespieRates *rates)
{
  int transcriptrule=-1, oldstate=-1, numactive=-1;

  transcriptrule = ready_to_transcribe(gene_id, gene_copy, 
                                       state->tf_bound_indexes, 
                                       state->tf_bound_num,
                                       genotype->all_binding_sites,
                                       genotype->activating,
                                       &numactive);

  /* get last state of transcription initiation */
  oldstate = state->active[gene_id][gene_copy];
  //printf("gene id = %d\n", gene_id);
  //if(gene_id == 0){
            // system("PAUSE");
  //}
  /*
   * first set of rules:
   * ACTIVATING TFs exceed REPRESSING TFs 
   */
   LOG_ERROR("transcript rule = %d, oldstate = %d, transcriptrule & oldstate = %d, OFF_FULL = %d\n", transcriptrule, oldstate, (transcriptrule) && oldstate, OFF_FULL);
  /* OFF_FULL -> ON_WITH_NUCLEOSOME */
  if ((transcriptrule) && oldstate==OFF_FULL){
                       LOG_ERROR("UP ACE NUM!\n");
                       //printf("A\n");
    state->active[gene_id][gene_copy] = ON_WITH_NUCLEOSOME;
    state->state_change_ids[ACETYLATION_STATE][gene_copy][rates->acetylation_num[gene_copy]] = gene_id;
    (rates->acetylation_num[gene_copy])++; // rates will be updated in calc_dt
  }
  /* OFF_NO_PIC -> ON_NO_PIC */
  if ((transcriptrule) && oldstate==OFF_NO_PIC) {
                       //printf("B\n");
    state->active[gene_id][gene_copy] = ON_NO_PIC;
    remove_from_array(gene_id, DEACETYLATION_STATE,  state->state_change_ids[DEACETYLATION_STATE][gene_copy], 
                      &(rates->deacetylation_num[gene_copy]), (int) 1);
    if (numactive){
      state->state_change_ids[PICASSEMBLY_STATE][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id;
      (rates->pic_assembly_num[gene_copy])++;
      LOG_ERROR("PIC NUM = %d\n", rates->pic_assembly_num[gene_copy]);
    }
  }
  /* OFF_PIC -> ON_FULL */
  if ((transcriptrule) && oldstate==OFF_PIC) {
                      // printf("C\n");
    state->active[gene_id][gene_copy] = ON_FULL;
  }

  /*
   * second set of rules:
   * REPRESSING TFs exceed ACTIVATING TFs 
   */

  /* ON_WITH_NUCLEOSOME -> OFF_FULL */
  if (!(transcriptrule) && oldstate==ON_WITH_NUCLEOSOME) {
                       // printf("D\n");
    state->active[gene_id][gene_copy] = OFF_FULL;
    LOG_VERBOSE("removing gene=%d, copy=%d from state_change_ids[ACETYLATION_STATE][%d]\n", gene_id, gene_copy, gene_copy);
    remove_from_array(gene_id, ACETYLATION_STATE, state->state_change_ids[ACETYLATION_STATE][gene_copy], 
                      &(rates->acetylation_num[gene_copy]), (int) 1);
  }
  
  /* ON_NO_PIC -> OFF_NO_PIC */
  if (!(transcriptrule) && oldstate==ON_NO_PIC){    
                       // printf("E\n");      
    state->active[gene_id][gene_copy] = OFF_NO_PIC;
    remove_from_array(gene_id, PICASSEMBLY_STATE, state->state_change_ids[PICASSEMBLY_STATE][gene_copy], 
                      &(rates->pic_assembly_num[gene_copy]), (int) 0);
    state->state_change_ids[DEACETYLATION_STATE][gene_copy][rates->deacetylation_num[gene_copy]] = gene_id;
    (rates->deacetylation_num[gene_copy])++;
  }
  /* ON_FULL -> OFF_PIC  */
  if (!(transcriptrule) && oldstate==ON_FULL) {
                       // printf("F\n");
    state->active[gene_id][gene_copy] = OFF_PIC;
  }
//system("PAUSE");
  /* do remaining transitions:
   * OFF_PIC -> OFF_NO_PIC
   * ON_FULL -> ON_NO_PIC 
   */
  if ((state->active[gene_id][gene_copy]==OFF_PIC || state->active[gene_id][gene_copy]==ON_FULL) && numactive==0)
    disassemble_PIC(state, genotype, gene_id, gene_copy, rates);

  if (verbose && (oldstate!=state->active[gene_id][gene_copy])) {
    LOG_VERBOSE("state change from %d to %d in gene %d, copy %d\n", oldstate, state->active[gene_id][gene_copy], gene_id, gene_copy);
  }
}

/*
 * remove the bound transcription factor from the specified
 * site_id_to_unbind
 */
void remove_tf_binding(Genotype *genotype,
                       CellState *state,
                       GillespieRates *rates,
                       KonStates *kon_states,
                       int site_id_to_unbind,
                       float koffvalues[],
                       float t)
{
  int i, j, k, bound, site_id, gene_id, gene_copy;

  i = 0;

  /* given site 'site_id_to_unbind', look for the index in the list of bound sites */
  while ((state->tf_bound_indexes[i] != site_id_to_unbind) && (i < state->tf_bound_num)) 
    i++;
  if (i == state->tf_bound_num) {  /* couldn't find the site */
    LOG_ERROR("t=%g could not find site %d with %d possibilities\n Bound sites are\n",
              t, site_id_to_unbind, state->tf_bound_num);
    for (j = 0; j < state->tf_bound_num; j++)  {
      LOG_NOFUNC("%d\n", state->tf_bound_indexes[j]);
    }
  }
  else {
    j = 0;
    /* loop through the sterically hindered sites */
    while (j < state->tf_hindered_num) {

      /* check all sites hindered by binding to location 'site_id_to_unbind' */
      if (state->tf_hindered_indexes[j][1] == site_id_to_unbind) {
        k = bound = 0;

        /* is anything else hindering the same site? */
        while (bound == 0 && k < state->tf_hindered_num) {
          if (state->tf_hindered_indexes[j][0] == state->tf_hindered_indexes[k][0] && j != k) 
            bound=1;
          k++;
        }

        /* if nothing else is hindering this site then allow site_id_to_unbind to be (re-)bound */
        if (bound==0) {
          site_id = state->tf_hindered_indexes[j][0];
          LOG_VERBOSE("Site %d left_edge_pos %d on gene %d freed from steric hindrance\n",
                      site_id, genotype->all_binding_sites[site_id].left_edge_pos, genotype->all_binding_sites[site_id].cisreg_id);

          /* adjust rates by returning kon to pool */
          add_kon(state->protein_conc[genotype->all_binding_sites[site_id].tf_id],
                  kon_states->konvalues[genotype->all_binding_sites[site_id].tf_id][KON_SALPHC_INDEX],
                  genotype->all_binding_sites[site_id].tf_id,
                  site_id,
                  rates,
                  kon_states);
        }

        /* now we have one less sterically hindered site */
        (state->tf_hindered_num)--;

        if (j < state->tf_hindered_num) {
          /* shorten array by moving the end of the array to the hole opened up by removed site */
          state->tf_hindered_indexes[j][0] = state->tf_hindered_indexes[state->tf_hindered_num][0];
          state->tf_hindered_indexes[j][1] = state->tf_hindered_indexes[state->tf_hindered_num][1];
        }
      } else { /* only increment if we haven't shortened array */
        j++;
      }
    }    

    /* reduce the koff rate by the amount */
    rates->koff -= koffvalues[i];
    rates->koff_operations++;

    /* one less bound site */
    (state->tf_bound_num)--;

    /* shift end of array to hole opened up */
    state->tf_bound_indexes[i] = state->tf_bound_indexes[state->tf_bound_num];

    /* likewise with koffvalues */
    koffvalues[i] = koffvalues[state->tf_bound_num];

    /* find the gene and copy whose cisreg region has an unbinding event */
    gene_id = genotype->all_binding_sites[site_id_to_unbind].cisreg_id;
    gene_copy = genotype->all_binding_sites[site_id_to_unbind].gene_copy;
    LOG_VERBOSE("Add site %d at left_edge_pos %d on gene %d copy %d freed by unbinding\n",
                site_id_to_unbind, genotype->all_binding_sites[site_id_to_unbind].left_edge_pos, gene_id, gene_copy);

    /* adjust kon */
    add_kon(state->protein_conc[genotype->all_binding_sites[site_id_to_unbind].tf_id],
            kon_states->konvalues[genotype->all_binding_sites[site_id_to_unbind].tf_id][KON_SALPHC_INDEX],
            genotype->all_binding_sites[site_id_to_unbind].tf_id,
            site_id_to_unbind,
            rates,
            kon_states);

    /* adjust the state of the gene */
    revise_activity_state(gene_id, gene_copy, genotype, state, rates);

    /* when TF unbinds adjust the co-operativity at close sites */
    scan_nearby_sites(site_id_to_unbind, genotype->all_binding_sites, state, rates, koffvalues, t);
  }
}

/*
 * do transcription factor binding at specified site_id
 */
void attempt_tf_binding(Genotype *genotype,
                        CellState *state,
                        GillespieRates *rates,
                        float **koffvalues,
                        KonStates *kon_states,
                        int *maxbound2,
                        int *maxbound3,
                        int site_id,
                        float t)
{
  int gene_id, gene_copy, k, posdiff,i;

  LOG_VERBOSE("kon1 event at site %d out of %d possible, %d TFs previously bound binding_sites_num=%d\n",
              site_id, kon_states->nkon, state->tf_bound_num, genotype->binding_sites_num);

  /* if we have run out of space, double memory  */
  if (state->tf_bound_num >= *maxbound2){
    (*maxbound2) *= 2;
    state->tf_bound_indexes = realloc(state->tf_bound_indexes, (*maxbound2)*sizeof(int));
    
    /* do the copy */
    *koffvalues = realloc(*koffvalues, (*maxbound2)*sizeof(float));
    
    for(i=*maxbound2/2;i<*maxbound2;i++){
        state->tf_bound_indexes[i]=-1;
        (*koffvalues)[i]=0.0;
    }
    /* check return value */
    if (!state->tf_bound_indexes || !(*koffvalues)) {
      LOG_ERROR("memory allocation error resetting maxbound2=%d\n", *maxbound2);
      exit(1);
    }
  }

  /* append the site to end of indexes */
  state->tf_bound_indexes[state->tf_bound_num] = site_id;
  LOG_VERBOSE("remove site %3d on gene %2d (copy %d)\n", 
              site_id, genotype->all_binding_sites[site_id].cisreg_id, genotype->all_binding_sites[site_id].gene_copy);

  /* remove the site_id from the kon pool */
  remove_kon(site_id,
             genotype->all_binding_sites[site_id].tf_id,
             rates, 
             kon_states->konvalues[genotype->all_binding_sites[site_id].tf_id][KON_SALPHC_INDEX],
             kon_states,
             state->protein_conc[genotype->all_binding_sites[site_id].tf_id]);

  /* recompute the koffvalues */
  calc_koff(site_id, genotype->all_binding_sites, state, &((*koffvalues)[state->tf_bound_num]), t);

  LOG_VERBOSE("new koff = %g is number %d\n",
              (*koffvalues)[state->tf_bound_num], (state->tf_bound_num+1));

  /* adjust rates by adding the new koffvalue to rates->koff */
  rates->koff += (*koffvalues)[state->tf_bound_num];
  rates->koff_operations++;
  /* append site_id to list of bound sites */
  state->tf_bound_indexes[state->tf_bound_num] = site_id;
  /* increment number of bound TFs */
  (state->tf_bound_num)++;
  /* adjust co-operative binding in context of new TF */
  scan_nearby_sites(site_id, genotype->all_binding_sites, state, rates, *koffvalues, t);
  /* get the gene that the TF is binding to */
  gene_id = genotype->all_binding_sites[site_id].cisreg_id;
  /* get the copy that the TF is binding to */
  gene_copy = genotype->all_binding_sites[site_id].gene_copy;
  
  /* update steric hindrance data structures */
  /* JM: this cycles over all sites, not just bound ones, in order to
     record redundancy in steric hindrance*/
  for (k = 0; k < genotype->binding_sites_num; k++) {
    /* if we are on the same gene and not the same binding site */
    if (gene_id == genotype->all_binding_sites[k].cisreg_id &&
        gene_copy == genotype->all_binding_sites[k].gene_copy &&
        !(k==site_id)) {

      /* check distance from current binding site (k) to the original (site_id) */
      LOG_VERBOSE("site_id=%d k=%d\n", genotype->all_binding_sites[site_id].left_edge_pos, genotype->all_binding_sites[k].left_edge_pos);
      posdiff = genotype->all_binding_sites[site_id].left_edge_pos - genotype->all_binding_sites[k].left_edge_pos;

      /* if within HIND_LENGTH, we prevent future binding by adding to steric hindrance */
      if (abs(posdiff) < HIND_LENGTH) {
        /* if not enough memory, reallocate */
        if (state->tf_hindered_num >= *maxbound3) {
          (*maxbound3) *= 2;
          state->tf_hindered_indexes = realloc(state->tf_hindered_indexes,2*(*maxbound3)*sizeof(int));
         
          for(i=*maxbound3/2;i<*maxbound3;i++){
            state->tf_hindered_indexes[1][i]=-1;
            state->tf_hindered_indexes[0][i]=-1;
          }
        }
        /* record hindrance: 'site_id' blocks 'k' */
        state->tf_hindered_indexes[state->tf_hindered_num][0] = k;
        state->tf_hindered_indexes[state->tf_hindered_num][1] = site_id;
        /* update list of hindered count */
        (state->tf_hindered_num)++;
        LOG_VERBOSE("%d steric hindrance sites after %d blocks site %d\n", state->tf_hindered_num, site_id, k);

        /* remove the kon from pool */
        remove_kon(k,
                   genotype->all_binding_sites[k].tf_id,
                   rates,
                   kon_states->konvalues[genotype->all_binding_sites[k].tf_id][KON_SALPHC_INDEX],
                   kon_states,
                   state->protein_conc[genotype->all_binding_sites[k].tf_id]);
      }
    }
  }
  LOG_VERBOSE("tf_bound_num=%d tf_hindered_num=%d maxbound2=%d maxbound3=%d\n",
              state->tf_bound_num, state->tf_hindered_num, *maxbound2, *maxbound3);

  /* gene activity may change as a result of binding */
  revise_activity_state(gene_id, gene_copy, genotype, state, rates);
}

/*
 * time course of [TF]s represented as array of TimeCourse lists.
 */
void add_time_points(float time,
                     float protein_conc[NPROTEINS],
                     TimeCourse **timecoursestart,
                     TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i < NPROTEINS; i++)
    add_time_point(time, protein_conc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

void add_integer_time_points(float time,
                             int protein_conc[NPROTEINS],
                             TimeCourse **timecoursestart,
                             TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i < NPROTEINS; i++)
    add_time_point(time, (float) protein_conc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

/*
 * compute tprime factor used in the integration of growth rate
 */
float compute_tprime(float c, 
                     float P, 
                     float alpha, 
                     float s_mRNA) 
{
  return (1/c) * log((c*P - alpha*s_mRNA)/(c*P - alpha*s_mRNA));
}

/*
 * get integral for growth rate
 */
float compute_integral(float alpha, 
                       float c, 
                       float gmax, 
                       float deltat, 
                       float s_mRNA, 
                       float P, 
                       float Pp,                       
                       float ect1) 
{
  return 1.0/(pow(c,2)*Pp) * gmax * (-alpha*ect1*s_mRNA + c*(P*ect1 + alpha*deltat*s_mRNA));
}

/*
 * return the instantaneous growth rate given the current cell state,
 * also return the integrated growth rate as a pointer
 */
float compute_growth_rate_dimer(float *integrated_growth_rate,
                                float alpha_a, 
                                float s_mRNA_a,
                                float alpha_b, 
                                float s_mRNA_b,
                                float all_alpha[NGENES],
                                int all_s_mRNA[NGENES],
                                float P_a,
                                float P_b,
                                float P_next_a,
                                float P_next_b,
                                float t, 
                                float deltat,
                                float c_a,                                
                                float ect1_a,
                                float c_b,                                
                                float ect1_b,
				int env) 
{
  int i;
  float instantaneous_growth_rate;  /* this is returned from the function */
  float total_alpha_s = 0.0;
  float deltatprime, deltatrest;



//  LOG_VERBOSE_NOCELLID("P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//                       P, P_next, c, t, t+deltat, s_mRNA);
                       
  switch (env){ // indicator of environment
  
	  case 1:{ // protein A is necessary!
	  
			 /* choose the appropriate piecewise linear integral */
		  if (((P_a >= Pp_a) && (P_next_a >= P_a)) || ((P_next_a >= Pp_a) && (P_a >= P_next_a))) {          /* P > Pp throughout */
//		    if (verbose)
//		      LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);
		
		    *integrated_growth_rate = gmax_a * deltat;	    
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);	     
		  } 
		  
		  else if (((P_next_a <= Pp_a) && (P_next_a >= P_a)) || ((P_a <= Pp_a) && (P_a >= P_next_a))) {   /* P < Pp throughout */
//		    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_a, c_a, gmax_a, deltat, s_mRNA_a, P_a, Pp_a, ect1_a);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);

		  } 
		  
		  
		  else if ((Pp_a > P_a) && (P_next_a > Pp_a)) {    /* P < Pp up until t' then P > Pp */
		    deltatprime = compute_tprime(c_a, P_a, alpha_a, s_mRNA_a);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_a, c_a, gmax_a, deltatprime, s_mRNA_a, P_a, Pp_a, ect1_a);
		    *integrated_growth_rate += gmax_a * deltatrest;
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);

		  } 
		  
		  
		  else if ((P_a > Pp_a) && (Pp_a > P_next_a)) {   /* P > Pp up until t' then P < Pp */
		    deltatprime = compute_tprime(c_a, P_a, alpha_a, s_mRNA_a);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = gmax_a * deltatprime;
		    *integrated_growth_rate += compute_integral(alpha_a, c_a, gmax_a, deltatrest, s_mRNA_a, P_a, Pp_a, ect1_a);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);

		  } 
		  
		  
		  else {
//		    LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//		                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
//		    LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
//		                       cell_id);		
			
		    exit(1);
		  }
		
		  /* compute instantaneous growth rate at t */
		  if (P_next_a < Pp_a)
		    instantaneous_growth_rate = gmax_a*P_next_a/Pp_a - penalty*P_next_b;
		  else
		    instantaneous_growth_rate = gmax_a - penalty*P_next_b;
			
			break;
	}
	case 0:{ // protein b is necessary!
		 /* choose the appropriate piecewise linear integral */
		  if (((P_b >= Pp_b) && (P_next_b >= P_b)) || ((P_next_b >= Pp_b) && (P_b >= P_next_b))) {          /* P > Pp throughout */
//		    if (verbose)
//		      LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);
		    *integrated_growth_rate = gmax_b * deltat;
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  else if (((P_next_b <= Pp_b) && (P_next_b >= P_b)) || ((P_b <= Pp_b) && (P_b >= P_next_b))) {   /* P < Pp throughout */
//		    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_b, c_b, gmax_b, deltat, s_mRNA_b, P_b, Pp_b, ect1_b);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  
		  else if ((Pp_b > P_b) && (P_next_b > Pp_b)) {    /* P < Pp up until t' then P > Pp */
		    deltatprime = compute_tprime(c_b, P_b, alpha_b, s_mRNA_b);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_b, c_b, gmax_b, deltatprime, s_mRNA_b, P_b, Pp_b, ect1_b);
		    *integrated_growth_rate += gmax_b * deltatrest;
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  
		  else if ((P_b > Pp_b) && (Pp_b > P_next_b)) {   /* P > Pp up until t' then P < Pp */
		    deltatprime = compute_tprime(c_b, P_b, alpha_b, s_mRNA_b);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = gmax_b * deltatprime;
		    *integrated_growth_rate += compute_integral(alpha_b, c_b, gmax_b, deltatrest, s_mRNA_b, P_b, Pp_b, ect1_b);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  
		  else {
//		    LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//		                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
//		    LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
//		                       cell_id);
		    exit(1);
		  }
		
		  /* compute instantaneous growth rate at t */
		  if (P_next_b < Pp_b)
		    instantaneous_growth_rate = gmax_b*P_next_b/Pp_b - penalty*P_next_a;
		  else
		    instantaneous_growth_rate = gmax_b - penalty*P_next_a;
		
		  break;
	}
  }
//  /* choose the appropriate piecewise linear integral */
//  if (((P > Pp) && (P_next >= P)) || ((P_next > Pp) && (P >= P_next))) {          /* P > Pp throughout */
//    if (verbose)
//      LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);
//    *integrated_growth_rate = gmax * deltat;
//  } 
//  
//  else if (((P_next < Pp) && (P_next >= P)) || ((P < Pp) && (P >= P_next))) {   /* P < Pp throughout */
//    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
//    *integrated_growth_rate = compute_integral(alpha, c, gmax, deltat, s_mRNA, P, Pp, ect, ect1);
//  } 
//  
//  
//  else if ((Pp > P) && (P_next > P)) {    /* P < Pp up until t' then P > Pp */
//    deltatprime = compute_tprime(c, P, alpha, s_mRNA);
//    deltatrest = deltat - deltatprime;
//    LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
//                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
//    *integrated_growth_rate = compute_integral(alpha, c, gmax, deltatprime, s_mRNA, P, Pp, ect, ect1);
//    *integrated_growth_rate += gmax * deltatrest;
//  } 
//  
//  
//  else if ((P > Pp) && (P > P_next)) {   /* P > Pp up until t' then P < Pp */
//    deltatprime = compute_tprime(c, P, alpha, s_mRNA);
//    deltatrest = deltat - deltatprime;
//    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
//                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
//    *integrated_growth_rate = gmax * deltatprime;
//    *integrated_growth_rate += compute_integral(alpha, c, gmax, deltatrest, s_mRNA, P, Pp, ect, ect1);
//  } 
//  
//  
//  else {
//    LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
//    LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
//                       cell_id);
//    exit(1);
//  }
//
//  /* compute instantaneous growth rate at t */
//  if (P < Pp)
//    instantaneous_growth_rate = gmax*P/Pp;
//  else
//    instantaneous_growth_rate = gmax;

  LOG_VERBOSE_NOCELLID("growth rate (variable %g)-", *integrated_growth_rate);

  /* compute the total cost of translation across all genes  */
  for (i=0; i < NGENES; i++) {
    total_alpha_s += all_alpha[i] * all_s_mRNA[i];
  }

  /* add constant term for integrated rate */
  *integrated_growth_rate += -h * deltat * (total_alpha_s);

  /* and instantaneous integrated rate */
  instantaneous_growth_rate += -h * (total_alpha_s);

  LOG_VERBOSE_NOFUNC("(constant %g) = (total %g)\n", (h*deltat*total_alpha_s), *integrated_growth_rate);

  /* make sure growth rates can't be negative */
  if (*integrated_growth_rate < 0.0)
    *integrated_growth_rate = 0.0;

  if (instantaneous_growth_rate < 0.0)
    instantaneous_growth_rate = 0.0;
  
#if 0 /* currently disable printing out growth rate information */
  fprintf(fp_growthrate[0], "%g %g %g %g %g %g\n", t, instantaneous_growth_rate, *integrated_growth_rate, P, s_mRNA, c);
#endif


  /* return the instantaneous growth rate */
  return (instantaneous_growth_rate);
}

/* 
 * update both the protein concentration and current cell size
 *
 * need some sort of control in case it declines to essentially zero.
 * Add in discrete, stochastic and/or zero values, but this may create
 * false attractor if time steps are small and rising tide keeps
 * getting rounded down
 */
void update_protein_conc_cell_size(float protein_conc[],
                                   CellState *state,
                                   Genotype *genotype,
                                   float dt,
                                   GillespieRates *rates,
                                   KonStates *kon_states,
                                   float t,
//                                   TimeCourse **timecoursestart,
//                                   TimeCourse **timecourselast,
                                   float otherdata[],
				   int *env)
{
  int i;
  float ct, ect, ect1,ect_a,ect1_a,ect_b,ect1_b;
  float L_a, L_b;
  float instantaneous_growth_rate = 0.0;
  float integrated_growth_rate = 0.0;
  float adjusted_decay;


 
  rates->max_salphc = rates->min_salphc = 0.0;
  for (i=0; i < NPROTEINS; i++) {
//   printf("p%i=%f\n",i,protein_conc[i]);
    if (i == SELECTION_GENE_A)  /* if we are looking at the selection gene, record protein concentration before update */
      L_a = protein_conc[i];
    if (i == SELECTION_GENE_B)
	  L_b = protein_conc[i];
	
    /* update protein decay rates due to dilution caused by growth */
    adjusted_decay = genotype->proteindecay[i] + state->growth_rate;
 
    /* if this results in a very small or zero decay rate, use protein aging term */
    if (adjusted_decay > protein_aging)
      kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX] = adjusted_decay;
    else 
      kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX] = protein_aging;

    /* print out warning if decay rates get too low */
    if (kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX] < 1e-10) {
      LOG_WARNING("protein=%02d, protein_decay=%g, genotype->proteindecay=%g, protein_aging=%g\n", i, adjusted_decay, 
                  genotype->proteindecay[i], protein_aging);
    }

    LOG_VERBOSE("prot decay[%d]=%g\n", i, kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX]);

    ct = kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX]*dt;
    ect = exp(-ct);
    if (fabs(ct)<EPSILON) ect1=ct;
    else ect1 = 1-ect;  
	
	if (i == SELECTION_GENE_A) 
	{
		ect_a=ect;
		ect1_a=ect1;		
	}
	if (i == SELECTION_GENE_B)
	{
		ect_b=ect;
		ect1_b=ect1;		
	}
	
    /* get the new protein concentration for this gene */
    protein_conc[i] = kon_states->konvalues[i][KON_SALPHC_INDEX]*ect1 + ect*protein_conc[i];

    /* recompute the konvalues and max and min salphc */
    kon_states->konvalues[i][KON_DIFF_INDEX] = (protein_conc[i] - kon_states->konvalues[i][KON_SALPHC_INDEX]) / kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX];
    rates->max_salphc += ((float) kon_states->kon_list[i]->site_count) * fmaxf(protein_conc[i], kon_states->konvalues[i][KON_SALPHC_INDEX]);
    rates->min_salphc += ((float) kon_states->kon_list[i]->site_count) * fminf(protein_conc[i], kon_states->konvalues[i][KON_SALPHC_INDEX]);
   
  }
 
// printf("La=%f,Lb=%f\n",L_a,L_b);

  if (*env==1) {  /* now find out the protein concentration at end of dt interval and compute growth rate */ // a is the required
      
      instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genotype->translation[SELECTION_GENE_A], state->mRNA_cyto_num[SELECTION_GENE_A], 
                                                            genotype->translation[SELECTION_GENE_B], state->mRNA_cyto_num[SELECTION_GENE_B],
                                                            genotype->translation, state->mRNA_cyto_num, 
                                                            L_a, L_b, protein_conc[SELECTION_GENE_A], protein_conc[SELECTION_GENE_B],t, dt, 
                                                            kon_states->konvalues[SELECTION_GENE_A][KON_PROTEIN_DECAY_INDEX], ect1_a,
                                                            kon_states->konvalues[SELECTION_GENE_B][KON_PROTEIN_DECAY_INDEX], ect1_b,1);
      /* us the integrated growth rate to compute the cell size in the next timestep */
       
      state->cell_size = (state->cell_size)*exp(integrated_growth_rate);
      
      fprintf(fp_cellsize[0], "%g %g\n", t, state->cell_size);
    }
    
    if (*env==0) {  /* now find out the protein concentration at end of dt interval and compute growth rate */ // b is the required
         
      instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genotype->translation[SELECTION_GENE_A], state->mRNA_cyto_num[SELECTION_GENE_A], 
                                                            genotype->translation[SELECTION_GENE_B], state->mRNA_cyto_num[SELECTION_GENE_B], 
                                                            genotype->translation, state->mRNA_cyto_num, 
                                                            L_a, L_b, protein_conc[SELECTION_GENE_A],protein_conc[SELECTION_GENE_B], t, dt, 
                                                            kon_states->konvalues[SELECTION_GENE_A][KON_PROTEIN_DECAY_INDEX], ect1_a,
                                                            kon_states->konvalues[SELECTION_GENE_B][KON_PROTEIN_DECAY_INDEX], ect1_b,0);
      /* us the integrated growth rate to compute the cell size in the next timestep */
       
      state->cell_size = (state->cell_size)*exp(integrated_growth_rate);
      
      fprintf(fp_cellsize[0], "%g %g\n", t, state->cell_size);
    }
// 	printf("dt = %f\n", dt);
// 	printf(" integrated GR = %f\n", integrated_growth_rate);
// 

 
 
 
  /* scale up the rates using kon global */
  rates->max_salphc *= kon;
  rates->min_salphc *= kon;
//  if ((output) && (*timecourselast)->time < t+dt-0.1) 
//    add_time_points(t+dt, otherdata, timecoursestart, timecourselast);

  /* update the instantaneous growth rate for the beginning of the next timestep */
  state->growth_rate = instantaneous_growth_rate;
  
}

/*
 * compute the number of transcription factors expected to be bound
 */
void calc_num_bound(float protein_conc[],
                    int tf_bound_num)
{
  float sum;
  int i;
  
  sum = 0.0;
  for (i=0; i < TFGENES; i++) 
    sum += protein_conc[i];
  /* if this is wrong for random sequences adjust Kr accordingly */
  /* fprintf(fperrors, "%d bound %g expected\n", tf_bound_num, 0.0003*sum);*/

  LOG_VERBOSE_NOCELLID("%d bound %g expected\n", tf_bound_num, (CISREG_LEN*TFGENES*sum)/NUMSITESINGENOME);
}

/*
 * sums the rate across all copies of the gene
 */
int sum_rate_counts(int rate_array[MAX_COPIES])
{
  int i;
  float retval = 0.0;

  for (i = 0; i < MAX_COPIES; i++) {
    retval += rate_array[i];
  }
  return retval;
}

/* 
 * given a rate array and position within the array, get the the gene
 * copy and location
 */
void get_gene(int rate_array[MAX_COPIES], 
              int pos, 
              int *gene_loc, 
              int *gene_copy)
{
  int i = 0;
  int total_rate = 0;
  *gene_copy = -1;   /* haven't found the copy yet */

  while (i < MAX_COPIES && *gene_copy < 0) {
    if (pos < (total_rate + rate_array[i])) {
      *gene_copy = i;
      *gene_loc = pos - total_rate;
    } else {
      total_rate += rate_array[i];
      i++;
    }
  }
}


/*
 * START
 * Functions that handle each possible Gillespie event 
 *
 */
void transport_event(GillespieRates *rates,
                     CellState *state,
                     Genotype *genotype,
                     KonStates *kon_states,
                     float transport_rate[NGENES],
//                     TimeCourse **timecoursestart, 
//                     TimeCourse **timecourselast, 
                     float dt,
                     float t,
                     float x,
					 int *env)
{
  int i;
  float temp_rate;
  float endtime = t + dt + TTRANSLATION; 

//for(i=0;i<NGENES;i++)
//{
//	printf("rate%d=%f the address is %i\n",i,transport_rate[i],&transport_rate[i]);
//	
//}
  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, 
                                rates, kon_states, 
                                t, 
//								timecoursestart, timecourselast, 
                                state->protein_conc,env);
  
  i = -1;
  temp_rate = 0.0;  
 
  /* choose gene product (mRNA) that gets transported to cytoplasm
     based on weighting in transport[] array */
  while (i < NGENES && x > temp_rate) {
    i++;
    x -= transport_rate[i];
//    printf("mRNA%d=%d, rate=%f\n",i,state->mRNA_nuclear_num[i],transport_rate[i]);
//    printf("%f\n",transport_rate[i]);
  }
 
  if (i >= NGENES) {
    LOG_ERROR("[cell %03d] attempted to choose mRNA for gene=%d which doesn't exist\n", state->cell_id, i);
    exit(0);
  } 

  LOG_VERBOSE("do transport event mRNA from gene=%d from %d copies (x=%g)\n", i, state->mRNA_nuclear_num[i], x);


  (state->mRNA_nuclear_num[i])--;   /* one less mRNAs in nucleus */
  (state->mRNA_transl_cyto_num[i])++;   /* it has just arrived in cytoplasm, ready to be translated */

  /* add the endtime for translation */
  LOG_VERBOSE("add translation event endtime=%f for mRNA encoded by gene=%d \n", endtime, i);
  add_fixed_event_end(i, -1, endtime, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));

  transport_rate[i] -= KRNA;   /* decrease transport frequency */

  /* if below a threshold, make zero */
  if (transport_rate[i] < 0.1*KRNA) 
    transport_rate[i]=0.0;

  /* adjust rates */
  rates->transport -= KRNA;
  rates->transport_operations++;

  /* do similar threshold check */
  if (rates->transport < 0.1*KRNA) 
    rates->transport=0.0;
}


void tf_binding_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                      KonStates *kon_states, float *koffvalues, 
//					  TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float konrate, float dt, float t, int maxbound2, int maxbound3, int update_protein_flag,int *env)
{
  float x = ran1(&seed) * (rates->salphc + konrate)/kon;
  int k, l = -1;  /* new */
  float total_konrate, konrate_for_TF = 0.0;     
  int site_id = -1;
  int tracker=0;

  /* loop through all TFs, then choose a particular binding site */
  for (k=0; k < TFGENES; k++) {
 
    /* if no sites available for this TF, skip to next gene */
    if (kon_states->kon_list[k]->site_count == 0) {
      LOG_VERBOSE("looking at TF: %d, has %d sites available, skipping\n", k, kon_states->kon_list[k]->site_count);
      continue;
    }
  
    /* compute the total rate for all available binding sites for this
     * particular TF: see if we are in the right range
     */

    /* TODO: commented-out code that may help fix numerical issues with 1-exp(), but needs further testing */
    /* float c = kon_states->konvalues[k][KON_PROTEIN_DECAY_INDEX];
     float ectdt;
     if (fabs(c*dt)<EPSILON) ectdt=c;
     else ectdt = (1-exp(-c*dt))/dt;    
     konrate_for_TF = kon_states->konvalues[k][KON_SALPHC_INDEX] + kon_states->konvalues[k][KON_DIFF_INDEX] * ectdt;  */

    /* first, cache the konrate for this particular gene */
    konrate_for_TF = kon_states->konvalues[k][KON_SALPHC_INDEX] + 
       kon_states->konvalues[k][KON_DIFF_INDEX] * (1-exp(-kon_states->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt))/dt;  
 
    LOG_VERBOSE("TF:%d [KON_SALPHC: %g, KON_DIFF: %g, KON_PROTEIN_DECAY: %g]\nTF:%d [1-exp(-ct): %g, (1-exp(-ct)/dt): %g, konrate_for_TF: %g]\n", 
                k,
                kon_states->konvalues[k][KON_SALPHC_INDEX],
                kon_states->konvalues[k][KON_DIFF_INDEX],
                kon_states->konvalues[k][KON_PROTEIN_DECAY_INDEX],
                k,
                1-exp(-kon_states->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt),
                (1-exp(-kon_states->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt))/dt,
                konrate_for_TF); 

    /* compute the *total* kon rate for all unbound sites for this TF  */
    total_konrate = ((kon_states->kon_list[k]->site_count)) * konrate_for_TF;

//    LOG_VERBOSE("looking at TF: %d, has %d sites available [konrate: %g, total_konrate: %g, x: %g]\n", 
//                k, kon_states->kon_list[k]->site_count, konrate_for_TF, total_konrate, x); 

    /* if we are already in the appropriate TF, now choose a binding site amongst this set */
    if (!(x > total_konrate) || (k == TFGENES - 1)) {
        	
      float konrate = 0.0;
      LOG_VERBOSE("selecting TF: %d, konrate: %g, total_konrate: %g, x: %g\n", k, konrate_for_TF, total_konrate, x); 
      
      /* look through the list of available sites using index 'l' */
      while (l < (kon_states->kon_list[k]->site_count - 1) && x > konrate) {
      	
        l++;
        site_id = kon_states->kon_list[k]->available_sites[l];         /* get ID of site */
        
        LOG_VERBOSE("l: %d, site: %d, binds to TF: %d, x = %g (site_count=%d)\n", l, site_id, k, x, kon_states->kon_list[k]->site_count); 
        konrate = konrate_for_TF;
               
        x -= konrate;         /* adjust random number */
       
      }
      /* found it, so break out of the outer for loop */
      
      break;
    } else {
    	x -= total_konrate; 
    }
     
  }
 
  /* print error if site not found */
  if (site_id == -1) {
    LOG_ERROR("no binding site could be found  TF: total_konrate: %g, x: %g\n", total_konrate, x);
  }
  else {
    LOG_VERBOSE("found a binding site l: %d, site: %d, binds to TF: %d, konrate: %g, x: %g\n", l, site_id, k, konrate_for_TF, x);  
  }

  if (update_protein_flag) {
    /* update protein concentration before doing the binding */
    
    update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, 
                                  rates, kon_states, 
                                  t,
//								   timecoursestart, timecourselast,
                                  state->protein_conc,env);
  }
   
  /* bind site_id, only if found */
  if (site_id != -1) {
    attempt_tf_binding(genotype, state, rates, &koffvalues, kon_states, &maxbound2, &maxbound3, site_id, t);
  }

  /* calculate the number of TFs bound */
  calc_num_bound(state->protein_conc, state->tf_bound_num);
}

void tf_unbinding_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                       KonStates *kon_states, float *koffvalues, 
//					   TimeCourse **timecoursestart, TimeCourse **timecourselast,
                       float konrate, float dt, float t, float x, int update_protein_flag, int *ignore_event,int *env)
{
  int i, j = -1;
  int site;
  
  *ignore_event = 0;
  
  /* if an attempt is made to do an unbinding event when there is
     nothing bound then we recompute kon and koff */
  if (state->tf_bound_num <= 0) {
    LOG_ERROR("tf_bound_num=%d, t=%g, dt=%g, x=%g, rates->koff=%g\n", state->tf_bound_num, t, dt, x, rates->koff);

    for (i=0; i < NGENES; i++) {
      int k;
      LOG_ERROR("[gene %2d]: mRNA nuclear=%d, cyto=%d, transl=%d, transcribing ", 
                i, state->mRNA_nuclear_num[i], state->mRNA_cyto_num[i], state->mRNA_transl_cyto_num[i]);
      for (k=0; k < MAX_COPIES; k++)
        LOG_NOFUNC("copy%d=%d ", k, state->mRNA_transcr_num[i][k]);
      LOG_NOFUNC("\n");
    }
    LOG_ERROR("attempting to unbind when nothing bound, recomputing koff and kon rates\n");
    recompute_koff_rates(rates, state, genotype, koffvalues, t);    
    recompute_kon_rates(rates, state, genotype, kon_states, 0);
    *ignore_event = 1; /* ignore this event, don't advance t */
    return;
  }

  while (j < state->tf_bound_num && x > 0) {
    j++;
    x -= koffvalues[j];
  }
  if (j==state->tf_bound_num) {
    float temp_rate = 0.0;
    for (i = 0; i < state->tf_bound_num; i++) 
      temp_rate += koffvalues[i];
    LOG_WARNING("t=%g koffvalues add up to %g instead of rates->koff=%g\n",
                t, temp_rate, rates->koff);
    rates->koff = temp_rate;
    j--; /* a bit of a fudge for rounding error, really should move on to rates->transport, but too complicated for something so minor */
  } 
  site = state->tf_bound_indexes[j];
  LOG_VERBOSE("t=%g koff event %d of %d at site %d\n", t, j, state->tf_bound_num,site);
  if (j < 0) { LOG_ERROR("t=%g (x=%g) koff event %d of %d at site %d\n", t, x, j, state->tf_bound_num,site); exit(-1); }
  
  if (update_protein_flag) 
    /* update protein concentration before removing TF */
    update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, 
                                  rates, kon_states, 
                                  t, 
//								  timecoursestart, timecourselast, 
                                  state->protein_conc,env);
  
  /* remove TF binding from 'site' */
  remove_tf_binding(genotype, state, rates, kon_states, site, koffvalues, t);
  calc_num_bound(state->protein_conc, state->tf_bound_num);
}

void mRNA_decay_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                      KonStates *kon_states, float *mRNAdecay, 
//					  TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t, float x, int *env)
{
  int i = -1, j;
  float temp_rate = 0.0;

  /* loop through mRNA products, to choose the mRNA with the
     proportionally higher decay rate */
  while (i < NGENES-1 && x > temp_rate) {
    i++;
    temp_rate += mRNAdecay[i];
  }
  if (x > temp_rate) { /* JM: had some rounding errors with rates->mRNAdecay. Calculate in calc_dt, hopefully fixed now */
    LOG_WARNING("x=%g > temp_rate=%g out of rates->mRNAdecay=%g\n",
                x, temp_rate, rates->mRNAdecay);
  }

  /* assume mRNA cytoplasm transport events equally likely */
  x = ran1(&seed)*((float) (state->mRNA_cyto_num[i] + state->mRNA_transl_cyto_num[i]));

  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt,
                                rates, kon_states, 
                                t, 
//								timecoursestart, timecourselast,
                                state->protein_conc, env);
  /* decay mRNA in cytoplasm */
  if (x < (float)state->mRNA_cyto_num[i]) {
    LOG_VERBOSE("mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
                i, state->mRNA_cyto_num[i], state->mRNA_transl_cyto_num[i]);
    
    /* remove the mRNA from the cytoplasm count */
    (state->mRNA_cyto_num[i])--;  
    change_mRNA_cytoplasm(i, genotype, state, rates, kon_states); 
    
  } else {
    /* decay mRNA in process of translating */
    x = ran1(&seed)*((float) state->mRNA_transl_cyto_num[i]);
    LOG_VERBOSE("mRNA decay event gene %d not from %d copies in cytoplasm but %f from %d copies translating\n",
                i, state->mRNA_cyto_num[i], trunc(x), state->mRNA_transl_cyto_num[i]);
    
    /* delete this fixed event: this mRNA will never be translated */
    LOG_VERBOSE("delete fixed TRANSLATION EVENT at time =%d for gene=%d\n", (int) trunc(x), i);
    delete_fixed_event(i, -1, (int) trunc(x), &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
    
    /* remove the mRNA from the count */
    (state->mRNA_transl_cyto_num[i])--;
    if (verbose) 
      for (j=0; j < NGENES; j++) {
        LOG_VERBOSE("%d copies of gene %d translating\n", state->mRNA_transl_cyto_num[j], j);
      }
  }
}

void histone_acteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                               KonStates *kon_states, 
//							   TimeCourse **timecoursestart, TimeCourse **timecourselast,
                               float dt, float t, int *env)
{
  int gene_loc, gene_copy;
  float x = ran1(&seed)*((float) sum_rate_counts(rates->acetylation_num));

  /* choose a particular gene to change state */
  get_gene(rates->acetylation_num, (int)trunc(x), &gene_loc, &gene_copy);
  int gene_id = state->state_change_ids[ACETYLATION_STATE][gene_copy][gene_loc];

  LOG_VERBOSE("acetylation event gene %d (copy %d)\nstate change from %d to 4\n",
              gene_id, gene_copy, state->active[gene_id][gene_copy]);
  if (state->active[gene_id][gene_copy] != ON_WITH_NUCLEOSOME) {
    LOG_ERROR("acetylation event on gene %d (copy %d) attempted from state %d\n", gene_loc, gene_copy, state->active[gene_id][gene_copy]);
  }

  /* update protein concentration and cell size */
  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt,
                                rates, kon_states, 
                                t,
//								 timecoursestart, timecourselast,
                                state->protein_conc,env);
  
  /* set state: eject nucleosome, but there is no PIC yet */
  state->active[gene_id][gene_copy] = ON_NO_PIC;
  remove_from_array(gene_id, ACETYLATION_STATE, state->state_change_ids[ACETYLATION_STATE][gene_copy], &(rates->acetylation_num[gene_copy]), (int) 1);
  if (is_one_activator(gene_id, gene_copy, state->tf_bound_indexes, state->tf_bound_num, 
                       genotype->all_binding_sites, genotype->activating)) {
    state->state_change_ids[PICASSEMBLY_STATE][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id; 
    (rates->pic_assembly_num[gene_copy])++;
  }
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                                 KonStates *kon_states, 
//								 TimeCourse **timecoursestart, TimeCourse **timecourselast,
                                 float dt, float t, int *env)
{
  int gene_copy, gene_loc; 
  float x = ran1(&seed)*((float) sum_rate_counts(rates->deacetylation_num));

  /* choose a particular gene and copy to change state */
  get_gene(rates->deacetylation_num, (int)trunc(x), &gene_loc, &gene_copy);
  int gene_id = state->state_change_ids[DEACETYLATION_STATE][gene_copy][gene_loc];

  LOG_VERBOSE("deacetylation event gene %d (copy %d)\nstate change from %d to 1\n",
              gene_id, gene_copy, state->active[gene_id][gene_copy]);
  if (state->active[gene_id][gene_copy] != OFF_NO_PIC) {
    LOG_ERROR("deacetylation event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, 
                                rates, kon_states, 
                                t,
//								 timecoursestart, timecourselast,
                                state->protein_conc,env);
  /* set state: nucleosome returns */
  state->active[gene_id][gene_copy] = OFF_FULL;
  remove_from_array(gene_id, DEACETYLATION_STATE, state->state_change_ids[DEACETYLATION_STATE][gene_copy], &(rates->deacetylation_num[gene_copy]), (int) 1);
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                      KonStates *kon_states, 
//					  TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t, int *env)
{
  float x = ran1(&seed)*((float) sum_rate_counts(rates->pic_assembly_num));

  int gene_copy; 
  int gene_loc; 

  /* choose a particular gene and copy to change state */
  get_gene(rates->pic_assembly_num, (int)trunc(x), &gene_loc, &gene_copy);
  int gene_id = state->state_change_ids[PICASSEMBLY_STATE][gene_copy][gene_loc];

  LOG_VERBOSE("PIC assembly event gene %d copy %d\nstate change from %d to 6\n",
              gene_id, gene_copy, state->active[gene_id][gene_copy]);

  if (state->active[gene_id][gene_copy] != ON_NO_PIC) {
    LOG_ERROR("PIC assembly event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt,
                                rates, kon_states, 
                                t, 
//								timecoursestart, timecourselast,
                                state->protein_conc, env);
  
  /* turn gene fully on: ready for transcription and adjust rates */
  state->active[gene_id][gene_copy] = ON_FULL;
  remove_from_array(gene_id, PICASSEMBLY_STATE, state->state_change_ids[PICASSEMBLY_STATE][gene_copy], 
                    &(rates->pic_assembly_num[gene_copy]), (int) 1);
  state->state_change_ids[TRANSCRIPTINIT_STATE][gene_copy][rates->transcript_init_num[gene_copy]] = gene_id;
  (rates->transcript_init_num[gene_copy])++;
  state->state_change_ids[PICDISASSEMBLY_STATE][gene_copy][rates->pic_disassembly_num[gene_copy]] = gene_id;

  (rates->pic_disassembly_num[gene_copy])++;

  rates->pic_disassembly += genotype->pic_disassembly[gene_id][gene_copy];
  rates->pic_disassembly_operations++;
}

void disassemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                           KonStates *kon_states,
//						    TimeCourse **timecoursestart, TimeCourse **timecourselast,
                           float dt, float t, float x)
{
  int gene_copy, gene_loc, gene_id;
  int j=-1;

  /* choose an appropriate gene copy to disassemble the PIC from */
  while (j < NGENES*current_ploidy && x>0) {
    j++;
    get_gene(rates->pic_disassembly_num, j, &gene_loc, &gene_copy);
    x -= genotype->pic_disassembly[state->state_change_ids[PICDISASSEMBLY_STATE][gene_copy][gene_loc]][gene_copy];
  }
  if (j==NGENES*current_ploidy) { LOG_ERROR("error in PIC disassembly\n"); }
  /* now get the gene_id */
  gene_id = state->state_change_ids[PICDISASSEMBLY_STATE][gene_copy][gene_loc];
  LOG_VERBOSE("PIC disassembly event in copy %d of gene %d\n", gene_copy, gene_id);
  /* do the disassembling */
  disassemble_PIC(state, genotype, gene_id, gene_copy, rates);
}

void transcription_init_event(GillespieRates *rates, CellState *state, Genotype *genotype,
                              KonStates *kon_states, 
//							  TimeCourse **timecoursestart, TimeCourse **timecourselast,
                              float dt, float t, float x, int *env)
{
  int gene_id;
  int gene_copy; 
  int gene_loc; 

  x /= TRANSCRIPTINIT;
LOG_ERROR("TRANSCRIPT INIT\n");
  /* choose the gene and copy that gets transcribed */
  get_gene(rates->transcript_init_num, (int)trunc(x), &gene_loc, &gene_copy);
  gene_id = state->state_change_ids[TRANSCRIPTINIT_STATE][gene_copy][gene_loc];
  LOG_VERBOSE("transcription event gene %d, copy %d\n", gene_id, gene_copy);
  LOG_ERROR("Gene id = %d\n", gene_id);

  if (state->active[gene_id][gene_copy] != ON_FULL && state->active[gene_id][gene_copy] != OFF_PIC) {
    LOG_ERROR("transcription event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt,
                                rates, kon_states, 
                                t, 
//								timecoursestart, timecourselast,
                                state->protein_conc, env);

  /* now that transcription of gene has been initiated, 
   * we add the time it will end transcription, 
   * which is dt+time of transcription from now */
  add_fixed_event_end(gene_id, gene_copy, t+dt+TTRANSCRIPTION, 
                      &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));

  /* increase the number mRNAs being transcribed */
  (state->mRNA_transcr_num[gene_id][gene_copy])++;                      
}
/*
 * END
 * Functions that handle each possible Gillespie event 
 */

/* Helper function that shifts site_ids in KonStates, tf_hindered_indexes
   and tf_bound_indexes after a position by the specified offset */
void shift_binding_site_ids(CellState *state, 
                            KonStates *kon_states,
                            int end,
                            int offset)
{
  int i, j, k, site_id;

  /* shift all sites in kon_list */
  for (i=0; i < TFGENES; i++) {
    k = 0;
    while (k < kon_states->kon_list[i]->site_count) {
      site_id = kon_states->kon_list[i]->available_sites[k];
      if (site_id >= end)
        kon_states->kon_list[i]->available_sites[k] += offset;
      k++;
    }
  }

  /* shift all sites in tf_hindered_indexes */
  j = 0;
  while (j < state->tf_hindered_num) {
    site_id = state->tf_hindered_indexes[j][0];
    if (site_id >= end)
      state->tf_hindered_indexes[j][0] += offset;
    site_id = state->tf_hindered_indexes[j][1];
    if (site_id >= end)
      state->tf_hindered_indexes[j][1] += offset;
    j++;
  }

  /* shift all sites in tf_bound_indexes */
  for (j = 0; j < state->tf_bound_num; j++) {
    site_id = state->tf_bound_indexes[j];
    if (site_id >= end)
      state->tf_bound_indexes[j] += offset;
  }
}

/*
 * recomputes all kon rates from the full cell state
 *
 * if the "recalibrate" parameter is set to true (1), this will also
 * regenerate the cached KonState data structure
 */
void recompute_kon_rates(GillespieRates *rates,
                         CellState *state,
                         Genotype *genotype,
                         KonStates *kon_states,
                         int recalibrate) 
{
  int j, k;
  float salphc;
  float orig_rates = rates->salphc;

  LOG("BEFORE rates->subtotal=%g, rates->salphc=%g\n",   rates->subtotal, rates->salphc);
  // subtract off current salphc from total
  rates->subtotal -= rates->salphc;

  // reset rates
  rates->salphc=0.0;
  rates->salphc_operations=0;

  rates->max_salphc=0.0;
  rates->max_salphc_operations=0;

  rates->min_salphc=0.0;
  rates->min_salphc_operations=0;

  // go over all the currently unbound and increment salphc
  for (k=0; k < genotype->binding_sites_num; k++) {
    int site_id;
    int tf_id = genotype->all_binding_sites[k].tf_id;
    
    int notbound = 1;
    j = 0;
    while (j < state->tf_bound_num && notbound) {
      site_id = state->tf_bound_indexes[j];
      if (site_id == k)  // site is not available, don't add to kon 
        notbound = 0;
      j++;
    }
    int nothindered = 1;
    j = 0;
    while (j < state->tf_hindered_num && nothindered) {
      site_id = state->tf_hindered_indexes[j][0];
      if (site_id == k)  // site is not available, don't add to kon 
        nothindered = 0;
      j++;
    }

    if (notbound && nothindered) {

      salphc = kon_states->konvalues[tf_id][KON_SALPHC_INDEX];

      LOG_VERBOSE("for unbound site %03d [tf_id=%02d] adding salphc=%g to rates->salphc=%g\n", k, tf_id, salphc, rates->salphc);

      rates->salphc += salphc;
      rates->salphc_operations++;

      rates->max_salphc += fmaxf(state->protein_conc[tf_id], salphc);
      rates->max_salphc_operations++;

      rates->min_salphc += fminf(state->protein_conc[tf_id], salphc);
      rates->min_salphc_operations++;

      /* only do if recalibrating state of cell from scratch */
      if (recalibrate) {
        /* update the list of sites that bind for a particular TF */
        kon_states->kon_list[tf_id]->available_sites[kon_states->kon_list[tf_id]->site_count] = k;
        (kon_states->kon_list[tf_id]->site_count)++;
        kon_states->nkon++;
      }
    } 
  }
  rates->salphc = kon*rates->salphc;
  rates->max_salphc = kon*rates->max_salphc;
  rates->min_salphc = kon*rates->min_salphc;

  /* now that it is recomputed, add salphc back to total */
  rates->subtotal += rates->salphc;
  LOG("AFTER rates->subtotal=%g, rates->salphc=%g, percent difference=%g\n",   
      rates->subtotal, rates->salphc, 100.0*(fabs(rates->salphc-orig_rates)/rates->salphc));
}

/*
 * recomputes all koff rates from the full cell state
 */
void recompute_koff_rates(GillespieRates *rates,
                          CellState *state,
                          Genotype *genotype,
                          float *koffvalues,
                          float t) 
{
  int i;
  float orig_rates = rates->koff;
  float temprate = 0.0;

  LOG("BEFORE rates->subtotal=%g, rates->koff=%g\n",   rates->subtotal, rates->koff);
  /* subtract off current koff from total */
  rates->subtotal -= rates->koff;
  
  for (i = 0; i < state->tf_bound_num; i++) {
    int site = state->tf_bound_indexes[i];
    calc_koff(site, genotype->all_binding_sites, state, &(koffvalues[i]), t);
    temprate += koffvalues[i];
  }

#if 0
  fprintf(fp_koff[state->cell_id], "%g %g\n", t, rates->koff - temprate);
#endif

  rates->koff = temprate;
  rates->koff_operations = 0;

  /* add back new koff */
  rates->subtotal += rates->koff;
  LOG("AFTER rates->subtotal=%g, rates->koff=%g, percent difference=%g\n",   
      rates->subtotal, rates->koff, 100.0*(fabs(rates->koff-orig_rates)/rates->koff));
}

/*
 * recalibrate the rates and cached data structures (KonStates) from
 * the current cell state
 */
void recalibrate_cell(GillespieRates *rates,
                      CellState *state,
                      Genotype *genotype,
                      KonStates *kon_states,
                      float **koffvalues,
                      float mRNAdecay[NGENES],
                      float transport_rate[NGENES],
                      float dt) 
{
  int i, j; 
  float protein_decay;
  float salphc = 0.0;

  /* reset the total rate for current step */
  rates->subtotal=0.0;
  
  /* reset all rates and operations */
  rates->koff=0.0;
  rates->koff_operations=0;

  rates->transport=0.0;
  rates->transport_operations=0;

  rates->mRNAdecay=0.0;
  rates->mRNAdecay_operations=0;

  rates->pic_disassembly=0.0;
  rates->pic_disassembly_operations=0;

  rates->salphc=0.0;
  rates->salphc_operations=0;

  rates->max_salphc=0.0;
  rates->max_salphc_operations=0;

  rates->min_salphc=0.0;
  rates->min_salphc_operations=0;

  /* regenerate kon_states */
  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to fall below cut-off, use aging term */
    protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;
    salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);
    kon_states->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / (protein_decay);
    kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
    kon_states->konvalues[i][KON_SALPHC_INDEX] = salphc;
    kon_states->kon_list[i]->site_count = 0;
  }
  kon_states->nkon = 0;   /* initialize to zero */

  /* regenerate kon_states and rates->{salphc,max_salphc,min_salphc} */
  recompute_kon_rates(rates, state, genotype, kon_states, 1);

  for (i=0; i < NGENES; i++) {
    /* transport rates */
    transport_rate[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
    rates->transport += transport_rate[i];
    rates->transport_operations++;

    /* regenerate decay rates */
    mRNAdecay[i] = genotype->mRNAdecay[i] * ((float) state->mRNA_cyto_num[i] + (float) state->mRNA_transl_cyto_num[i]);
    rates->mRNAdecay += mRNAdecay[i];
    rates->mRNAdecay_operations++;

  }
  /* recompute koffvalues for all sites: note that a time is required
     by the function only for diagnostic purposes, it isn't used by
     the function, and since we don't have available a current time in
     this function, we simply feed it a dummy time of 0.0  */
  recompute_koff_rates(rates, state, genotype, *koffvalues, 0.0);

  /* recompute and cache the total rate in data structure */
  rates->subtotal += rates->transport;
  rates->subtotal += rates->mRNAdecay;
  rates->subtotal += rates->pic_disassembly;

  /* 
   * convert the counts back into rates using the constants 
   */
  for (j=0; j < MAX_COPIES; j++) {
    rates->subtotal += (float) rates->acetylation_num[j] * ACETYLATE;
    rates->subtotal += (float) rates->deacetylation_num[j] * DEACETYLATE;
    rates->subtotal += (float) rates->pic_assembly_num[j] * PICASSEMBLY;
    rates->subtotal += (float) rates->transcript_init_num[j] * TRANSCRIPTINIT;    
  } 
}

/*
 * once all TF hindered indexes are computed, reallocate the
 * appropriate dynamic memory
 */
void realloc_cell_memory(CellState *state, float **koffvalues) 
{
  // TODO: check to see if we go back to exact computation of memory size needed
  //state->tf_bound_indexes = realloc(state->tf_bound_indexes, 2*(state->tf_bound_num+1)*sizeof(int));
  state->tf_bound_indexes = realloc(state->tf_bound_indexes, 10*MAXBOUND*sizeof(int));

  //state->tf_hindered_indexes = realloc(state->tf_hindered_indexes, 10*state->tf_hindered_num*sizeof(int));
  state->tf_hindered_indexes = realloc(state->tf_hindered_indexes, 100*MAXBOUND*sizeof(int));
  
  // reallocate koffvalues memory
  //*koffvalues = realloc(*koffvalues, 2*(state->tf_bound_num+1)* sizeof(float)); 
  *koffvalues = realloc(*koffvalues, 100*MAXBOUND*sizeof(float)); 

  if (!state->tf_bound_indexes || !state->tf_hindered_indexes) {
    LOG_ERROR_NOCELLID("memory allocation error cell\n");
    exit(1);
  }
}

void clone_cell(Genotype *genotype_orig,                
                Genotype *genotype_clone)
{
  int i, p;

  /* copy Genotype */

  /* clone the cis-regulatory sequence */
  memcpy(genotype_clone->cisreg_seq, genotype_orig->cisreg_seq, sizeof(char [NGENES][MAX_COPIES][CISREG_LEN]));

  /* clone all other aspects of genotype */
  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      genotype_clone->site_id_pos[i][p][0] =  genotype_orig->site_id_pos[i][p][0];
      genotype_clone->site_id_pos[i][p][1] =  genotype_orig->site_id_pos[i][p][1];
      genotype_clone->activating[i][p]=  genotype_orig->activating[i][p];
      genotype_clone->pic_disassembly[i][p]=  genotype_orig->pic_disassembly[i][p];
    }
    genotype_clone->sites_per_gene[i] = genotype_orig->sites_per_gene[i];
    genotype_clone->copies[i] = genotype_orig->copies[i];
//    genotype_clone->replication_time[i] =  genotype_orig->replication_time[i];
    genotype_clone->mRNAdecay[i] =  genotype_orig->mRNAdecay[i];
    genotype_clone->proteindecay[i] =  genotype_orig->proteindecay[i];
    genotype_clone->translation[i] =  genotype_orig->translation[i];
  }

  /* clone the TF sequence */
  memcpy(genotype_clone->tf_seq, genotype_orig->tf_seq, 
         sizeof(char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN]));

  /* clone the hindrance positions */
  for (i=0; i < TFGENES; i++) {
    genotype_clone->hindrance_positions[i] = genotype_orig->hindrance_positions[i];
  }

  /* make a pointer to original all_binding_sites, we don't modify it yet */
  genotype_clone->all_binding_sites =  genotype_orig->all_binding_sites;

  /* end Genotype copy */
}

/*
 * diagnostic function to dump a copy of the GillespieRates and some
 * of the CellState to the current error file
 */
void log_snapshot(GillespieRates *rates,
                  CellState *state,
                  Genotype *genotype,
                  KonStates *kon_states,
                  float **koffvalues,
                  float mRNAdecay[NGENES],
                  float transport[NGENES],
                  float konrate,
                  float x,
                  float t)
{
  int i, p, nkon = 0;

  LOG("snapshot at time=%g:\n x=%g, koff=%g = %d (tf_bound_num) * %g (koff/tf_bound_num)\n transport=%g\n decay=%g\n",
      t, x, rates->koff, state->tf_bound_num, rates->koff/(float)state->tf_bound_num, 
      rates->transport, rates->mRNAdecay);
  LOG_NOFUNC(" rates->salphc=%g\n rates->max_salphc=%g rates->min_salphc=%g\n", rates->salphc, rates->max_salphc, rates->min_salphc);
  LOG_NOFUNC(" konrate=%g\n", konrate);
  LOG_NOFUNC(" pic_disassembly=%g\n kon=%g = %d * %g\n",
             rates->pic_disassembly, rates->salphc+(konrate), kon_states->nkon, (rates->salphc+(konrate))/(float)kon_states->nkon);
  
  for (p=0; p < MAX_COPIES; p++) {
    LOG_NOFUNC(" acetylation=%g (copy %d)\n deacetylation=%g (copy %d)\n PIC assembly=%g (copy %d)\n transcriptinit=%g (copy %d)\n",
               (float)rates->acetylation_num[p]*ACETYLATE, p, (float)rates->deacetylation_num[p]*DEACETYLATE, p, 
               (float)rates->pic_assembly_num[p]*PICASSEMBLY, p, (float)rates->transcript_init_num[p]*TRANSCRIPTINIT, p);
  }
  LOG_NOFUNC(" total rates=%g=%g+%g\n", rates->subtotal + (konrate), rates->subtotal, konrate);
  LOG_NOFUNC(" total free=%d + total bound=%d + total hindered=%d = total sites=%d\n", 
             kon_states->nkon, state->tf_bound_num, state->tf_hindered_num, genotype->binding_sites_num);
  for (i = 0; i < TFGENES; i++) {
    nkon += kon_states->kon_list[i]->site_count;
    LOG_NOFUNC(" unoccupied binding sites=%d available for TF=%d \n", kon_states->kon_list[i]->site_count, i);
  }
  LOG_NOFUNC(" nkon recomputed=%d\n", nkon);
  LOG_NOFUNC("\n");
}

/*
 * run the model for a specified cell for a single timestep:
 *  - returns 0  if cell is not "dead" (i.e. rates haven't deteroriated to zero)
 *  - returns -1 if cell is "dead"
 */
int do_single_timestep(Genotype *genotype, 
                       CellState *state, 
                       KonStates *kon_states, 
                       GillespieRates *rates, 
                       float *t,
                       float *koffvalues,
                       float *transport_rate,
                       float mRNAdecay[NGENES],
                       float *x,
                       float *dt,
                       float *konrate,
//                       TimeCourse *timecoursestart[NPROTEINS],
//                       TimeCourse *timecourselast[NPROTEINS],
                       int maxbound2,
                       int maxbound3,
                       int no_fixed_dev_time, 
                       int burn_in,
		       int *env) 
{
  int i,j;
  int event;     /* boolean to keep track of whether FixedEvent has ended */
  int total;     /* total possible translation events */
  float fixed_time;
  float temp_t;
  
//  if (verbose) //VERBOSE JU CHANGE
//    for (j=0; j < MAX_COPIES; j++) {
//      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
////      LOG_ERROR("FDSTS rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
//    }
//    
//    for (j=0; j < MAX_COPIES; j++) {
//      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
//      //LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
////       LOG_ERROR("rates->deacetylation_num[%d]=%d\n", j, rates->deacetylation_num[j]);
//    }

  /* if enabled, check for rounding error and recompute rates */
  /* check drift from koff every RATE_OPERATIONS operations or if they go negative */
  if (recompute_koff && (rates->koff_operations >= RATE_OPERATIONS || rates->koff < 0.0)) {
    LOG_WARNING("(t=%g) after %d operations recompute koff rates\n", *t, rates->koff_operations);
    recompute_koff_rates(rates, state, genotype, koffvalues, *t);
  }
  /* check drift from various kon rates every RATE_OPERATIONS operations or if they go negative */
  if (recompute_kon && (rates->salphc_operations >= RATE_OPERATIONS || rates->salphc < 0.0)) {
    LOG_WARNING("(t=%g) after %d operations recompute kon rates\n", *t, rates->salphc_operations);
    recompute_kon_rates(rates, state, genotype, kon_states, 0);
  } 

  /* if burn-in is set and we have completed the 0.003 sec (= 0.00005
     min) burn in, then recompute the rates for the new burn-in.  The
     0.003 sec length was found by multiple runs of the model with
     kon=0.225 and finding the length of time it took for the number
     of bound TFs to equilibrate  */
  if (*t > 0.00005 && state->burn_in) {
//    printf("recalibrating cell %3d after burn-in!\n", state->cell_id);
    LOG("recalibrating cell %3d after burn-in!\n", state->cell_id);
    log_snapshot(rates, state, genotype, kon_states, &koffvalues,
                 mRNAdecay, transport_rate, *konrate, *x, *t);

    kon = kon_after_burnin;     /* reset kon to the post-burn-in value */
    recalibrate_cell(rates,
                     state,
                     genotype,
                     kon_states,
                     &koffvalues,
                     mRNAdecay,
                     transport_rate,
                     *dt); 

    log_snapshot(rates, state, genotype, kon_states, &koffvalues,
                 mRNAdecay, transport_rate, *konrate, *x, *t);
    state->burn_in = 0;    /* now disable burn-in */
  } 
  
#if 0 /* currently disable printing out the TF occupancy and amount of rounding */
  print_tf_occupancy(state, genotype->all_binding_sites, *t);
  print_rounding(state, rates, *t);
#endif
  
  *x = expdev(&seed);        /* draw random number */

  /* compute the initial dt for the next event */
  calc_dt(x, dt, rates, kon_states, mRNAdecay, genotype->mRNAdecay,
          state->mRNA_cyto_num, state->mRNA_transl_cyto_num, state->cell_id);

  if (*dt < 0.0) {
    LOG_ERROR("dt=%g is negative after calc_dt, t=%g\n", *dt, *t);
    exit(-1);
  }

	
  LOG_VERBOSE("next stochastic event due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);
  
  /* first check to see if a fixed event occurs in current t->dt window,
     or in tdevelopment if running for a fixed development time */
  fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);

  event = does_fixed_event_end(state->mRNA_transl_time_end,
                               state->mRNA_transcr_time_end,
                               state->env0_time_end,
                               state->env1_time_end,
                               fixed_time);

  /* while there are either transcription or translation events
     occuring in current t->dt window */
  while (event > 0) {
    *konrate = (*x)/(*dt);
    
    switch (event) {
    case 1:   /* if a transcription event ends */
      end_transcription(dt, *t, state, transport_rate, rates);
      LOG_ERROR("END TRANSCRIPTION\n");
      
      update_protein_conc_cell_size(state->protein_conc, state, genotype, *dt,
                                    rates, kon_states, *t,
//                                    timecoursestart, timecourselast,
                                    state->protein_conc, env);
                                    
      break;
    case 2:            /* if a translation event ends */
      *dt = state->mRNA_transl_time_end->time - *t;         /* make dt window smaller */
      total=0;  /* number of translation events */
      
      /* count current number of mRNAs that have recently arrived in cytoplasm */
      for (i=0; i<NGENES; i++) total += state->mRNA_transl_cyto_num[i];

      LOG_VERBOSE("translation event finishes out of %d possible t=%g dt=%g\n",
                  total, *t, *dt); /* bug: dt can be negative */
      
      /* get identity of gene that has just finished translating */
      i=state->mRNA_transl_time_end->gene_id;   
      
      /* there is one less mRNA that has just finished translation */
      (state->mRNA_transl_cyto_num[i])--;   
      
      /* delete the event that just happened */
      LOG_VERBOSE("delete translation event that just happened at time=%g", *t);
      delete_fixed_event_start(&(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
      
      /* there is one more mRNA that is now in cytoplasm */
      (state->mRNA_cyto_num[i])++;
      
      /* update protein concentration */
      update_protein_conc_cell_size(state->protein_conc, state, genotype, *dt,
                                    rates, kon_states, *t,
//                                    timecoursestart, timecourselast,
                                    state->protein_conc, env);
      
      /* the number of mRNAs in cytoplasm affects binding */
      change_mRNA_cytoplasm(i, genotype, state, rates, kon_states);
     
      
      break;
    case 3: // env0->env1
//    printf("time=%f, dt=%f, env= %d\n",*t, *dt, *env);
	   *dt = state->env0_time_end->time - *t;
	   *env = 1;
	   delete_fixed_event_start(&(state->env0_time_end),&(state->env0_time_end_last));	 
	   add_fixed_event(0,0,*t+*dt+duration_env1,&(state->env1_time_end),&(state->env1_time_end_last)); 
//	   printf("env1_time_end=%f\n", state->env1_time_end->time);

	   break;
	case 4:	// env1->env0
//	printf("env1_end_time=%f\n", state->env1_time_end->time);
//	 printf("time=%f, dt=%f, env= %d\n",*t, *dt, *env);
      *dt = state->env1_time_end->time - *t;
	   *env = 0;
	   delete_fixed_event_start(&(state->env1_time_end),&(state->env1_time_end_last));	 
	   add_fixed_event(0,0,*t+*dt+duration_env0,&(state->env0_time_end),&(state->env0_time_end_last)); 
	  
   	   break;
    default:
      LOG_ERROR("event=%d should never get here\n", event);
      printf("here3.5\n");
      exit(-1);
      break;
    }
    
    *t += *dt;                  /* advance time by the dt */
    *x -= (*dt)*(*konrate);
    
     
    
    LOG_VERBOSE("dt=%g t=%g fixed event old x=%g new x=%g\n", *dt, *t, (*x)+(*dt)*(*konrate), *x);
    
    /* re-compute a new dt */
    calc_dt(x, dt, rates, kon_states, mRNAdecay, 
            genotype->mRNAdecay, state->mRNA_cyto_num, state->mRNA_transl_cyto_num, state->cell_id);
    
    LOG_VERBOSE("next stochastic event (2) due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

    fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);    

    /* check to see there aren't more fixed events to do */
    event = does_fixed_event_end(state->mRNA_transl_time_end, 
                                 state->mRNA_transcr_time_end, 
//                                 state->replication_time_end,
								 state->env0_time_end,
								 state->env1_time_end,
                                 fixed_time);
  } 

  /* no remaining fixed events to do in dt, now do stochastic events */
  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */
          
  if (*t+*dt < tdevelopment || no_fixed_dev_time) {
    
    /* compute konrate */
    if (kon_states->nkon==0) {
      *konrate = (-rates->salphc);    /* all binding sites are occupied, total rate of binding is zero */
      LOG_WARNING("kon_states->nkon=%d, konrate=%g\n", kon_states->nkon, *konrate);
    } else  {
      calc_kon_rate(*dt, kon_states, konrate);  /* otherwise compute konrate */
    }

    /* if the total rates falls below zero, we do an emergency recalibration of cell */
    if (!(rates->subtotal + *konrate > 0.0)) {
        
      log_snapshot(rates, state, genotype, kon_states, &koffvalues, mRNAdecay, transport_rate, *konrate, *x, *t);
      LOG_ERROR("x should always be >0 t=%g (x=%g) rates->subtotal=%g, konrate=%g, recalibrate cell!\n", *t, *x, rates->subtotal, *konrate); 
      recalibrate_cell(rates, state, genotype, kon_states, &koffvalues, mRNAdecay, transport_rate, *dt); 
      log_snapshot(rates, state, genotype, kon_states, &koffvalues, mRNAdecay, transport_rate, *konrate, *x, *t);

      /* if this still results in either zero or negative total rates,
         this most likely due the cell being "dead" no TFs bound, no
         activity etc.  We mark cell as "dead" in this case, and
         remove from queue. */
      if (!(rates->subtotal + *konrate > 0.0)) {  
        LOG_ERROR("cell is effectively dead\n"); 
        return -1;        /* return cell status as "dead" */
      }
    }
   
    /* 
     * choose a new uniform random number weighted by the rate of all
     * Gillespie events, note that konrate is *not* included in
     * rates->subtotal, so it needs to be added here
     */
//      LOG_VERBOSE("rates subtotal = %f\n", rates->subtotal);
//      LOG_VERBOSE("koff = %f, transport = %f, mrnadecay = %f, pic dis = %f, salphc = %f, max s = %f, min s = %f\n", 
//    rates->koff, rates->transport, rates->mRNAdecay, rates->pic_disassembly, rates->salphc, rates->max_salphc, rates->min_salphc);
    
//    for(i=0; i<MAX_COPIES; i++){
//             LOG_VERBOSE("ace_num[%d] = %d, deace_num[%d] = %d, pic[%d] = %d, trans_num[%d] = %d, pic_dis[%d] = %d\n",
//             i,rates->acetylation_num[i],i,rates->deacetylation_num[i],i, rates->pic_assembly_num[i],i, rates->transcript_init_num[i],i, rates->pic_disassembly_num[i]);
//    }
    *x = ran1(&seed)*(rates->subtotal + *konrate);  
    //LOG_ERROR("*x = %f, rates + kon  = %f \n", *x, (rates->subtotal+ *konrate ));
   // LOG_ERROR("rates koff = %f\n", rates->koff);
   // LOG_ERROR("rates pic = %f\n sum rate ACE = %f\n sum rate DEACE = %f\n sum rate PIC = %f\n sum rate trans = %f\n", rates->pic_disassembly, sum_rate_counts(rates->acetylation_num) * ACETYLATE, sum_rate_counts(rates->deacetylation_num) * DEACETYLATE, sum_rate_counts(rates->pic_assembly_num) * PICASSEMBLY, sum_rate_counts(rates->transcript_init_num) * TRANSCRIPTINIT);
    //printf("*x = %f\n", *x);
    //system("PAUSE");
    if (verbose) {
      log_snapshot(rates,
                   state,
                   genotype,
                   kon_states,
                   &koffvalues,
                   mRNAdecay,
                   transport_rate, 
                   *konrate,
                   *x,
                   *t);
    }
    /* 
     * STOCHASTIC EVENT: a TF unbinds (koff) 
     */
     //printf("*x = %f, rates transport= %f, rates koff= %f\n", *x, rates->transport, rates->koff);
     LOG_ERROR("*x = %f\n", *x);
     LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );

    if (*x < rates->koff) {     // printf("here4\n");
           LOG_ERROR("unbinding event\n");
      int ignore_event;
      tf_unbinding_event(rates, state, genotype, kon_states, koffvalues,
//                         timecoursestart, timecourselast, 
						 (*konrate), *dt, *t, *x, 1, &ignore_event, env);
      if (ignore_event) { /* if the TF unbinding event failed, we return without updating *t  */
     // printf("here5.1\n");
        return 0;
      }
    } else {
           LOG_ERROR("x = %f\n", *x);
      *x -= rates->koff;  
      LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
      //printf("*x = %f\n", *x);
      /* 
       * STOCHASTIC EVENT: a transport event
       */
//       printf("here5.2\n");
      if (*x < rates->transport) {    // printf("here5.3\n"); printf("%f\n",rates->transport);
             LOG_ERROR("transport event\n");
        transport_event(rates, state, genotype, kon_states, transport_rate, 
//                        timecoursestart, timecourselast,
						 *dt, *t, *x, env);
                        
               // printf("here5.3.1\n");      
       for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("2rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);}
    
      } else {
        LOG_ERROR("x = %f\n", *x); 
        *x -= rates->transport;
        LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
        /* 
         * STOCHASTIC EVENT: an mRNA decay event
         */
        if (*x < rates->mRNAdecay) {  //printf("here5.4\n");
               LOG_ERROR("decay event\n");
          mRNA_decay_event(rates, state, genotype, kon_states, mRNAdecay,
//                           timecoursestart, timecourselast,
						    *dt, *t, *x, env);
        } else {
               LOG_ERROR("x = %f\n", *x); 
          *x -= rates->mRNAdecay;
          LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
          /* 
           * STOCHASTIC EVENT: PIC disassembly
           */
          if (*x < rates->pic_disassembly) {//printf("here5.5\n");
                 LOG_ERROR("pic disassembly event\n");
            disassemble_PIC_event(rates, state, genotype, kon_states, 
//                                  timecoursestart,  timecourselast, 
								  *dt, *t, *x);
          } else {
                 LOG_ERROR("x = %f\n", *x);
            *x -= rates->pic_disassembly;
            LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
            /* 
             * STOCHASTIC EVENT: TF binding event
             */
            if (*x < rates->salphc + (*konrate)) {   /* get total on rate = sum of (salphc) and (konrate) *///printf("here5.6\n");
              LOG_ERROR("binding event\n");
              tf_binding_event(rates, state, genotype, kon_states, koffvalues,
//                               timecoursestart, timecourselast, 
							   (*konrate), *dt, *t, 
                               maxbound2, maxbound3, 1, env);  
							  
            } else {
                   LOG_ERROR("Before x = %f\n", *x); 
              *x -= (rates->salphc + (*konrate));
              LOG_ERROR("after x = %f\n", *x);
              LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
              /* 
               * STOCHASTIC EVENT: histone acetylation
               */
                LOG_ERROR("ACHECK! x = %f,  sum_rate_counts = %d, acetylate = %f, PIC num = %f\n", *x, sum_rate_counts(rates->acetylation_num), ACETYLATE, sum_rate_counts(rates->pic_assembly_num) * PICASSEMBLY );
              if (*x < (float) sum_rate_counts(rates->acetylation_num) * ACETYLATE) {//printf("here5.7\n");
                LOG_ERROR("hist act event\n");
                histone_acteylation_event(rates, state, genotype, kon_states, 
//                                          timecoursestart, timecourselast,
										   *dt, *t, env);
              } else { 
                
                *x -= (float) sum_rate_counts(rates->acetylation_num) * ACETYLATE;
                /* 
                 * STOCHASTIC EVENT: histone deacetylation
                 */
                 LOG_ERROR("CHECK! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
                if (*x < (float) sum_rate_counts(rates->deacetylation_num) * DEACETYLATE) {
                  LOG_ERROR("deact event\n"); 
                  histone_deacteylation_event(rates, state, genotype, kon_states, 
//                                              timecoursestart, timecourselast, 
											  *dt, *t, env);
                } else {
                  *x -= (float) sum_rate_counts(rates->deacetylation_num) * DEACETYLATE; 
                  /* 
                   * STOCHASTIC EVENT: PIC assembly
                   */
                  if (*x < (float) sum_rate_counts(rates->pic_assembly_num) * PICASSEMBLY) {
                         LOG_ERROR("pic assembly event\n");
                    assemble_PIC_event(rates, state, genotype, kon_states, 
//                                       timecoursestart, timecourselast, 
									   *dt, *t, env);
                  } else {
                    *x -= (float) sum_rate_counts(rates->pic_assembly_num) * PICASSEMBLY;
                    /* 
                     * STOCHASTIC EVENT: transcription initiation
                     */
                    if (*x < (float) sum_rate_counts(rates->transcript_init_num) * TRANSCRIPTINIT) {
                           LOG_ERROR("transcript init event time = %f\n", *t);
                      transcription_init_event(rates, state, genotype, kon_states, 
//                                               timecoursestart, timecourselast, 
											   *dt, *t, *x, env);
                    } else {
                      /*
                       * FALLBACK: shouldn't get here, previous
                       * events should be exhaustive
                       */
                      
                      LOG_ERROR("[cell %03d] t=%g no event assigned: x=%g, rates->subtotal+konrate=%g, recalibrate cell\n", 
                                state->cell_id, *t, *x, rates->subtotal + *konrate);

                      log_snapshot(rates, state, genotype, kon_states, &koffvalues, mRNAdecay,
                                   transport_rate, *konrate, *x, *t);
                      recalibrate_cell(rates, state, genotype, kon_states, &koffvalues,
                                       mRNAdecay, transport_rate, *dt); 
                      log_snapshot(rates, state, genotype, kon_states, &koffvalues, mRNAdecay,
                                   transport_rate, *konrate, *x, *t);
                      if (*t > 1000.0) {
                        LOG_ERROR("cell has been running for %g which is > 1000 mins\n", *t);
                      }
                      /* ignore this timestep, we throw out the dt by returning without updating t */
                      return 0;
                    }
                  }
                }
              }
            }
          }       
        }
      }
    }
    
    
    
    /* Gillespie step: advance time to next event at dt */
    *t += *dt;
    LOG_VERBOSE("dt=%g t=%g\n", *dt, *t);
  } else {
     
         
    /* we will reach the end of development in dt */
    LOG_VERBOSE("finish at t=%g dt=%g\n", *t, *dt);
    
    /* do remaining dt */
    *dt = tdevelopment - *t;
    
    /* final update of protein concentration */
    update_protein_conc_cell_size(state->protein_conc, state, genotype, *dt,
                                  rates, kon_states,*t,
//								   timecoursestart, timecourselast,
                                  state->protein_conc, env);
    /* advance to end of development (this exits the outer while loop) */
    *t = tdevelopment;
    
    
  }
  
   
  return 0;
}


void free_fixedevent(CellState *state)
{
	FixedEvent *temp1, *temp2;
	
	temp1=state->env1_time_end;
	
	while(temp1){		
		temp1=temp1->next;
		temp2=temp1;
		free(temp2);	
	}
	
	temp1=state->env0_time_end;
	while(temp1){		
		temp1=temp1->next;
		temp2=temp1;
		free(temp2);	
	}
	
	temp1=state->mRNA_transcr_time_end;
	while(temp1){		
		temp1=temp1->next;
		temp2=temp1;
		free(temp2);	
	}
	
	temp1=state->mRNA_transl_time_end;
	while(temp1){		
		temp1=temp1->next;
		temp2=temp1;
		free(temp2);	
	}
}

float calc_avg_growth_rate(int current_genotype,
                           Genotype *genotype, 
                           CellState *state, 
                           float init_mRNA[NGENES],
                           float init_protein_conc[NGENES],
                           float *t,
                           float *koffvalues,
                           float transport_rate[NGENES],
                           float mRNAdecay[NGENES],
                           float *x,
                           float *dt,
                           float *konrate,
                           KonStates *kon_states,
                           GillespieRates *rates,
                           float maxbound2,
                           float maxbound3,
                           float temperature,
                           int no_fixed_dev_time)
{
	
	int cell_status=0; 
	int i,j;
	int env=0;
	float avg_GR=0;
	
	state->tf_bound_indexes=NULL; //moved from init_cell
	state->tf_hindered_indexes = NULL;

	initialize_cell_cache(state, *genotype, kon_states, &koffvalues, maxbound2, maxbound3);
	
	for(i=0;i<N_replicates;i++)
	{	 
            env=0;
		
            printf("%d\n",i);
			
		
            for (j=0; j < NGENES; j++) {
		init_protein_conc[j] = exp(1.25759*gasdev(&seed)+7.25669);
		init_mRNA[j] = exp(0.91966*gasdev(&seed)-0.465902);
            }
		
            initialize_cell(state, current_genotype, genotype->copies, genotype->mRNAdecay, init_mRNA, init_protein_conc, burn_in);
	
            state->temperature = temperature;
	    state->RTlnKr = GASCONSTANT * temperature * log(KR);
	    	    
	    calc_from_state(genotype, state, rates, kon_states, transport_rate, mRNAdecay);
	    	
	    *t = 0.0; 
                     
            
            
            do_single_timestep(genotype, 
	                       state, 
	                       kon_states, 
	                       rates, 
	                       t,
	                       koffvalues,
	                       transport_rate,
	                       mRNAdecay,
	                       x,
	                       dt,
	                       konrate,
	                       maxbound2,
	                       maxbound3,
	                       no_fixed_dev_time,
	                       burn_in,
                                &env);
						   
		
				   
		while(*t<tdevelopment){  		                                 
	    	   
		    cell_status = do_single_timestep(genotype, 
		                                     state, 
		                                     kon_states, 
		                                     rates, 
		                                     t,
		                                     koffvalues,
		                                     transport_rate,
		                                     mRNAdecay,
		                                     x,
		                                     dt,
		                                     konrate,
		                                     maxbound2,
		                                     maxbound3,
		                                     no_fixed_dev_time,
		                                     burn_in,
						     &env);
											 
			if(cell_status==-1) {
				printf("tdevelopment=%f, t=%f\n",tdevelopment, *t);
				printf("dead!\n");
				break;				
			}
		 		    
	   }
	   
	   printf("t=%f\n",*t);
	   
	   avg_GR += log(state->cell_size)/tdevelopment; 
	   	   
	   free_fixedevent(state);
	   			   
	}
	
	free(state->tf_bound_indexes);
    free(koffvalues);
    free(state->tf_hindered_indexes);
	
	state->tf_bound_indexes=NULL;
	koffvalues=NULL;
	state->tf_hindered_indexes=NULL; 
	
	for(i=0;i<NGENES;i++){
		
		free(kon_states->kon_list[i]->available_sites);
		kon_states->kon_list[i]->available_sites= NULL;
		free(kon_states->kon_list[i]);
		kon_states->kon_list[i]=NULL;
		
	}
	
	return avg_GR/N_replicates;

}
/*
 * initialize the population of cells, then run for a given
 * number of divisions or run cell(s) for a specific length of time
 */
void init_run_pop(Genotype genotype[N_para_threads],
                  CellState state[N_para_threads],
//                  TimeCourse *timecoursestart[2][NPROTEINS],
//                  TimeCourse *timecourselast[2][NPROTEINS],
                  float temperature,   /* in Kelvin */
                  float kdis[NUM_K_DISASSEMBLY],
                  int output_binding_sites,
                  int no_fixed_dev_time)
//                  int max_divisions)
{
  
  int i, j;
  int current_genotype = 0;
  float current_fitness;
  int maxbound2, maxbound3; 
  float init_mRNA[N_para_threads][NGENES], init_protein_conc[N_para_threads][NGENES];
  float t[N_para_threads];             /* time of last event */
  float *koffvalues[N_para_threads];   /* rates of unbinding */
  float transport_rate[N_para_threads][NGENES];  /* transport rates of each mRNA */
  float mRNAdecay[N_para_threads][NGENES];  /* mRNA decay rates */
  float x[N_para_threads];                  /* random number */
  float dt[N_para_threads];                 /* delta-t */
  float konrate[N_para_threads];
	
  KonStates kon_states[N_para_threads];
  GillespieRates rates[N_para_threads];
	
  maxbound2 = MAXBOUND;
  maxbound3 = 10*MAXBOUND;
 
  output=1;
    /* for all cells in this replicate initialize all parts of the
       genotype, *except* for the cis-reg sequence using the initial
       genes[0]  */
  initialize_genotype(&genotype[current_genotype], &genotype[current_genotype], kdis, current_genotype); // checked
  
  for(i=1;i<N_para_threads;i++) // the current genotype is always the first element in genotype[N_para_threads]
  {
  	clone_cell(&genotype[current_genotype],&genotype[i]);  	
  }
  
  current_fitness = calc_avg_growth_rate(current_genotype,
                                         &genotype[current_genotype],
                                         &state[current_genotype],
                                         &init_mRNA[current_genotype][0],
                                         &init_protein_conc[current_genotype][0],
                                         &t[current_genotype],
                                         koffvalues[current_genotype],
                                         &transport_rate[current_genotype][0],
                                         &mRNAdecay[current_genotype][0],
                                         &x[current_genotype],
                                         &dt[current_genotype],
                                         &konrate[current_genotype],
                                         &kon_states[current_genotype],
                                         &rates[current_genotype],
                                         maxbound2,
                                         maxbound3,
                                         temperature,
                                         no_fixed_dev_time);										 
	
	 									 
	  printf("current_fitness=%f\n", current_fitness);
//  omp_set_num_threads(N_para_threads);	
//  #pragma omp parallel
//  {	
//  	  int ID = omp_get_thread_num();
//	  current_fitness = calc_avg_growth_rate(ID,
//	  									     &genotype[ID],
//	  										 &state[ID],
//											 &init_mRNA[ID][NGENES],
//											 &init_protein_conc[ID][NGENES],
//											 &t[ID],
//											 koffvalues[ID],
//											 &transport[ID][NGENES],
//											 &mRNAdecay[ID][NGENES],
//											 &x[ID],
//											 &dt[ID],
//											 &konrate[ID],
//											 &kon_states[ID],
//											 &rates[ID],
//											 maxbound2,
//											 maxbound3,
//											 temperature,
//											 no_fixed_dev_time);
//
////	    free(koffvalues[ID]);
////	    
//	    for (i=0; i < NPROTEINS; i++) {
//	      free(kon_states[ID].kon_list[i]->available_sites);
//	      free(kon_states[ID].kon_list[i]);
//	    }
//	  									 
//	  printf("ID %d, current_fitness=%f\n", ID, current_fitness);
//  }
}

//void print_time_course(TimeCourse *start,
//                       int i,
//                       int j)
//{
//  FILE *fpout;
//  char filename[80];
//  
//  /* do the normal thing on the first cell */
//  if (POP_SIZE == 1)
//    sprintf(filename, "%s/p%d.dat", output_directory, i);
//  else
//    sprintf(filename, "%s/p%d-%d.dat", output_directory, j, i);
//  if ((fpout = fopen(filename,"w"))==NULL) {
//    LOG_ERROR_NOCELLID("error: Can't open %s file\n", filename);
//  }
//  while (start) {
//    fprintf(fpout,"%g %g\n", start->time, start->concentration);
//    start = start->next;
//  }
//  fclose(fpout);  
//}
////
//void print_all_protein_time_courses(TimeCourse *timecoursestart[2][NPROTEINS],
//                                    TimeCourse *timecourselast[2][NPROTEINS])
//{
////  int i, j;
////  for (j = 0; j < 1; j++) {
////    for (i=0; i < NPROTEINS; i++) {
////      if ((output)) print_time_course(timecoursestart[j][i], i, j);
////    }
////  }
//print_time_course(timecoursestart[0][SELECTION_GENE_A],1,0);
//print_time_course(timecoursestart[0][SELECTION_GENE_B],2,0);
//}

