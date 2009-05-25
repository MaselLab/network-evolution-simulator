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
#include <sys/stat.h>
#include <errno.h>

/* local includes */
#include "random.h"
#include "lib.h"
#include "priority-queue.h"
#include "netsim.h"

const int maxelements=500*MAX_COPIES; 
/* start by allocating maxelements when initializing a genotype, double as needed, reduce at end */
const int maxbound=500*MAX_COPIES;
const int nmin=4;
const float kRNA=618.0;
const float ttranslation=1.0;
const float ttranscription=1.0;
const float pact=0.62;
const float transcriptinit=8.5; /* replace betaon and betaoff */
const float deacetylate=0.462;
const float acetylate=0.1155;
const float PICassembly=0.0277;
const float startnucleus=0.1;
const float Kr=10.0;    /* don't put this less than 1, weird things happen to koff calculation */
const float GasConstant=8.31447;
const float cooperativity=1.0;/* dGibbs, relative to 1 additional specific nt */
const float cooperative_distance=11;  /* distance co-operativity operates, changed from 20 */ 
const float NumSitesInGenome = 1.3e+6; /* updated from 1.8e+6 */
const float selection = 1.0;

const float mN = 0.1;  // TODO: mN: integrate with mutation rate
const int Generations=5;

float kon=1e-4; /* lower value is so things run faster */
float kon_after_burnin=1e-4; /* lower value is so things run faster */
/* kon=0.2225 is based on 1 molecule taking 240seconds=4 minutes
   and 89% of the proteins being in the nucleus*/
int burn_in = 0;

float tdevelopment = 120.0;/* default  development time: can be changed at runtime */
float timemax = -1.0;      /* set an upper limit to development time (default to -1.0=no limit) */
int current_ploidy = 2;    /* ploidy can be changed at run-time: 1 = haploid, 2 = diploid */
int output = 0;
long seed = 28121;         /* something is wrong here: changing seed changes nothing */
int dummyrun = 4;          /* used to change seed */
int recompute_koff = 0;    /* toggle whether to recompute certain features at each time to avoid
                              compounding rounding error (off by default) */
int recompute_kon = 0;
float critical_size = 1.0; /* critical size at which cell divides, 
                              set to negative to prevent division  */
float growth_rate_scaling = 2.0; /* set growth rate scaling factor */
float time_s_phase = 30.0;  /* length of S phase (default: 30 mins) */
float time_g2_phase = 30.0; /* length of G2/mitosis phase (default: 30 mins) */
int random_replication_time = 0; /* toggle replication times in S phase being random */

/* file output parameters */
char *output_directory = "output";
int verbose = 0;

/* protein aging term: used when c=c'+g=0, set to 1e-4 < mean-3*sd of
   Belle et al. (2006) and small with respect to most growth rates  */
float protein_aging = 1e-4;

/* initialize the growth rate parameters: 
 * do computations here so that we can easily change the scaling factor and Pp */
void initialize_growth_rate_parameters() {
  float hc, gpeak, Ltf;
  gpeak = 0.005776*growth_rate_scaling;  /* in min^-1 based on doubling time of 120 min: ln(2)/(120 min)=0.005776 */
  Pp = 12000;              /* mean gene expression of all proteins is 12064.28 */
  Ltf= 1418;               /* mean gene expression of only TFs is 1418 */
  hc = (gpeak/Pp)*(1-(log(2-2*0.2)/log(2)));      /* in min^-1 cost of doubling gene expression, based on Wagner (2005) 
                                                   * using {s=0.2, N=500} matches {s=10^-5, N=10^7} combination (both Ns=100) */
  h = hc/0.023;            /* using c=0.023/min from mean of distribution from Belle et al (2006)*/
  gmax = gpeak + hc*(Pp+(TFGENES*Ltf));    /* compute the gmax coefficient based on gpeak and other parameters */
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

  //printf("len=%d, NGENES=%d, ploidy=%d, current_element=%d\n", len, NGENES, ploidy, current_element); 
  for (i=0; i<len/ploidy; i++) {
    first = (i / current_element)*ploidy*current_element + i % current_element;
    second = first + current_element;
    third = second + current_element;
    fourth = third + current_element;
    //printf("first=%d, second=%d, third=%d, fourth=%d\n", first, second, third, fourth); 
    x = ran1(&seed);
    
    Seq[first] = set_base_pair(x);
    /* clone the randomly chosen sequence for all other sites */
    Seq[second] = Seq[first];
    Seq[third] = Seq[first];
    Seq[fourth] = Seq[first];
  }
  //printf("length: %d, sequence is %s\n", strlen(Seq), Seq);
}

void print_genotype(Genotype *genotype, int genotypeID) {
  int i, p;

  printf("[genotype %03d] hind_pos: ", genotypeID);
  for (i=0; i < TFGENES; i++) {
    printf("%d ", genotype->hindrance_positions[i]);
  }
  printf("\n");

  for (i=0; i < NGENES; i++) {
    printf("[genotype %03d gene %02d] ", genotypeID, i);
    printf("repl=%g, ", genotype->replication_time[i]);
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

  // TODO: currently comment-out generation of TF bound information
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

void initialize_genotype_fixed(Genotype *indiv, 
                               float kdis[],
                               int genotypeID)
{
  int i, j, p;

  LOG_NOCELLID("[genotype %03d] activators vs repressors ", genotypeID);

  /* initialize hindrance for all TFGENES */
  for (p=0; p < TFGENES; p++) {
    if (HIND_LENGTH == TF_ELEMENT_LEN) {
      indiv->hindrance_positions[p]=0;
    } else  {
      // TODO: only keep until new regression output is generated
#ifdef USE_RAND
      indiv->hindrance_positions[p]=rand()%10;
#else
      indiv->hindrance_positions[p]=rint(ran1(&seed)*(HIND_LENGTH - TF_ELEMENT_LEN));
#endif

    }
  } 
  
  for (i=0; i < NGENES; i++) {
    indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    while (indiv->mRNAdecay[i]<0.0)
      indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    indiv->proteindecay[i]=-1.0;
    while (indiv->proteindecay[i] < 0.0) {
      if (ran1(&seed) < 0.08421)
        indiv->proteindecay[i] = 0.0;
      else indiv->proteindecay[i] = exp(0.7874*gasdev(&seed)-3.7665);
    }
    /* dilution no longer done here, because it is now variable (function of instantaneous growth rate) */
    indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    while (indiv->translation[i] < 0.0)
      indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);

    /* make the activations the same in each copy */
    if (ran1(&seed)<pact) {
      for (p=0; p < MAX_COPIES; p++) 
        indiv->activating[i][p] = 1;
    } else {
      for (p=0; p < MAX_COPIES; p++) 
        indiv->activating[i][p] = 0;
    }

    for (p=0; p < MAX_COPIES; p++) 
      LOG_NOFUNC("%d ", indiv->activating[i][p]);

    j = trunc(NUM_K_DISASSEMBLY * ran1(&seed));
    
    for (p=0; p < MAX_COPIES; p++) 
      indiv->pic_disassembly[i][p] = kdis[j];
  }
  LOG_NOFUNC("\n");

  /* for each gene determine a time point during the [0, time_s_phase min] S-phase interval */
  for (i=0; i < NGENES; i++) {
    if (random_replication_time) 
      indiv->replication_time[i] = time_s_phase*ran1(&seed); 
    else
      indiv->replication_time[i] = time_s_phase*(i/(float)NGENES);
    LOG_VERBOSE_NOCELLID("[genotype %03d] offset for replication time after S-phase starts: %g\n", genotypeID, indiv->replication_time[i]);
  }
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, but the same random TF sequence,
 * hindrance positions, replication times, etc.  (full list below)
 */
void initialize_genotype(Genotype *indiv, 
                         Genotype *clone,
                         float kdis[],
                         int genotypeID)
{
  int p;
  
  initialize_sequence((char *)indiv->cisreg_seq, CISREG_LEN*MAX_COPIES*NGENES, MAX_COPIES, NGENES);

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
  if (genotypeID == 0) {
    initialize_sequence((char *)indiv->tf_seq, TF_ELEMENT_LEN*MAX_COPIES*TFGENES, MAX_COPIES, TFGENES);
    initialize_genotype_fixed(indiv, kdis, genotypeID);
  // otherwise clone it from the first instance
  } else {
    /* copy the TF sequence data */
    memcpy(indiv->tf_seq, clone->tf_seq, sizeof(char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN]));
    initialize_new_cell_genotype(indiv, clone);
  }

  /* start number of copies of gene at current_ploidy */
  for (p=0; p < NGENES; p++) {
    indiv->copies[p] = current_ploidy;
  }

  calc_all_binding_sites(indiv->copies, indiv->cisreg_seq, indiv->tf_seq, 
                         &(indiv->binding_sites_num), &(indiv->all_binding_sites), indiv->hindrance_positions,
                         indiv->sites_per_gene, indiv->site_id_pos);
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
                                int *maxAlloc,
                                int gene_id,
                                int gene_copy,
                                int hind_pos[TFGENES])
{
  int i, j, tfind, match, maxBindingSiteAlloc;

  maxBindingSiteAlloc = *maxAlloc;

  for (i=0; i < CISREG_LEN-TF_ELEMENT_LEN; i++) {  /* scan forwards */
    for (tfind=0; tfind < TFGENES; tfind++) {      /* only loop through TF genes */
      match=0;
      for (j=i; j < i+TF_ELEMENT_LEN; j++) {
        if (cisreg_seq[gene_id][gene_copy][j] == tf_seq[tfind][gene_copy][j-i])
          match++;
      }
      if (match >= nmin){
        if (binding_sites_num + 1 >= maxBindingSiteAlloc) {
          maxBindingSiteAlloc = 2*maxBindingSiteAlloc;
          *all_binding_sites = realloc(*all_binding_sites, maxBindingSiteAlloc*sizeof(AllTFBindingSites));
          if (!(*all_binding_sites)) {
            LOG_ERROR_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d failed.\n", maxBindingSiteAlloc);
            exit(1);
          }
          else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d succeeded\n", maxBindingSiteAlloc);
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
      if (match >= nmin){
        if (binding_sites_num + 1 >= maxBindingSiteAlloc){
          maxBindingSiteAlloc = 2*maxBindingSiteAlloc;
          *all_binding_sites = realloc(*all_binding_sites, maxBindingSiteAlloc*sizeof(AllTFBindingSites));
          if (!(*all_binding_sites)){
            LOG_ERROR_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d failed.\n", maxBindingSiteAlloc);
            exit(1);
          }
          else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d succeeded\n", maxBindingSiteAlloc);
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
  *maxAlloc = maxBindingSiteAlloc;
  return binding_sites_num;
}

/*
 * compute the list of binding sites for the specified number of gene
 * copies
 */
void calc_all_binding_sites(int copies[NGENES],
                            char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN],
                            char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                            int *newBindSiteCount,
                            AllTFBindingSites **all_binding_sites,
                            int hind_pos[TFGENES],
                            int sites_per_gene[NGENES],
                            int site_id_pos[NGENES][MAX_COPIES][2])
{
  int p, maxBindingSiteAlloc, binding_sites_num;
  int gene_id;

  maxBindingSiteAlloc = maxelements;
  *all_binding_sites = malloc(maxBindingSiteAlloc*sizeof(AllTFBindingSites));
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
                                                    &maxBindingSiteAlloc,
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
  *newBindSiteCount = binding_sites_num;
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

void add_fixed_event_end(int i,
                         int p,
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
  newtime->gene_id = i;
  newtime->copy = p;
  newtime->time = t;
  LOG_VERBOSE_NOCELLID("adding event end at time=%f for gene=%d (copy=%d)\n", t, i, p);
  sls_store_end(newtime, start, last);
}

void delete_fixed_event(int gene_id,
                        int p,
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
    if ((info->gene_id==gene_id && info->copy==p)) {
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
                       j+1, i, gene_id, p);
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
                     int cell_id,
                     int copies[NGENES],
                     float mRNAdecay[NGENES],
                     float meanmRNA[NGENES],
                     float initProteinConc[NPROTEINS],
                     int burn_in)
{
  int i, j, k, totalmRNA;
  float t;

  /* initialize the ID of the cell */
  state->cell_id = cell_id;

  /* initialize the founder cell using the current ID */
  state->founder_id = cell_id;

  /* don't start in S phase */
  state->in_s_phase = 0;

  /* initialize whether to do kon burn-in or not */
  state->burn_in = burn_in;

  /* start cell size at 0.5 */
  state->cell_size = 0.5;

  /* start cell with no divisions */
  state->divisions = 0;

  /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
  state->growth_rate = 0.0;

  state->mRNA_transcr_time_end = NULL;
  state->mRNA_transcr_time_end_last = NULL;
  state->mRNA_transl_time_end = NULL;
  state->mRNA_transl_time_end_last = NULL;
  state->replication_time_end = NULL;
  state->replication_time_end_last = NULL;

  state->tf_bound_num = 0;  /* initialize with nothing bound */
  state->tf_hindered_num = 0;
  state->tf_bound_indexes = NULL;
  state->tf_hindered_indexes = NULL;

  for (i=0; i < NGENES; i++) {

    for (j=0; j < MAX_COPIES; j++) {
      state->active[i][j] = ON_WITH_NUCLEOSOME;
    }

    totalmRNA = (int) poidev(meanmRNA[i],&seed);
    state->mRNA_nuclear_num[i] = (int) bnldev(startnucleus, totalmRNA, &seed);
    state->mRNA_cyto_num[i] = totalmRNA - state->mRNA_nuclear_num[i];
    state->mRNA_transl_cyto_num[i] = 0;

    for (k=0; k<state->mRNA_cyto_num[i]; k++) {
      t = expdev(&seed) / mRNAdecay[i];
      if (t < ttranslation) {
        (state->mRNA_cyto_num[i])--;
        (state->mRNA_transl_cyto_num[i])++;
        LOG_VERBOSE("add translation event time=%g for gene=%d\n", (ttranslation-t), i);
        add_fixed_event(i, -1, ttranslation-t, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
      }
    } 

    int total_mRNA_transcribing = (int) poidev(meanmRNA[i]*ttranscription*mRNAdecay[i], &seed);
    
    /* split it up evenly between the copies */
    int mRNA_copy1 = trunc(total_mRNA_transcribing/current_ploidy);
    int mRNA_copy2 = total_mRNA_transcribing - mRNA_copy1;

    for (j=0; j < MAX_COPIES; j++) {
      if (j < current_ploidy)  {
        state->mRNA_transcr_num[i][j] = (j==0) ? mRNA_copy1 : mRNA_copy2;
        LOG_VERBOSE_NOCELLID("initializing state->mRNA_transcr_num[%2d][%2d]=%d\n", i, j, state->mRNA_transcr_num[i][j]);
        for (k=0; k < state->mRNA_transcr_num[i][j]; k++)
          add_fixed_event(i, j, ran1(&seed)*ttranscription, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));
      } else {
        state->mRNA_transcr_num[i][j] = 0;
      }
    }
  }
  for (i=0; i < NPROTEINS; i++) {
    state->protein_conc[i] = initProteinConc[i];
  }
}

/*
 * initialize the cell state that are "cached", i.e. that are
 * maintained for code efficiency and are not part of the underlying
 * molecular biology.  All of these data structures can be regenerated
 * from the cell state at any point in the simulation
 */
void initialize_cell_cache(CellState *state,
                           Genotype genes,
                           KonStates *konStates,
                           float **koffvalues,
                           int maxbound2,
                           int maxbound3)
{
  int i;
  /* number of possible binding sites */
  // TODO: currently create for all proteins, not just TFs, as we need info for protein decay
  for (i=0; i < NPROTEINS; i++){
    konStates->kon_list[i] = malloc(sizeof(KonList));
    konStates->kon_list[i]->available_sites = malloc(genes.binding_sites_num*sizeof(int));
  }
  
  state->tf_bound_indexes = realloc(state->tf_bound_indexes, maxbound2*sizeof(int));
  *koffvalues = malloc(maxbound2*sizeof(float)); 
  state->tf_hindered_indexes = realloc(state->tf_hindered_indexes, 2*maxbound3*sizeof(int));
  
  if (!konStates->konvalues || !state->tf_bound_indexes || !koffvalues ||
      !state->tf_hindered_indexes || !konStates->kon_list) {
    LOG_ERROR("memory allocation error at start of develop\n");
    exit(1);
  }
}

/* could perhaps be a little faster with option to skip *df calculation for first 2 calls */
void calc_time (float t, 
                float x, 
                GillespieRates *rates,
                KonStates *konStates,
                float *f, 
                float *df)
{
  float r, denom, numer, ct, ect;
  int i;
  
  r = numer = 0.0;

  /* loop over all proteins (TFs and non-TFs) */
  for (i=0; i < NPROTEINS; i++) {

    LOG_VERBOSE_NOCELLID("t=%g, konStates->kon_list[%d]->site_count=%d, konStates->konvalues[%d][KON_DIFF_INDEX]=%g\n", 
                         t, i, konStates->kon_list[i]->site_count, i, konStates->konvalues[i][KON_DIFF_INDEX]);

    /* if this particular TF is bound somewhere */
    if (konStates->kon_list[i]->site_count > 0) {
      ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] * t;
      if (fabs(ct)<EPSILON) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) konStates->kon_list[i]->site_count) * konStates->konvalues[i][KON_DIFF_INDEX] * ect;
      numer += ((float) konStates->kon_list[i]->site_count) * konStates->konvalues[i][KON_DIFF_INDEX] * (ect-ct*exp(-ct));
    }
  }
  numer *= kon;
  r *= kon;
  denom = r + t*(rates->total + rates->salphc);
  denom = denom * denom;
  r /= t;

  /* compute delta-t */
  *f = x/(r + rates->total + rates->salphc) - t;

  LOG_VERBOSE_NOCELLID("x=%g, r=%g, rates->total=%g, rates->salphc=%g, f=%g\n", x, r, rates->total, rates->salphc, *f);

  /* compute derivative of equation */
  *df = x*numer/denom - 1.0;

  LOG_VERBOSE_NOCELLID("t=%g f=%g df=%g\n", t, *f, *df);
}

void calc_kon_rate(float t,
                   KonStates *konStates,
                   float *konrate)
{
  float r,ct,ect;
  int i;
  
  r = 0.0;
  /* loop through all TFs */
  for (i=0; i < TFGENES; i++) {  
    if (konStates->kon_list[i]->site_count > 0) {
      ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]*t;
      if (fabs(ct)<EPSILON) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) konStates->kon_list[i]->site_count)*konStates->konvalues[i][KON_DIFF_INDEX]*ect;
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
                           Genotype *genes,
                           CellState *state,
                           GillespieRates *rates,
                           KonStates *konStates)
{
  float salphc; 
  
  /* number of mRNAs in cytoplasm affects kon rates */
  salphc = (float) (state->mRNA_cyto_num[i]) * genes->translation[i] / konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  
  LOG_VERBOSE("change_mRNA_cytoplasm[%d]: mRNA=%d, transl rate=%g, protein decay=%g, salphc=%g\n", 
              i, state->mRNA_cyto_num[i], genes->translation[i], konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX], salphc);
  
  rates->salphc += konStates->kon_list[i]->site_count*kon*(salphc - konStates->konvalues[i][KON_SALPHC_INDEX]);
  rates->salphc_operations++;

  rates->max_salphc += konStates->kon_list[i]->site_count*kon*(fmaxf(state->protein_conc[i], salphc) - fmaxf(state->protein_conc[i], konStates->konvalues[i][KON_SALPHC_INDEX]));
  rates->max_salphc_operations++;
  rates->min_salphc += konStates->kon_list[i]->site_count*kon*(fminf(state->protein_conc[i], salphc) - fminf(state->protein_conc[i], konStates->konvalues[i][KON_SALPHC_INDEX]));    
  rates->min_salphc_operations++;

  konStates->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
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
  for (j=0; j < state->tf_bound_num; j++) {
    if (all_binding_sites[k].cisreg_id==all_binding_sites[state->tf_bound_indexes[j]].cisreg_id &&
        all_binding_sites[k].gene_copy==all_binding_sites[state->tf_bound_indexes[j]].gene_copy &&
        !(k==state->tf_bound_indexes[j])) {
      posdiff = all_binding_sites[k].left_edge_pos - all_binding_sites[state->tf_bound_indexes[j]].left_edge_pos;
      //printf("diff=%d\n", posdiff);
      if (abs(posdiff) < HIND_LENGTH) {/*Phey*/
        LOG_ERROR("t=%g steric hindrance breached with site %d (at pos %d, strand %d copy %d gene %d), %d away from site %d \
                   (at pos %d, strand %d copy %d of gene %d)\n",
                  t, k, all_binding_sites[k].left_edge_pos, all_binding_sites[k].strand, all_binding_sites[k].gene_copy, all_binding_sites[k].cisreg_id,
                  posdiff, state->tf_bound_indexes[j], all_binding_sites[state->tf_bound_indexes[j]].left_edge_pos, 
                  all_binding_sites[state->tf_bound_indexes[j]].strand, all_binding_sites[state->tf_bound_indexes[j]].gene_copy, 
                  all_binding_sites[state->tf_bound_indexes[j]].cisreg_id);
        exit(-1);
      }
      if (abs(posdiff) < cooperative_distance) {
        if (posdiff>0) front++; else back++;
      }
    }
    if ((front) && (back)) 
      j=state->tf_bound_num;
  }
  if (front>0) Gibbs -= cooperativity*state->RTlnKr/3;
  if (back>0) Gibbs -= cooperativity*state->RTlnKr/3;
  *koff = NumSitesInGenome*kon*0.25/exp(-Gibbs/(GasConstant*state->temperature));
  /*  fprintf(fperrors,"state->RTlnKr=%g front=%d back=%d H=%d Gibbs=%g koff=%g\n",state->RTlnKr,front,back,all_binding_sites[k][4],Gibbs,*koff); */
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
      //printf("diff3=%d\n", posdiff);
      if (abs(posdiff) < HIND_LENGTH) { /* within HIND_LENGTH: bad: shouldn't happen Phey*/
        LOG_ERROR("t=%g steric hindrance 2 has been breached with site %d %d away from site %d\n",
                  t, indexChanged, posdiff, state->tf_bound_indexes[j]);
      }
      if (abs(posdiff) < cooperative_distance) {  /* within 20, adjust koff */

        /* save old value */
        diff = -koffvalues[j];

        /* recompute koffvalues */
        calc_koff(state->tf_bound_indexes[j], all_binding_sites, state, &(koffvalues[j]), t);

        /* calculating how koff changes  */
        diff += koffvalues[j];

        /* adjust rates by difference */
        rates->koff += diff;
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
                KonStates *konStates,
                float protein_concTFID)
{
  int k;
  
  k = 0;
  
  while (!(konStates->kon_list[TFID]->available_sites[k] == site_id) && k < konStates->kon_list[TFID]->site_count) {
    k++;
  }

  LOG_VERBOSE_NOCELLID(">>> remove site %d kon_list (k=%d of %d total sites for TF %d, grandtotal=%d)\n", 
                       site_id, k, konStates->kon_list[TFID]->site_count, TFID, konStates->nkon);

  /* make sure that we have enough unoccupied sites left */
  if (k < konStates->kon_list[TFID]->site_count && k < konStates->nkon) { 
    /* adjust rates */
    rates->salphc -= kon*salphc;
    rates->salphc_operations++;

    rates->max_salphc -= kon*fmaxf(protein_concTFID, salphc);
    rates->max_salphc_operations++;

    rates->min_salphc -= kon*fminf(protein_concTFID, salphc);
    rates->min_salphc_operations++;

    /* one less site available for binding of total */
    (konStates->nkon)--;

    /* also per gene */
    (konStates->kon_list[TFID]->site_count)--;

    /* move the last element end of array into space vacated by site k */
    konStates->kon_list[TFID]->available_sites[k] = konStates->kon_list[TFID]->available_sites[konStates->kon_list[TFID]->site_count];
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
             KonStates *konStates)
{

  /* update rates because new site is now available */
  rates->salphc += kon*salphc;
  rates->salphc_operations++;

  rates->max_salphc += fmaxf(protein_concTFID, salphc);
  rates->max_salphc_operations++;

  rates->min_salphc += fminf(protein_concTFID, salphc);
  rates->min_salphc_operations++;

  /* add back site_id to pool of available sites */
  konStates->kon_list[TFID]->available_sites[konStates->kon_list[TFID]->site_count] = site_id;

  /* one more site available */
  (konStates->kon_list[TFID]->site_count)++;
  (konStates->nkon)++;
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
  int i,off;
  
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
void calc_from_state(Genotype *genes,
                     CellState *state,
                     GillespieRates *rates,
                     KonStates *konStates,
                     float transport[],
                     float mRNAdecay[])
{
  int i, j, k;
  float salphc, protein_concTFID; 
  float protein_decay;

  // TODO FIXME: this updates konStates for non-TFs, even though this isn't used
  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to fall below cut-off, use aging term */
    protein_decay = genes->proteindecay[i] >= protein_aging ? genes->proteindecay[i] : protein_aging;
    salphc = (float) (state->mRNA_cyto_num[i]) * genes->translation[i] / (protein_decay);
    konStates->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / (protein_decay);
    konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
    konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
    konStates->kon_list[i]->site_count = 0;
    LOG_VERBOSE("protein decay[%d]=%g\n", i, konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]);
  }
  state->tf_bound_num=0;

  rates->koff=0.0;
  rates->koff_operations = 0;

  rates->transport=0.0;
  rates->transport_operations = 0;

  rates->mRNAdecay=0.0;
  rates->mRNAdecay_operations = 0;

  rates->pic_disassembly=0.0;
  rates->pic_disassembly_operations = 0;

  rates->salphc=0.0;
  rates->salphc_operations=0;

  rates->max_salphc=0.0;
  rates->max_salphc_operations=0;

  rates->min_salphc=0.0;
  rates->min_salphc_operations=0;

  for (k=0; k < genes->binding_sites_num; k++) {
    i = genes->all_binding_sites[k].tf_id;
    protein_concTFID = state->protein_conc[i];
    salphc = konStates->konvalues[i][KON_SALPHC_INDEX];

    rates->salphc += salphc;
    rates->salphc_operations++;

    rates->max_salphc += fmaxf(protein_concTFID, salphc);
    rates->max_salphc_operations++;

    rates->min_salphc += fminf(protein_concTFID, salphc);
    rates->min_salphc_operations++;

    /* update the list of sites that bind for a particular TF, i */
    konStates->kon_list[i]->available_sites[konStates->kon_list[i]->site_count] = k;
    (konStates->kon_list[i]->site_count)++;
  }

  /* initialize konStates->nkon as the total number of binding sites */
  konStates->nkon = genes->binding_sites_num;

  for (i=0; i < NGENES; i++) {

    // FIXME: this updates konStates for non-TFs, even though this isn't used
    LOG("after initializing konStates for gene=%d, site_count=%d, sites_per_gene=%d, nkon=%d\n", 
        i, konStates->kon_list[i]->site_count, genes->sites_per_gene[i], konStates->nkon);

    transport[i] = kRNA * (float) (state->mRNA_nuclear_num[i]);
    rates->transport += transport[i];
    rates->transport_operations++;
    LOG_VERBOSE("mRNA_nuclear_num=%d, initializing transport[%d]=%g\n", state->mRNA_nuclear_num[i], i, transport[i]);
  }
  LOG_VERBOSE("initializing rates->transport=%g\n", rates->transport);

  /* start all genes in acteylated state */
  for (j=0; j < MAX_COPIES; j++) {  
    int pos = 0;
    for (i=0; i < NGENES; i++) {
      if (genes->copies[i] > j) {  
        state->state_change_ids[ACETYLATION][j][pos] = i;
        LOG_VERBOSE("Initializing statechange gene=%d, ploidy=%d state_change_ids[%d][%d]=%d\n", i, j, j, 
                    pos, state->state_change_ids[ACETYLATION][j][pos]);
        pos++;
      }
    }
  }

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
    for (j=0; j < genes->copies[i]; j++) {
      rates->acetylation_num[j]++;
    }
    
  }

  if (verbose) 
    for (j=0; j < MAX_COPIES; j++) {
      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
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
                         FixedEvent *replication_time_end,
                         float t)
{
  int retval;
  float t1;
  float t2;
  float t3;

  if (mRNA_transcr_time_end == NULL && mRNA_transl_time_end==NULL && replication_time_end==NULL) {
    retval = 0;
  } else {
    // TODO: rewrite this to avoid use of magic number
    t1 = mRNA_transcr_time_end ? mRNA_transcr_time_end->time : TIME_INFINITY;
    t2 = mRNA_transl_time_end ? mRNA_transl_time_end->time : TIME_INFINITY;
    t3 = replication_time_end ? replication_time_end->time : TIME_INFINITY;

    LOG_VERBOSE_NOCELLID("check fixed event: t1=%g, t2=%g, t3=%f [t=%g] ", t1, t2, t3, t);

    if ((t1 < t2) && (t1 < t3) && (t1 < t)) { 
      if (mRNA_transcr_time_end == NULL) retval = 0;
      else retval = 1;
    } else
      if ((t2 < t3) && (t2 < t1) && (t2 < t)) {
        if (mRNA_transl_time_end == NULL)  retval = 0;
        else retval = 2;
      }
      else 
      if ((t3 < t1) && (t3 < t2) && (t3 < t)) {
        if (replication_time_end == NULL) retval = 0;
        else retval = 3;
      } 
      else {
        retval = 0;
      }
  }
  LOG_VERBOSE_NOCELLID("event=%d\n", retval);
  return retval;  
}

void calc_dt(float *x,
             float *dt,
             GillespieRates *rates,
             KonStates *konStates,
             float mRNAdecay[],
             float mRNAdecayrates[],
             int mRNA_cyto_num[],
             int mRNA_transl_cyto_num[],
             int cell_id)
{
  float tbound1, tbound2;
  int i, j;

  /* reset the total rate for current step */
  rates->total=0.0;
  
  /* reset mRNA decay rate */
  rates->mRNAdecay=0.0;

  /* hence reset the number of rounding operations (TODO: check!) */
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
  rates->total += rates->koff;
  rates->total += rates->transport;
  rates->total += rates->mRNAdecay;
  rates->total += rates->pic_disassembly;
  rates->total += rates->salphc;

  /* 
   * convert the counts back into rates using the constants 
   */
  for (j=0; j < MAX_COPIES; j++) {
    rates->total += (float) rates->acetylation_num[j] * acetylate;
    rates->total += (float) rates->deacetylation_num[j] * deacetylate;
    rates->total += (float) rates->pic_assembly_num[j] * PICassembly;
    rates->total += (float) rates->transcript_init_num[j] * transcriptinit;    
  } 

  tbound1 = *x/(rates->total + rates->max_salphc);
  tbound2 = *x/(rates->total + rates->min_salphc);
  LOG_VERBOSE_NOCELLID("[cell %03d] bounds %g %g\n", cell_id, tbound1, tbound2);

  /* if bounds are the same, simply choose tbound1 */
  if (tbound1==tbound2){
    if (konStates->nkon!=0) {
      LOG_ERROR_NOCELLID("[cell %03d] nkon=%d when it should be zero x=%f rates->max_salphc=%g rates->min_salphc=%g rates->total=%g\n",
                         cell_id, konStates->nkon, *x, rates->max_salphc, rates->min_salphc, rates->total);
    }
    *dt = tbound1;
  } else {
    /* otherwise get delta t by solving the equation using Newton-Raphson method */
    *dt = rtsafe(&calc_time, *x, rates, konStates, tbound1, tbound2, (float) RT_SAFE_EPSILON); 
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

  /* add rate kRNA to transport and Gillespie rates */
  transport[i] += kRNA;
  rates->transport += kRNA;
  rates->transport_operations++;

  LOG_VERBOSE("add one new mRNA in nucleus, updating transport[%d]=%g, rates->transport=%g\n", i, transport[i], rates->transport);

}

// TODO: move to lib.c?  utility function
void remove_from_array(int toberemoved,
                       int type,
                       int a[],
                       int *len,
                       int force)
{
  int i;
  i = 0;

  /* check the range of i first so we don't access an uninitialized
     array position in 'a'  */
  while ((i < *len) && !(a[i]==toberemoved)) { 
    i++;
  }
  if (i < *len) {  
    (*len)--;
    a[i]=a[*len];
  }
  else 
    if (force)  {
      /* don't always print because with a 4->3 transition PIC assembly is not there to be removed */
      LOG_ERROR_NOCELLID("error removing %d from array of length %d, type=%d\n", toberemoved, *len, type);
    }
}

void disassemble_PIC(CellState *state,
                     Genotype *genes,
                     int gene_id,
                     int gene_copy,
                     GillespieRates *rates)
{
  float disassembly = genes->pic_disassembly[gene_id][gene_copy];
  remove_from_array(gene_id, TRANSCRIPTINIT, state->state_change_ids[TRANSCRIPTINIT][gene_copy], &(rates->transcript_init_num[gene_copy]), (int) 1);
  remove_from_array(gene_id, PICDISASSEMBLY, state->state_change_ids[PICDISASSEMBLY][gene_copy], &(rates->pic_disassembly_num[gene_copy]), (int) 1);
  rates->pic_disassembly -= disassembly;
  rates->pic_disassembly_operations++;
  
  /* disassemble PIC in OFF state */
  if (state->active[gene_id][gene_copy] == OFF_PIC) {
    (state->active[gene_id][gene_copy]) = OFF_NO_PIC;
    state->state_change_ids[DEACETYLATION][gene_copy][rates->deacetylation_num[gene_copy]] = gene_id;
    (rates->deacetylation_num[gene_copy])++;
  }
  /* disassemble PIC in ON state */
  if (state->active[gene_id][gene_copy] == ON_FULL) {
    (state->active[gene_id][gene_copy]) = ON_NO_PIC;
  }
}

void revise_activity_state(int gene_id,
                           int gene_copy,
                           Genotype *genes,
                           CellState *state,
                           GillespieRates *rates)
{
  int transcriptrule, oldstate, numactive;

  transcriptrule = ready_to_transcribe(gene_id, gene_copy, 
                                       state->tf_bound_indexes, 
                                       state->tf_bound_num,
                                       genes->all_binding_sites,
                                       genes->activating,
                                       &numactive);

  /* get last state of transcription initiation */
  oldstate = state->active[gene_id][gene_copy];

  /*
   * first set of rules:
   * ACTIVATING TFs exceed REPRESSING TFs 
   */

  /* OFF_FULL -> ON_WITH_NUCLEOSOME */
  if ((transcriptrule) && oldstate==OFF_FULL){
    state->active[gene_id][gene_copy] = ON_WITH_NUCLEOSOME;
    state->state_change_ids[ACETYLATION][gene_copy][rates->acetylation_num[gene_copy]] = gene_id;
    (rates->acetylation_num[gene_copy])++;
  }
  /* OFF_NO_PIC -> ON_NO_PIC */
  if ((transcriptrule) && oldstate==OFF_NO_PIC) {
    state->active[gene_id][gene_copy] = ON_NO_PIC;
    remove_from_array(gene_id, DEACETYLATION,  state->state_change_ids[DEACETYLATION][gene_copy], &(rates->deacetylation_num[gene_copy]), (int) 1);
    if (numactive){
      state->state_change_ids[PICASSEMBLY][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id;
      (rates->pic_assembly_num[gene_copy])++;
    }
  }
  /* OFF_PIC -> ON_FULL */
  if ((transcriptrule) && oldstate==OFF_PIC) {
    state->active[gene_id][gene_copy] = ON_FULL;
  }

  /*
   * second set of rules:
   * REPRESSING TFs exceed ACTIVATING TFs 
   */

  /* ON_WITH_NUCLEOSOME -> OFF_FULL */
  if (!(transcriptrule) && oldstate==ON_WITH_NUCLEOSOME) {
    state->active[gene_id][gene_copy] = OFF_FULL;
    LOG_VERBOSE("removing gene=%d, copy=%d from state_change_ids[ACETYLATION][%d]\n", gene_id, gene_copy, gene_copy);
    remove_from_array(gene_id, ACETYLATION, state->state_change_ids[ACETYLATION][gene_copy], &(rates->acetylation_num[gene_copy]), (int) 1);
  }
  
  /* ON_NO_PIC -> OFF_NO_PIC */
  if (!(transcriptrule) && oldstate==ON_NO_PIC){          
    state->active[gene_id][gene_copy] = OFF_NO_PIC;
    remove_from_array(gene_id, PICASSEMBLY, state->state_change_ids[PICASSEMBLY][gene_copy], &(rates->pic_assembly_num[gene_copy]), (int) 0);
    state->state_change_ids[DEACETYLATION][gene_copy][rates->deacetylation_num[gene_copy]] = gene_id;
    (rates->deacetylation_num[gene_copy])++;
  }
  /* ON_FULL -> OFF_PIC  */
  if (!(transcriptrule) && oldstate==ON_FULL) {
    state->active[gene_id][gene_copy] = OFF_PIC;
  }

  /* do remaining transitions:
   * OFF_PIC -> OFF_NO_PIC
   * ON_FULL -> ON_NO_PIC 
   */
  if ((state->active[gene_id][gene_copy]==OFF_PIC || state->active[gene_id][gene_copy]==ON_FULL) && numactive==0)
    disassemble_PIC(state, genes, gene_id, gene_copy, rates);

  if (verbose && (oldstate!=state->active[gene_id][gene_copy])) {
    LOG_VERBOSE("state change from %d to %d in gene %d, copy %d\n", oldstate, state->active[gene_id][gene_copy], gene_id, gene_copy);
  }
}

void remove_tf_binding(Genotype *genes,
                       CellState *state,
                       GillespieRates *rates,
                       KonStates *konStates,
                       int site,
                       float koffvalues[],
                       float t)
{
  int i, j, k, bound, site_id, gene_id, gene_copy;

  i = 0;

  /* given site 'site', look for the index in the list of bound sites */
  while ((state->tf_bound_indexes[i] != site) && (i < state->tf_bound_num)) 
    i++;
  if (i == state->tf_bound_num) {  /* couldn't find the site */
    LOG_ERROR("t=%g could not find site %d with %d possibilities\n Bound sites are\n",
              t, site, state->tf_bound_num);
    for (j = 0; j < state->tf_bound_num; j++)  {
      LOG_NOFUNC("%d\n", state->tf_bound_indexes[j]);
    }
  }
  else {
    j = 0;
    /* loop through the sterically hindered sites */
    while (j < state->tf_hindered_num) {

      /* check all sites hindered by binding to location 'site' */
      if (state->tf_hindered_indexes[j][1] == site) {
        k = bound = 0;

        /* is anything else hindering the same site? */
        while (bound == 0 && k < state->tf_hindered_num) {
          if (state->tf_hindered_indexes[j][0] == state->tf_hindered_indexes[k][0] && j != k) 
            bound=1;
          k++;
        }

        /* if nothing else is hindering this site then allow site to be (re-)bound */
        if (bound==0) {
          site_id = state->tf_hindered_indexes[j][0];
          LOG_VERBOSE("Site %d left_edge_pos %d on gene %d freed from steric hindrance\n",
                      site_id, genes->all_binding_sites[site_id].left_edge_pos, genes->all_binding_sites[site_id].cisreg_id);

          /* adjust rates by returning kon to pool */
          add_kon(state->protein_conc[genes->all_binding_sites[site_id].tf_id],
                  konStates->konvalues[genes->all_binding_sites[site_id].tf_id][KON_SALPHC_INDEX],
                  genes->all_binding_sites[site_id].tf_id,
                  site_id,
                  rates,
                  konStates);
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
    gene_id = genes->all_binding_sites[site].cisreg_id;
    gene_copy = genes->all_binding_sites[site].gene_copy;
    LOG_VERBOSE("Add site %d at left_edge_pos %d on gene %d copy %d freed by unbinding\n",
                site, genes->all_binding_sites[site].left_edge_pos, gene_id, gene_copy);

    /* adjust kon */
    add_kon(state->protein_conc[genes->all_binding_sites[site].tf_id],
            konStates->konvalues[genes->all_binding_sites[site].tf_id][KON_SALPHC_INDEX],
            genes->all_binding_sites[site].tf_id,
            site,
            rates,
            konStates);

    /* adjust the state of the gene */
    revise_activity_state(gene_id, gene_copy, genes, state, rates);

    /* when TF unbinds adjust the co-operativity at close sites */
    scan_nearby_sites(site, genes->all_binding_sites, state, rates, koffvalues, t);
  }
}

void attempt_tf_binding(Genotype *genes,
                        CellState *state,
                        GillespieRates *rates,
                        float **koffvalues,
                        KonStates *konStates,
                        int *maxbound2,
                        int *maxbound3,
                        int site,
                        float t)
{
  int gene_id, gene_copy, k, posdiff;

  LOG_VERBOSE("kon1 event at site %d out of %d possible, %d TFs previously bound binding_sites_num=%d\n",
              site, konStates->nkon, state->tf_bound_num, genes->binding_sites_num);

  /* if we have run out of space, double memory  */
  if (state->tf_bound_num >= *maxbound2){
    (*maxbound2) *= 2;

    state->tf_bound_indexes = realloc(state->tf_bound_indexes, (*maxbound2)*sizeof(int));

    /* do the copy */
    *koffvalues = realloc(*koffvalues, (*maxbound2)*sizeof(float));

    /* check return value */
    if (!state->tf_bound_indexes || !(*koffvalues)) {
      LOG_ERROR("memory allocation error resetting maxbound2=%d\n", *maxbound2);
      exit(1);
    }
  }

  /* append the site to end of indexes */
  state->tf_bound_indexes[state->tf_bound_num] = site;
  LOG_VERBOSE("remove site %3d on gene %2d (copy %d)\n", 
              site, genes->all_binding_sites[site].cisreg_id, genes->all_binding_sites[site].gene_copy);

  /* remove the site from the kon pool */
  remove_kon(site,
             genes->all_binding_sites[site].tf_id,
             rates, 
             konStates->konvalues[genes->all_binding_sites[site].tf_id][KON_SALPHC_INDEX],
             konStates,
             state->protein_conc[genes->all_binding_sites[site].tf_id]);

  /* recompute the koffvalues */
  calc_koff(site, genes->all_binding_sites, state, &((*koffvalues)[state->tf_bound_num]), t);

  LOG_VERBOSE("new koff = %g is number %d\n",
              (*koffvalues)[state->tf_bound_num], (state->tf_bound_num+1));

  /* adjust rates by adding the new koffvalue to rates->koff */
  rates->koff += (*koffvalues)[state->tf_bound_num];
  rates->koff_operations++;

  state->tf_bound_indexes[state->tf_bound_num] = site;

  /* increment number of bound TFs */
  (state->tf_bound_num)++;

  /* adjust co-operative binding in context of new TF */
  scan_nearby_sites(site, genes->all_binding_sites, state, rates, *koffvalues, t);

  /* get the gene that the TF is binding to */
  gene_id = genes->all_binding_sites[site].cisreg_id;

  /* get the copy that the TF is binding to */
  gene_copy = genes->all_binding_sites[site].gene_copy;
  
  /* update steric hindrance data structures */
  /* JM: this cycles over all sites, not just bound ones, in order to
     record redundancy in steric hindrance*/
  for (k = 0; k < genes->binding_sites_num; k++) {

    /* if we are on the same gene and not the same binding site */
    if (gene_id == genes->all_binding_sites[k].cisreg_id &&
        gene_copy == genes->all_binding_sites[k].gene_copy &&
        !(k==site)) {

      /* check distance from current binding site (k) to the original (site) */
      //fprintf(fperrors, "site=%d k=%d\n", genes->all_binding_sites[site].left_edge_pos, genes->all_binding_sites[k].left_edge_pos);
      posdiff = genes->all_binding_sites[site].left_edge_pos - genes->all_binding_sites[k].left_edge_pos;

      /* if within HIND_LENGTH, we prevent binding by adding to steric hindrance */
      if (abs(posdiff) < HIND_LENGTH) {/* Phey */
        /* if not enough memory, reallocate */
        if (state->tf_hindered_num > *maxbound3 - 1) {
          (*maxbound3) *= 2;
          state->tf_hindered_indexes = realloc(state->tf_hindered_indexes,2*(*maxbound3)*sizeof(int));
        }
        
        /* record hindrance: 'site' blocks 'k' */
        state->tf_hindered_indexes[state->tf_hindered_num][0] = k;
        state->tf_hindered_indexes[state->tf_hindered_num][1] = site;

        /* update list of hindered count */
        (state->tf_hindered_num)++;

        LOG_VERBOSE("%d steric hindrance sites after %d blocks site %d\n", state->tf_hindered_num, site, k);

        /* remove the kon from pool */
        remove_kon(k,
                   genes->all_binding_sites[k].tf_id,
                   rates,
                   konStates->konvalues[genes->all_binding_sites[k].tf_id][KON_SALPHC_INDEX],
                   konStates,
                   state->protein_conc[genes->all_binding_sites[k].tf_id]);
      }
    }
  }
  LOG_VERBOSE("tf_bound_num=%d tf_hindered_num=%d maxbound2=%d maxbound3=%d\n",
              state->tf_bound_num, state->tf_hindered_num, *maxbound2, *maxbound3);

  /* gene activity may change as a result of binding */
  revise_activity_state(gene_id, gene_copy, genes, state, rates);
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

void reach_s_phase(CellState *state, Genotype *genes, float t) {
  int i;

  /* set the state as being in S phase */
  state->in_s_phase = 1;

  // only add replication events if S phase and G2 phase have non-zero length
  // otherwise just jump to division
  if (time_s_phase + time_g2_phase > 0.0) {
    for (i=0; i < NGENES; i++) {
      LOG("add replication for gene_id=%d at t=%g\n", i, t + genes->replication_time[i]);

      /* add the replication event to the queue */
      add_fixed_event(i, -1, t + genes->replication_time[i], 
                      &(state->replication_time_end), &(state->replication_time_end_last));
    }
  }
}

float compute_tprime(float c, float P, float alpha, float s_mRNA) {
  return (1/c) * log((c*P - alpha*s_mRNA)/(c*P - alpha*s_mRNA));
}

float compute_integral(float alpha, float c, float gmax, float deltat, float s_mRNA, float P, float Pp, float ect, float ect1) {
  return 1.0/(pow(c,2)*Pp) * gmax * (-alpha*ect1*s_mRNA + c*(P*ect1 + alpha*deltat*s_mRNA));
}

float compute_growth_rate_dimer(float *integrated_growth_rate,
                                float alpha, 
                                float s_mRNA,
                                float all_alpha[NGENES],
                                int all_s_mRNA[NGENES],
                                float P,
                                float P_next,
                                float t, 
                                float deltat,
                                float c,
                                float ect,
                                float ect1,
                                int cell_id) 
{
  int i;
  float instantaneous_growth_rate;
  float total_alpha_s = 0.0;
  float deltatprime, deltatrest;

  LOG_VERBOSE_NOCELLID("P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
                       P, P_next, c, t, t+deltat, s_mRNA);

  /* choose the appropriate piecewise linear integral */
  if (((P > Pp) && (P_next >= P)) || ((P_next > Pp) && (P >= P_next))) {          /* P > Pp throughout */
    if (verbose)
      LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);
    *integrated_growth_rate = gmax * deltat;
  } else if (((P_next < Pp) && (P_next >= P)) || ((P < Pp) && (P >= P_next))) {   /* P < Pp throughout */
    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
    *integrated_growth_rate = compute_integral(alpha, c, gmax, deltat, s_mRNA, P, Pp, ect, ect1);
  } else if ((Pp > P) && (P_next > P)) {    /* P < Pp up until t' then P > Pp */
    deltatprime = compute_tprime(c, P, alpha, s_mRNA);
    deltatrest = deltat - deltatprime;
    LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
    *integrated_growth_rate = compute_integral(alpha, c, gmax, deltatprime, s_mRNA, P, Pp, ect, ect1);
    *integrated_growth_rate += gmax * deltatrest;
  } else if ((P > Pp) && (P > P_next)) {   /* P > Pp up until t' then P < Pp */
    deltatprime = compute_tprime(c, P, alpha, s_mRNA);
    deltatrest = deltat - deltatprime;
    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
    *integrated_growth_rate = gmax * deltatprime;
    *integrated_growth_rate += compute_integral(alpha, c, gmax, deltatrest, s_mRNA, P, Pp, ect, ect1);
  } else {
    LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
    LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
                       cell_id);
    exit(1);
  }

  /* compute instantaneous growth rate at t */
  if (P < Pp)
    instantaneous_growth_rate = gmax*P/Pp;
  else
    instantaneous_growth_rate = gmax;

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
  
  // TODO: current disable printing out growth rate information
#if 0
  fprintf(fp_growthrate[cell_id], "%g %g %g %g %g %g\n", t, instantaneous_growth_rate, *integrated_growth_rate, P, s_mRNA, c);
#endif

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
                                   Genotype *genes,
                                   float dt,
                                   GillespieRates *rates,
                                   KonStates *konStates,
                                   float t,
                                   TimeCourse **timecoursestart,
                                   TimeCourse **timecourselast,
                                   float otherdata[])
{
  int i;
  float ct, ect, ect1;
  float L, L_next;
  float instantaneous_growth_rate = 0.0;
  float integrated_growth_rate = 0.0;
  float adjusted_decay;

  rates->max_salphc = rates->min_salphc = 0.0;
  for (i=0; i < NPROTEINS; i++) {
    if (i == SELECTION_GENE)  /* if we are looking at the selection gene, record protein concentration before update */
      L = protein_conc[i];

    /* update protein decay rates due to dilution caused by growth */
    adjusted_decay = genes->proteindecay[i] + state->growth_rate;

    /* if this results in a very small or zero decay rate, use protein aging term */
    // TODO: check 
    //if (adjusted_decay > 0.0)
    if (adjusted_decay > protein_aging)
      konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = adjusted_decay;
    else 
      konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = protein_aging;

    if (konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] < 1e-10) {
      LOG_WARNING("protein=%02d, protein_decay=%g, genes->proteindecay=%g, protein_aging=%g\n", i, adjusted_decay, genes->proteindecay[i], protein_aging);
    }

    LOG_VERBOSE("prot decay[%d]=%g\n", i, konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]);

    ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]*dt;
    ect = exp(-ct);
    if (fabs(ct)<EPSILON) ect1=ct;
    else ect1 = 1-ect;   
    protein_conc[i] = konStates->konvalues[i][KON_SALPHC_INDEX]*ect1 + ect*protein_conc[i];

    konStates->konvalues[i][KON_DIFF_INDEX] = (protein_conc[i] - konStates->konvalues[i][KON_SALPHC_INDEX]) / konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX];
    rates->max_salphc += ((float) konStates->kon_list[i]->site_count) * fmaxf(protein_conc[i], konStates->konvalues[i][KON_SALPHC_INDEX]);
    rates->min_salphc += ((float) konStates->kon_list[i]->site_count) * fminf(protein_conc[i], konStates->konvalues[i][KON_SALPHC_INDEX]);
    
    if (i == SELECTION_GENE) {  /* now find out the protein concentration at end of dt interval and compute growth rate */
      L_next = protein_conc[i];
      instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genes->translation[SELECTION_GENE], state->mRNA_cyto_num[SELECTION_GENE], 
                                                            genes->translation, state->mRNA_cyto_num, 
                                                            L, L_next, t, dt, 
                                                            konStates->konvalues[SELECTION_GENE][KON_PROTEIN_DECAY_INDEX], ect, ect1,
                                                            state->cell_id);
      state->cell_size = (state->cell_size)*exp(integrated_growth_rate);
      fprintf(fp_cellsize[state->cell_id], "%g %g\n", t, state->cell_size);
    }

  }
  rates->max_salphc *= kon;
  rates->min_salphc *= kon;
  if ((output) && (*timecourselast)->time < t+dt-0.1) 
    add_time_points(t+dt, otherdata, timecoursestart, timecourselast);

  /* update the growth rate for next timestep */
  state->growth_rate = instantaneous_growth_rate;
}

// TODO: this is diagnostic
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

  LOG_VERBOSE_NOCELLID("%d bound %g expected\n", tf_bound_num, (CISREG_LEN*TFGENES*sum)/NumSitesInGenome);
}

int sum_rate_counts(int rate_array[MAX_COPIES])
{
  int i;
  float retval = 0.0;

  for (i = 0; i < MAX_COPIES; i++) {
    retval += rate_array[i];
  }
  return retval;
}

void get_gene(int rate_array[MAX_COPIES], int pos, int *geneLoc, int *gene_copy)
{
  int i = 0;
  int total_rate = 0;
  *gene_copy = -1;   /* haven't found the copy yet */

  while (i < MAX_COPIES && *gene_copy < 0) {
    if (pos < (total_rate + rate_array[i])) {
      *gene_copy = i;
      *geneLoc = pos - total_rate;
    } else {
      total_rate += rate_array[i];
      i++;
    }
  }
}


/* -----------------------------------------------------
 * START
 * Functions that handle each possible Gillespie event 
 * ----------------------------------------------------- */
void transport_event(GillespieRates *rates,
                     CellState *state,
                     Genotype *genes,
                     KonStates *konStates,
                     float transport[NGENES],
                     TimeCourse **timecoursestart, 
                     TimeCourse **timecourselast, 
                     float tttranslation,
                     float dt,
                     float t,
                     float x)
{
  int i;
  float konrate2;
  float endtime = t+dt+ttranslation;

  update_protein_conc_cell_size(state->protein_conc, state, genes, dt, 
                                rates, konStates, 
                                t, timecoursestart, timecourselast, 
                                state->protein_conc);
  
  i = -1;
  konrate2 = 0.0;  

  /* choose gene product (mRNA) that gets transported to cytoplasm
     based on weighting in transport[] array */
  while (i < NGENES && x > konrate2) {
    i++;
    x -= transport[i];
  }

  if (i >= NGENES) {
    LOG_ERROR("[cell %03d] attempted to choose mRNA for gene=%d which doesn't exist\n", state->cell_id, i);
    exit(0);
  } 
  
  LOG_VERBOSE("do transport event mRNA from gene=%d from %d copies (x=%g)\n", i, state->mRNA_nuclear_num[i], x);

  /* one less mRNAs in nucleus */
  (state->mRNA_nuclear_num[i])--;

  /* it has just arrived in cytoplasm, ready to be translated */
  (state->mRNA_transl_cyto_num[i])++;
  
  /* add the endtime for translation */
  LOG_VERBOSE("add translation event endtime=%f for mRNA encoded by gene=%d \n", endtime, i);
  //add_fixed_event_end(i, endtime, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
  add_fixed_event_end(i, -1, endtime, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));

  /* decrease transport frequency */
  transport[i] -= kRNA;

  /* if below a threshold, make zero */
  if (transport[i] < 0.1*kRNA) 
    transport[i]=0.0;

  /* adjust rates */
  rates->transport -= kRNA;
  rates->transport_operations++;

  /* do similar threshold check */
  if (rates->transport < 0.1*kRNA) 
    rates->transport=0.0;
}


void tf_binding_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                      KonStates *konStates, float *koffvalues, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float konrate, float dt, float t, int maxbound2, int maxbound3, int update_protein_flag)
{
  float x = ran1(&seed) * (rates->salphc + konrate)/kon;
  int k, l = -1;  /* new */
  float total_konrate2, konrate2_for_TF = 0.0;     
  int site_id = -1;

  /* loop through all TFs, then choose a particular binding site */
  for (k=0; k < TFGENES; k++) {

    /* if no sites available for this TF, skip to next gene */
    if (konStates->kon_list[k]->site_count == 0) {
      LOG_VERBOSE("looking at TF: %d, has %d sites available, skipping\n", k, konStates->kon_list[k]->site_count);
      continue;
    }

    /* compute the total rate for all available binding sites for this
     * particular TF: see if we are in the right range
     */

    /* TODO: commented-out code that may help fix numerical issues with 1-exp(), but needs further testing */
    /* float c = konStates->konvalues[k][KON_PROTEIN_DECAY_INDEX];
     float ectdt;
     if (fabs(c*dt)<EPSILON) ectdt=c;
     else ectdt = (1-exp(-c*dt))/dt;    
     konrate2_for_TF = konStates->konvalues[k][KON_SALPHC_INDEX] + konStates->konvalues[k][KON_DIFF_INDEX] * ectdt;  */

    /* first, cache the konrate2 for this particular gene */
    konrate2_for_TF = konStates->konvalues[k][KON_SALPHC_INDEX] + 
       konStates->konvalues[k][KON_DIFF_INDEX] * (1-exp(-konStates->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt))/dt;  

    LOG_VERBOSE("TF:%d [KON_SALPHC: %g, KON_DIFF: %g, KON_PROTEIN_DECAY: %g]\nTF:%d [1-exp(-ct): %g, (1-exp(-ct)/dt): %g, konrate2_for_TF: %g]\n", 
                k,
                konStates->konvalues[k][KON_SALPHC_INDEX],
                konStates->konvalues[k][KON_DIFF_INDEX],
                konStates->konvalues[k][KON_PROTEIN_DECAY_INDEX],
                k,
                1-exp(-konStates->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt),
                (1-exp(-konStates->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt))/dt,
                konrate2_for_TF); 

    /* compute the *total* kon rate for all unbound sites for this TF  */
    total_konrate2 = ((konStates->kon_list[k]->site_count)) * konrate2_for_TF;

    LOG_VERBOSE("looking at TF: %d, has %d sites available [konrate2: %g, total_konrate2: %g, x: %g]\n", 
                k, konStates->kon_list[k]->site_count, konrate2_for_TF, total_konrate2, x); 

    /* if we are already in the appropriate TF, now choose a binding site */
    if (!(x > total_konrate2) || (k == TFGENES - 1)) {
      float konrate2 = 0.0;
      
      LOG_VERBOSE("selecting TF: %d, konrate2: %g, total_konrate2: %g, x: %g\n", k, konrate2_for_TF, total_konrate2, x); 
      
      while (l < (konStates->kon_list[k]->site_count - 1) && x > konrate2) {
        /* this will record the last binding site before we
           reach the appropriate binding site  */
        l++;
        
        /* get ID of site */
        site_id = konStates->kon_list[k]->available_sites[l];
        
        LOG_VERBOSE("l: %d, site: %d, binds to TF: %d, x = %g (site_count=%d)\n", l, site_id, k, x, konStates->kon_list[k]->site_count); 
        
        /* adjust random number */
        konrate2 = konrate2_for_TF;
        x -= konrate2;
      }
      /* found it, so break out of the outer for loop */
      break;
    } else {
      x -= total_konrate2; 
    }
    
  }
  
  /* print error if site not found */
  if (site_id == -1) {
    LOG_ERROR("no binding site could be found  TF: total_konrate2: %g, x: %g\n", total_konrate2, x);
  }
  else {
    LOG_VERBOSE("found a binding site l: %d, site: %d, binds to TF: %d, konrate2: %g, x: %g\n", l, site_id, k, konrate2_for_TF, x);  
  }

  if (update_protein_flag) 
    /* update protein concentration before doing the binding */
    update_protein_conc_cell_size(state->protein_conc, state, genes, dt, 
                                  rates, konStates, 
                                  t, timecoursestart, timecourselast,
                                  state->protein_conc);
  
  /* bind site_id, only if found */
  if (site_id != -1)
    attempt_tf_binding(genes, state, rates, &koffvalues, konStates, &maxbound2, &maxbound3, site_id, t);

  /* calculate the number of TFs bound */
  calc_num_bound(state->protein_conc, state->tf_bound_num);
}

void tf_unbinding_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                        KonStates *konStates, float *koffvalues, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                        float konrate, float dt, float t, float x, int update_protein_flag)
{
  int i, j = -1;
  int site;

  /* if an attempt is made to do an unbinding event when there is
     nothing bound, then we ignore this event and recompute kon and
     koff */
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
    recompute_koff_rates(rates, state, genes, koffvalues, t);    
    recompute_kon_rates(rates, state, genes, konStates, 0);

    return;  // TODO: should we throw out this dt? (check)
  }

  while (j < state->tf_bound_num && x > 0) {
    j++;
    x -= koffvalues[j];
  }
  if (j==state->tf_bound_num) {
    float konrate2 = 0.0;
    for (i = 0; i < state->tf_bound_num; i++) 
      konrate2 += koffvalues[i];
    LOG_WARNING("t=%g koffvalues add up to %g instead of rates->koff=%g\n",
                t, konrate2, rates->koff);
    rates->koff = konrate2;
    j--; /* a bit of a fudge for rounding error, really should move on to rates->transport, but too complicated for something so minor */
  } 
  site = state->tf_bound_indexes[j];
  LOG_VERBOSE("t=%g koff event %d of %d at site %d\n", t, j, state->tf_bound_num,site);
  if (j < 0) { LOG_ERROR("t=%g (x=%g) koff event %d of %d at site %d\n", t, x, j, state->tf_bound_num,site); exit(-1); }
  
  if (update_protein_flag) 
    /* update protein concentration before removing TF */
    update_protein_conc_cell_size(state->protein_conc, state, genes, dt, 
                                  rates, konStates, 
                                  t, timecoursestart, timecourselast, 
                                  state->protein_conc);
  
  /* remove TF binding from 'site' */
  remove_tf_binding(genes, state, rates, konStates, site, koffvalues, t);
  calc_num_bound(state->protein_conc, state->tf_bound_num);
}

void mRNA_decay_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                      KonStates *konStates, float *mRNAdecay, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t, float x)
{
  int i = -1, j;
  float konrate2 = 0.0;

  /* loop through mRNA products, to choose the mRNA with the
     proportionally higher decay rate */
  while (i < NGENES-1 && x > konrate2) {
    i++;
    konrate2 += mRNAdecay[i];
  }
  if (x > konrate2) { /* JM: had some rounding errors with rates->mRNAdecay. Calculate in calc_dt, hopefully fixed now */
    LOG_WARNING("x=%g > konrate2=%g out of rates->mRNAdecay=%g\n",
                x, konrate2, rates->mRNAdecay);
  }

  /* assume mRNA cytoplasm transport events equally likely */
  x = ran1(&seed)*((float) (state->mRNA_cyto_num[i] + state->mRNA_transl_cyto_num[i]));

  update_protein_conc_cell_size(state->protein_conc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  /* 
   * decay mRNA in cytoplasm 
   */
  if (x < (float)state->mRNA_cyto_num[i]) {
    LOG_VERBOSE("mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
                i, state->mRNA_cyto_num[i], state->mRNA_transl_cyto_num[i]);
    
    /* remove the mRNA from the cytoplasm count */
    (state->mRNA_cyto_num[i])--;  
    change_mRNA_cytoplasm(i, genes, state, rates, konStates); 
    
  } else {
    /* 
     * decay mRNA in process of translating
     */
    x = ran1(&seed)*((float) state->mRNA_transl_cyto_num[i]);
    LOG_VERBOSE("mRNA decay event gene %d not from %d copies in cytoplasm but %f from %d copies translating\n",
                i, state->mRNA_cyto_num[i], trunc(x), state->mRNA_transl_cyto_num[i]);
    
    /* delete this fixed event: this mRNA will never be translated */
    LOG_VERBOSE("delete fixed TRANSLATION EVENT at time =%d for gene=%d\n", (int) trunc(x), i);
    //delete_fixed_event(i, (int) trunc(x), &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
    delete_fixed_event(i, -1, (int) trunc(x), &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
    
    /* remove the mRNA from the count */
    (state->mRNA_transl_cyto_num[i])--;
    if (verbose) 
      for (j=0; j < NGENES; j++) {
        LOG_VERBOSE("%d copies of gene %d translating\n", state->mRNA_transl_cyto_num[j], j);
      }
  }
}

void histone_acteylation_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                               KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                               float dt, float t)
{
  int geneLoc, gene_copy;
  float x = ran1(&seed)*((float) sum_rate_counts(rates->acetylation_num));
  
  get_gene(rates->acetylation_num, (int)trunc(x), &geneLoc, &gene_copy);

  /* choose a particular gene to change state */
  int gene_id = state->state_change_ids[ACETYLATION][gene_copy][geneLoc];

  LOG_VERBOSE("acetylation event gene %d (copy %d)\nstate change from %d to 4\n",
              gene_id, gene_copy, state->active[gene_id][gene_copy]);
  if (state->active[gene_id][gene_copy] != ON_WITH_NUCLEOSOME) {
    LOG_ERROR("acetylation event on gene %d (copy %d) attempted from state %d\n", geneLoc, gene_copy, state->active[gene_id][gene_copy]);
  }

  /* update protein concentration and cell size */
  update_protein_conc_cell_size(state->protein_conc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  
  /* set state: eject nucleosome, but there is no PIC yet */
  state->active[gene_id][gene_copy] = ON_NO_PIC;
  remove_from_array(gene_id, ACETYLATION, state->state_change_ids[ACETYLATION][gene_copy], &(rates->acetylation_num[gene_copy]), (int) 1);
  if (is_one_activator(gene_id, gene_copy, state->tf_bound_indexes, state->tf_bound_num, 
                       genes->all_binding_sites, genes->activating)) {
    state->state_change_ids[PICASSEMBLY][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id; 
    (rates->pic_assembly_num[gene_copy])++;
  }
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                                 KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                                 float dt, float t)
{
  float x = ran1(&seed)*((float) sum_rate_counts(rates->deacetylation_num));

  int gene_copy; 
  int geneLoc; 

  get_gene(rates->deacetylation_num, (int)trunc(x), &geneLoc, &gene_copy);

  /* choose a particular gene and copy to change state */
  int gene_id = state->state_change_ids[DEACETYLATION][gene_copy][geneLoc];

  LOG_VERBOSE("deacetylation event gene %d (copy %d)\nstate change from %d to 1\n",
              gene_id, gene_copy, state->active[gene_id][gene_copy]);
  if (state->active[gene_id][gene_copy] != OFF_NO_PIC) {
    LOG_ERROR("deacetylation event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

  update_protein_conc_cell_size(state->protein_conc, state, genes, dt, 
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  /* set state: nucleosome returns */
  state->active[gene_id][gene_copy] = OFF_FULL;
  remove_from_array(gene_id, DEACETYLATION, state->state_change_ids[DEACETYLATION][gene_copy], &(rates->deacetylation_num[gene_copy]), (int) 1);
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                      KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t)
{

  float x = ran1(&seed)*((float) sum_rate_counts(rates->pic_assembly_num));

  int gene_copy; 
  int geneLoc; 

  get_gene(rates->pic_assembly_num, (int)trunc(x), &geneLoc, &gene_copy);

  /* choose a particular gene and copy to change state */
  int gene_id = state->state_change_ids[PICASSEMBLY][gene_copy][geneLoc];

  LOG_VERBOSE("PIC assembly event gene %d copy %d\nstate change from %d to 6\n",
              gene_id, gene_copy, state->active[gene_id][gene_copy]);

  if (state->active[gene_id][gene_copy] != ON_NO_PIC) {
    LOG_ERROR("PIC assembly event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

  update_protein_conc_cell_size(state->protein_conc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  
  /* turn gene fully on: ready for transcription */
  state->active[gene_id][gene_copy] = ON_FULL;
  remove_from_array(gene_id, PICASSEMBLY, state->state_change_ids[PICASSEMBLY][gene_copy], &(rates->pic_assembly_num[gene_copy]), (int) 1);
  state->state_change_ids[TRANSCRIPTINIT][gene_copy][rates->transcript_init_num[gene_copy]] = gene_id;

  (rates->transcript_init_num[gene_copy])++;

  state->state_change_ids[PICDISASSEMBLY][gene_copy][rates->pic_disassembly_num[gene_copy]] = gene_id;

  (rates->pic_disassembly_num[gene_copy])++;

  rates->pic_disassembly += genes->pic_disassembly[gene_id][gene_copy];
  rates->pic_disassembly_operations++;
}

void disassemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                           KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                           float dt, float t, float x)
{
  int gene_copy, geneLoc, gene_id;
  int j=-1;
  while (j < NGENES*current_ploidy && x>0) {
    j++;

    get_gene(rates->pic_disassembly_num, j, &geneLoc, &gene_copy);

    x -= genes->pic_disassembly[state->state_change_ids[PICDISASSEMBLY][gene_copy][geneLoc]][gene_copy];
  }
  if (j==NGENES*current_ploidy) { LOG_ERROR("error in PIC disassembly\n"); }
  gene_id = state->state_change_ids[PICDISASSEMBLY][gene_copy][geneLoc];
  LOG_VERBOSE("PIC disassembly event in copy %d of gene %d\n", gene_copy, gene_id);
  disassemble_PIC(state, genes, gene_id, gene_copy, rates);
}

void transcription_init_event(GillespieRates *rates, CellState *state, Genotype *genes,
                              KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                              float dt, float t, float x)
{
  int gene_id;

  x /= transcriptinit;

  int gene_copy; 
  int geneLoc; 

  /* choose the gene and copy that gets transcribed */
  get_gene(rates->transcript_init_num, (int)trunc(x), &geneLoc, &gene_copy);
  gene_id = state->state_change_ids[TRANSCRIPTINIT][gene_copy][geneLoc];
  LOG_VERBOSE("transcription event gene %d, copy %d\n", gene_id, gene_copy);

  if (state->active[gene_id][gene_copy] != ON_FULL && state->active[gene_id][gene_copy] != OFF_PIC) {
    LOG_ERROR("transcription event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

  update_protein_conc_cell_size(state->protein_conc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);

  /* now that transcription of gene has been initiated, 
   * we add the time it will end transcription, 
   * which is dt+time of transcription from now */
  //add_fixed_event_end(gene_id, t+dt+ttranscription, 
  //                    &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));
  add_fixed_event_end(gene_id, gene_copy, t+dt+ttranscription, 
                      &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));

  /* increase the number mRNAs being transcribed */
  (state->mRNA_transcr_num[gene_id][gene_copy])++;                      
}
/* -----------------------------------------------------
 * END
 * Functions that handle each possible Gillespie event 
 * ----------------------------------------------------- */

/* Helper function that shifts site_ids in KonStates, tf_hindered_indexes
   and tf_bound_indexes after a position by the specified offset */
void shift_binding_site_ids(CellState *state, 
                            KonStates *konStates,
                            int end,
                            int offset)
{
  int i, j, k, site_id;

  /* shift all sites in kon_list */
  for (i=0; i < TFGENES; i++) {
    k = 0;
    while (k < konStates->kon_list[i]->site_count) {
      site_id = konStates->kon_list[i]->available_sites[k];
      if (site_id >= end)
        konStates->kon_list[i]->available_sites[k] += offset;
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

/* eject TFs and replicate DNA */
void replicate_gene(CellState *state,
                    Genotype *genes,
                    GillespieRates *rates,
                    KonStates *konStates,
                    float *koffvalues,
                    int gene_id,
                    float t) 
{
  int i, k, l, p;
  int start_tfbs_pos, end_tfbs_pos, number_tfbs_pre_replication, offset;
  
  LOG("[gene %2d] duplicating at t=%g\n", gene_id, t);
  
  /* double the number of copies of gene being replicating */
  genes->copies[gene_id] = 2*current_ploidy;

  LOG("[gene %2d] before removing all TFs: site_count=%d, sites_per_gene=%d, nkon=%d\n", 
      gene_id, konStates->kon_list[gene_id]->site_count, genes->sites_per_gene[gene_id],
      konStates->nkon);
  LOG("[gene %2d] before removing all TFs: tf_bound_num=%d, tf_hindered_num=%d\n", gene_id, state->tf_bound_num, state->tf_hindered_num);
  
  /* first eject all TFs from this particular gene */
  
  int tfCount = state->tf_bound_num;
  for (k=0; k < tfCount; k++) {
    i = state->tf_bound_indexes[k];            /* first get the binding site ID */
    l = genes->all_binding_sites[i].cisreg_id;  /* now get the gene it belongs to */
    p = genes->all_binding_sites[i].gene_copy;
    
    // if we looking at the gene in question
    if (l == gene_id) {
      /* remove TF binding and update rates */
      LOG_VERBOSE("[gene %2d] ejecting TF on binding site_id=%d on copy=%d\n", i, l, p);
      remove_tf_binding(genes, state, rates, konStates, i, koffvalues, t);
    }
  }

  /* remove all PICs on that gene */
  for (p=0; p < current_ploidy; p++) 
    if ((state->active[gene_id][p]==OFF_PIC || state->active[gene_id][p]==ON_FULL))
      disassemble_PIC(state, genes, gene_id, p, rates);
  
  LOG_VERBOSE("[gene %2d] number of binding sites before adding new sites=%d at t=%g\n", 
              gene_id, genes->binding_sites_num, t);
  
  /* do mutation */
  // TODO: make mutation rate a parameter 
  for (p=0; p<2*current_ploidy; p++)
    mutate(genes, gene_id, p, 0.01);  

  /* record number of TFBS pre-replication */
  number_tfbs_pre_replication = genes->sites_per_gene[gene_id];
  start_tfbs_pos = 0;
  end_tfbs_pos = 0;

  /* record the beginning and end site_ids of the pre-replication list
     of binding sites  */
  for (i=0; i<=gene_id; i++) {  
    start_tfbs_pos = end_tfbs_pos;             
    end_tfbs_pos += genes->sites_per_gene[i];
  }
  end_tfbs_pos--;  /* always one less than the end point*/

  LOG("[gene %2d] has %d TFBS before replication [run from %d to %d]\n", 
      gene_id, number_tfbs_pre_replication, start_tfbs_pos, end_tfbs_pos);

  LOG("[gene %2d] after removing all TFs: site_count=%d, sites_per_gene=%d, nkon=%d\n", 
      gene_id, konStates->kon_list[gene_id]->site_count, genes->sites_per_gene[gene_id],
      konStates->nkon);

  LOG("[gene %2d] after removing all TFs: tf_bound_num=%d, tf_hindered_num=%d\n", gene_id, state->tf_bound_num, state->tf_hindered_num);  


  /* remove all of these old TFBS from konStates, some of them may no
     longer exist after mutation, we re-add them with add_kon() call
     after new binding sites computed */
  for (k=start_tfbs_pos; k < end_tfbs_pos + 1; k++) {
    remove_kon(k,
               genes->all_binding_sites[k].tf_id,
               rates, 
               konStates->konvalues[genes->all_binding_sites[k].tf_id][KON_SALPHC_INDEX],
               konStates,
               state->protein_conc[genes->all_binding_sites[k].tf_id]);
  }
  
  /* recompute *all* binding sites, then relabel sites offset by
     insertion (or deletion) of new sites created by replication */
  calc_all_binding_sites(genes->copies, 
                         genes->cisreg_seq, 
                         genes->tf_seq, 
                         &(genes->binding_sites_num),
                         &(genes->all_binding_sites),
                         genes->hindrance_positions,
                         genes->sites_per_gene,
                         genes->site_id_pos); 

  /* print_all_binding_sites(genes->copies, genes->all_binding_sites, genes->binding_sites_num, 
     genes->tf_seq, genes->cisreg_seq);  */

  /* use new sites_per_gene and pre-replication number to compute the
     offset to shift the site_ids */
  offset = genes->sites_per_gene[gene_id] - number_tfbs_pre_replication;

  LOG_VERBOSE("[gene %2d] has %d TFBS after replication [run from %d to %d]\n", 
              gene_id, genes->sites_per_gene[gene_id], start_tfbs_pos, end_tfbs_pos + offset);
  LOG_VERBOSE(" shift all TFBS starting at %d by %d\n", end_tfbs_pos + 1, offset);
  LOG_VERBOSE(" number of binding sites after adding new sites=%d at t=%g\n", genes->binding_sites_num, t);

  /* starting at the original ending point, move all sites along by
     'offset'.  Note this assumes that TFBS for a particular gene are
     always stored contiguously. */
  shift_binding_site_ids(state, konStates, end_tfbs_pos + 1, offset);
  
  /* update the konStates data structure to make available the newly
   * created TF binding sites in the full [start, end+offset] region
   */
  LOG("[gene %2d] adding new unbound sites from=%d to=%d\n", gene_id, start_tfbs_pos, (end_tfbs_pos + offset));

  for (k=start_tfbs_pos; k <= end_tfbs_pos + offset; k++) {
    add_kon(state->protein_conc[genes->all_binding_sites[k].tf_id],
            konStates->konvalues[genes->all_binding_sites[k].tf_id][KON_SALPHC_INDEX],
            genes->all_binding_sites[k].tf_id,
            k,
            rates,
            konStates);
  }
  
  for (i=0; i < current_ploidy; i++) {
    p = i + current_ploidy;
    
    /* set acetylation state in new gene copy  */
    state->state_change_ids[ACETYLATION][p][rates->acetylation_num[p]] = gene_id;
    
    /* update the counts for the acetylation */
    rates->acetylation_num[p]++;
    
    LOG_VERBOSE("[gene %2d] [clone acetylation]: ploidy=%d state_change_ids[%d][%d]=%d\n", gene_id, p, p, 
                rates->acetylation_num[p], state->state_change_ids[ACETYLATION][p][rates->acetylation_num[p]]);
    LOG_VERBOSE(" rates->acetylation_num[%d]=%d\n", p, rates->acetylation_num[p]);
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
                         Genotype *genes,
                         KonStates *konStates,
                         int recalibrate) 
{
  int j, k;
  float salphc;
  float orig_rates = rates->salphc;

  LOG("BEFORE rates->total=%g, rates->salphc=%g\n",   rates->total, rates->salphc);
  // subtract off current salphc from total
  rates->total -= rates->salphc;

  // reset rates
  rates->salphc=0.0;
  rates->salphc_operations=0;

  rates->max_salphc=0.0;
  rates->max_salphc_operations=0;

  rates->min_salphc=0.0;
  rates->min_salphc_operations=0;

  // go over all the currently unbound and increment salphc
  for (k=0; k < genes->binding_sites_num; k++) {
    int site_id;
    int tf_id = genes->all_binding_sites[k].tf_id;
    
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

      salphc = konStates->konvalues[tf_id][KON_SALPHC_INDEX];

      LOG_VERBOSE("for unbound site %03d [tf_id=%02d] adding salphc=%g to rates->salphc=%g\n", k, tf_id, salphc, rates->salphc);

      rates->salphc += salphc;
      rates->salphc_operations++;

      rates->max_salphc += fmaxf(state->protein_conc[tf_id], salphc);
      rates->max_salphc_operations++;

      rates->min_salphc += fminf(state->protein_conc[tf_id], salphc);
      rates->min_salphc_operations++;

      // TODO: check
      // only do if recalibrating state of cell from scratch
      if (recalibrate) {
        // update the list of sites that bind for a particular TF, i
        konStates->kon_list[tf_id]->available_sites[konStates->kon_list[tf_id]->site_count] = k;
        (konStates->kon_list[tf_id]->site_count)++;
        konStates->nkon++;
      }
    } 
  }
  rates->salphc = kon*rates->salphc;
  rates->max_salphc = kon*rates->max_salphc;
  rates->min_salphc = kon*rates->min_salphc;

  /* now that it is recomputed, add salphc back to total */
  rates->total += rates->salphc;
  LOG("AFTER rates->total=%g, rates->salphc=%g, percent difference=%g\n",   
      rates->total, rates->salphc, 100.0*(fabs(rates->salphc-orig_rates)/rates->salphc));
}

/*
 * recomputes all koff rates from the full cell state
 */
void recompute_koff_rates(GillespieRates *rates,
                          CellState *state,
                          Genotype *genes,
                          float *koffvalues,
                          float t) 
{
  int i;
  float orig_rates = rates->koff;
  float temprate = 0.0;

  LOG("BEFORE rates->total=%g, rates->koff=%g\n",   rates->total, rates->koff);
  /* subtract off current koff from total */
  rates->total -= rates->koff;
  
  for (i = 0; i < state->tf_bound_num; i++) {
    int site = state->tf_bound_indexes[i];
    calc_koff(site, genes->all_binding_sites, state, &(koffvalues[i]), t);
    temprate += koffvalues[i];
  }

#if 0
  fprintf(fp_koff[state->cell_id], "%g %g\n", t, rates->koff - temprate);
#endif

  rates->koff = temprate;
  rates->koff_operations = 0;

  /* add back new koff */
  rates->total += rates->koff;
  LOG("AFTER rates->total=%g, rates->koff=%g, percent difference=%g\n",   
      rates->total, rates->koff, 100.0*(fabs(rates->koff-orig_rates)/rates->koff));
}

/*
 * recalibrate the rates and cached data structures (KonStates) from
 * the current cell state
 */
void recalibrate_cell(GillespieRates *rates,
                      CellState *state,
                      Genotype *genes,
                      KonStates *konStates,
                      float **koffvalues,
                      float mRNAdecay[NGENES],
                      float transport[NGENES],
                      float dt) 
{
  int i, j; //, k;
  float protein_decay;
  float salphc = 0.0;
  //int site_id, tf_id;

  /* reset the total rate for current step */
  rates->total=0.0;
  
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

  /* regenerate konStates */
  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to fall below cut-off, use aging term */
    protein_decay = genes->proteindecay[i] >= protein_aging ? genes->proteindecay[i] : protein_aging;
    salphc = (float) (state->mRNA_cyto_num[i]) * genes->translation[i] / (protein_decay);
    konStates->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / (protein_decay);
    konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
    konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
    konStates->kon_list[i]->site_count = 0;
  }
  konStates->nkon = 0;   /* initialize to zero */

  /* regenerate konStates and rates->{salphc,max_salphc,min_salphc} */
  recompute_kon_rates(rates, state, genes, konStates, 1);

  for (i=0; i < NGENES; i++) {
    /* transport rates */
    transport[i] = kRNA * (float) (state->mRNA_nuclear_num[i]);
    rates->transport += transport[i];
    rates->transport_operations++;

    /* regenerate decay rates */
    mRNAdecay[i] = genes->mRNAdecay[i] * ((float) state->mRNA_cyto_num[i] + (float) state->mRNA_transl_cyto_num[i]);
    rates->mRNAdecay += mRNAdecay[i];
    rates->mRNAdecay_operations++;

  }

  /* recompute koffvalues for all sites */
  recompute_koff_rates(rates, state, genes, *koffvalues, TIME_INFINITY);

  /* recompute and cache the total rate in data structure */
  rates->total += rates->transport;
  rates->total += rates->mRNAdecay;
  rates->total += rates->pic_disassembly;

  /* 
   * convert the counts back into rates using the constants 
   */
  for (j=0; j < MAX_COPIES; j++) {
    rates->total += (float) rates->acetylation_num[j] * acetylate;
    rates->total += (float) rates->deacetylation_num[j] * deacetylate;
    rates->total += (float) rates->pic_assembly_num[j] * PICassembly;
    rates->total += (float) rates->transcript_init_num[j] * transcriptinit;    
  } 
}


int move_gene_copy(int from_copy,
                   int to_copy,
                   int gene,
                   Genotype *from_genotype,
                   Genotype *to_genotype,
                   CellState *from_state,
                   CellState *to_state,
                   int lastpos)
{
  int k;
  int site_id, new_site_id, gene_id, gene_copy;
  int start = from_genotype->site_id_pos[gene][from_copy][0];
  int end = from_genotype->site_id_pos[gene][from_copy][1];

  /* shift all sites in tf_bound_indexes */
  for (k = 0; k < from_state->tf_bound_num; k++) {
    site_id = from_state->tf_bound_indexes[k];
    gene_id = from_genotype->all_binding_sites[site_id].cisreg_id;
    gene_copy = from_genotype->all_binding_sites[site_id].gene_copy;

    /* check to see if TF is within the gene copy to be moved */
    if (site_id >= start && site_id <= end && gene_copy == from_copy && gene_id == gene) {
      /* fix  */
      if (start  == lastpos + 1)  /* if no offset required: keep site_id */
        new_site_id = site_id;
      else                        /* otherwise offset them by the difference */
        new_site_id = site_id  - (start - lastpos);
      LOG_VERBOSE_NOCELLID("gene=%d [copy=%d] orig site_id=%d, move to new [copy=%d] site_id=%d (new tf_bound_num=%d), lastpos=%d, offset=%d\n", 
                           gene, from_copy, site_id, to_copy, new_site_id, to_state->tf_bound_num, lastpos, (start - lastpos));
      to_state->tf_bound_indexes[to_state->tf_bound_num] = new_site_id;
      to_state->tf_bound_num++;
    }
  }

  /* shift TF hindered indexes */
  for (k = 0; k < from_state->tf_hindered_num; k++) {
    site_id = from_state->tf_hindered_indexes[k][0];
    gene_id = from_genotype->all_binding_sites[site_id].cisreg_id;
    gene_copy = from_genotype->all_binding_sites[site_id].gene_copy;
    if (site_id >= start && site_id <= end && gene_copy == from_copy && gene_id == gene) {
      new_site_id = site_id - (start - lastpos);
      to_state->tf_hindered_indexes[to_state->tf_hindered_num][0] = new_site_id;
      LOG_VERBOSE_NOCELLID("k=%d, gene=%d [copy=%d] orig hindering site_id=%d, move to new [copy=%d] site_id=%d (new tf_hindered_num=%d), lastpos= %d, offset=%d\n", 
                         k, gene, from_copy, site_id, to_copy, new_site_id, to_state->tf_hindered_num, lastpos, (start - lastpos));
      site_id = from_state->tf_hindered_indexes[k][1];
      new_site_id = site_id - (start - lastpos);
      to_state->tf_hindered_indexes[to_state->tf_hindered_num][1] = new_site_id;
      LOG_VERBOSE_NOCELLID("k=%d, gene=%d [copy=%d] orig hindered site_id=%d, move to new [copy=%d] site_id=%d (new tf_hindered_num=%d), lastpos= %d, offset=%d\n", 
                         k, gene, from_copy, site_id, to_copy, new_site_id, to_state->tf_hindered_num, lastpos, (start - lastpos));
      to_state->tf_hindered_num++;
    }
  }

  lastpos += (end - start) + 1;   /* update the last position in new binding site array */
  LOG_VERBOSE_NOCELLID("gene=%d [copy=%d] start=%d, end=%d, lastpos=%d\n", gene, from_copy, start, end, lastpos);
  return lastpos;
}

// TODO: move to helper lib?
void clone_queue(FixedEvent **start_orig,
                 FixedEvent **last_orig,
                 FixedEvent **start_clone,
                 FixedEvent **last_clone) 
{
  while (*start_orig != NULL) {
    LOG_NOCELLID("adding gene_id=%d time=%g to new clone, removing from orig\n", (*start_orig)->gene_id, (*start_orig)->time);
    //add_fixed_event((*start_orig)->gene_id, (*start_orig)->time, start_clone, last_clone);
    add_fixed_event((*start_orig)->gene_id, (*start_orig)->copy, (*start_orig)->time, start_clone, last_clone);
    delete_fixed_event_start(start_orig, last_orig);
  }

  if (*start_clone != NULL) {
    LOG_NOCELLID("clone queue is not empty: head is: gene_id=%d time=%g\n", (*start_clone)->gene_id, (*start_clone)->time);
  }
}

void split_mRNA(FixedEvent **start_clone,
                FixedEvent **last_clone,
                int count_clone[NGENES],
                FixedEvent **start_daughter,
                FixedEvent **last_daughter,
                int count_daughter[NGENES],
                FixedEvent **start_mother,
                FixedEvent **last_mother,
                int count_mother[NGENES],
                int i,
                float fraction) 
{
  // regenerate initial queue(s)
  int j;
  for (j=0; j < count_clone[i]; j++) {
    if (*start_clone != NULL) {
      int gene_id = (*start_clone)->gene_id;
      int copy = (*start_clone)->copy;
      float time = (*start_clone)->time;
      
      if (ran1(&seed) <= fraction) { /* fraction of time move to daughter */
        count_daughter[gene_id]++;
        add_fixed_event(gene_id, copy, time, start_daughter, last_daughter);
        LOG_NOCELLID("move event at time=%g on gene %d to daughter=%d of total count=%d\n", 
                     time, gene_id, count_daughter[i], count_clone[i]);
      } else {  /* (1-fraction of time move to mother */
        count_mother[gene_id]++;
        add_fixed_event(gene_id, copy, time, start_mother, last_mother);
        LOG_NOCELLID("move event at time=%g on gene %d to mother=%d of total count=%d\n", 
                     time, gene_id, count_daughter[i], count_clone[i]);
      }
      /* remove from original queue */
      LOG_NOCELLID("removing event at time=%g on gene %d from clone of queue\n", time, gene_id);
      delete_fixed_event_start(start_clone, last_clone);
    }
  }
}


void clone_cell(Genotype *genes_orig,
                CellState *state_orig,
                GillespieRates *rates_orig,
                Genotype *genes_clone,
                CellState *state_clone,
                GillespieRates *rates_clone)
{
  int i, k, p;

  state_clone->founder_id = state_orig->founder_id;
  state_clone->in_s_phase = state_orig->in_s_phase;
  state_clone->burn_in = state_orig->burn_in;
  state_clone->divisions = state_orig->divisions;

  state_clone->cell_size = state_orig->cell_size;
  state_clone->growth_rate =  state_orig->growth_rate;

  state_clone->mRNA_transcr_time_end = NULL;
  state_clone->mRNA_transcr_time_end_last = NULL;
  state_clone->mRNA_transl_time_end = NULL;
  state_clone->mRNA_transl_time_end_last = NULL;
  
  state_clone->RTlnKr = state_orig->RTlnKr;
  state_clone->temperature = state_orig->temperature;

  /* copy Genotype */
  // TODO: check to see how much of these states should be set for all 
  // genotypes or copied from all genotypes, either way it should be 
  // copied from the parent cell

  /* clone the cis-regulatory sequence */
  memcpy(genes_clone->cisreg_seq, genes_orig->cisreg_seq, sizeof(char [NGENES][MAX_COPIES][CISREG_LEN]));

  /* clone all other aspects of genotype */
  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      genes_clone->site_id_pos[i][p][0] =  genes_orig->site_id_pos[i][p][0];
      genes_clone->site_id_pos[i][p][1] =  genes_orig->site_id_pos[i][p][1];
      genes_clone->activating[i][p]=  genes_orig->activating[i][p];
      genes_clone->pic_disassembly[i][p]=  genes_orig->pic_disassembly[i][p];
    }
    genes_clone->sites_per_gene[i] = genes_orig->sites_per_gene[i];
    genes_clone->copies[i] = genes_orig->copies[i];
    genes_clone->replication_time[i] =  genes_orig->replication_time[i];
    genes_clone->mRNAdecay[i] =  genes_orig->mRNAdecay[i];
    genes_clone->proteindecay[i] =  genes_orig->proteindecay[i];
    genes_clone->translation[i] =  genes_orig->translation[i];
  }

  /* clone the TF sequence */
  memcpy(genes_clone->tf_seq, genes_orig->tf_seq, 
         sizeof(char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN]));

  /* clone the hindrance positions */
  for (i=0; i < TFGENES; i++) {
    genes_clone->hindrance_positions[i] = genes_orig->hindrance_positions[i];
  }

  // make a pointer to original all_binding_sites, we don't modify it yet
  genes_clone->all_binding_sites =  genes_orig->all_binding_sites;

  state_clone->tf_bound_num  = 0;  //state_orig->tf_bound_num;
  state_clone->tf_hindered_num = 0; // state_orig->tf_hindered_num;

  state_clone->tf_bound_indexes = NULL;
  state_clone->tf_hindered_indexes = NULL;

  // copy bound indexes if there are any
  if (state_orig->tf_bound_num > 0) {
    state_clone->tf_bound_indexes = malloc(state_orig->tf_bound_num*sizeof(int));
    for (k=0; k < state_orig->tf_bound_num; k++) {
      state_clone->tf_bound_indexes[k] = state_orig->tf_bound_indexes[k];
      state_clone->tf_bound_num++;
    }
  }

  // copy bound hindrance positions if there are any
  if (state_orig->tf_hindered_num > 0) {
    state_clone->tf_hindered_indexes = malloc(2*state_orig->tf_hindered_num*sizeof(int));
    for (k=0; k < state_orig->tf_hindered_num; k++) {
      state_clone->tf_hindered_indexes[k][0] = state_orig->tf_hindered_indexes[k][0];
      state_clone->tf_hindered_indexes[k][1] = state_orig->tf_hindered_indexes[k][1];
      state_clone->tf_hindered_num++;
    }
  }

  /* copy activation state */
  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      state_clone->active[i][p] = state_orig->active[i][p];
      state_clone->state_change_ids[ACETYLATION][p][i] = state_orig->state_change_ids[ACETYLATION][p][i];
      state_clone->state_change_ids[DEACETYLATION][p][i] = state_orig->state_change_ids[DEACETYLATION][p][i];
      state_clone->state_change_ids[PICASSEMBLY][p][i] = state_orig->state_change_ids[PICASSEMBLY][p][i];
      state_clone->state_change_ids[TRANSCRIPTINIT][p][i] = state_orig->state_change_ids[TRANSCRIPTINIT][p][i];
      state_clone->state_change_ids[PICDISASSEMBLY][p][i] = state_orig->state_change_ids[PICDISASSEMBLY][p][i];

      rates_clone->acetylation_num[p] = rates_orig->acetylation_num[p];
      rates_clone->deacetylation_num[p] =  rates_orig->deacetylation_num[p];
      rates_clone->pic_assembly_num[p] = rates_orig->pic_assembly_num[p];
      rates_clone->transcript_init_num[p] = rates_orig->transcript_init_num[p];
      rates_clone->pic_disassembly_num[p] = rates_orig->pic_disassembly_num[p];
      state_clone->mRNA_transcr_num[i][p] = state_orig->mRNA_transcr_num[i][p];
    }
    state_clone->protein_conc[i] = state_orig->protein_conc[i];;
    state_clone->mRNA_cyto_num[i] = state_orig->mRNA_cyto_num[i];
    state_clone->mRNA_nuclear_num[i] = state_orig->mRNA_nuclear_num[i];
    state_clone->mRNA_transl_cyto_num[i] = state_orig->mRNA_transl_cyto_num[i];
  }

  /* clone queue and empty original */
  clone_queue(&(state_orig->mRNA_transcr_time_end), &(state_orig->mRNA_transcr_time_end_last),
              &(state_clone->mRNA_transcr_time_end), &(state_clone->mRNA_transcr_time_end_last));
  clone_queue(&(state_orig->mRNA_transl_time_end), &(state_orig->mRNA_transl_time_end_last),
              &(state_clone->mRNA_transl_time_end), &(state_clone->mRNA_transl_time_end_last));
}

void initialize_new_cell_genotype(Genotype *genes, Genotype *genes_clone)
{
  int i, p;

  /* initialize hindrance for all TFGENES */
  for (p=0; p < TFGENES; p++) {
    genes->hindrance_positions[p] = genes_clone->hindrance_positions[p];
  }

  /* initialize the non-cisregulatory parts of the genotype */
  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      genes->activating[i][p]=  genes_clone->activating[i][p];
      genes->pic_disassembly[i][p]=  genes_clone->pic_disassembly[i][p];
    }
    genes->replication_time[i] =  genes_clone->replication_time[i];
    genes->mRNAdecay[i] =  genes_clone->mRNAdecay[i];
    genes->proteindecay[i] =  genes_clone->proteindecay[i];
    genes->translation[i] =  genes_clone->translation[i];
  }
}

void initialize_new_cell_state(CellState *state, CellState state_clone, 
                               GillespieRates *rates, double fraction)
{
  int i, j;
  // reset pointers
  state->mRNA_transcr_time_end = NULL;
  state->mRNA_transcr_time_end_last = NULL;
  state->mRNA_transl_time_end = NULL;
  state->mRNA_transl_time_end_last = NULL;
  state->replication_time_end = NULL;
  state->replication_time_end_last = NULL;

  fflush(stdout);

  if (state->mRNA_transcr_time_end != NULL) {
    printf("after delete_queues  mRNA_transcr_time_end is not NULL\n");
    printf("time=%g, gene_id=%d\n", state->mRNA_transcr_time_end->time, state->mRNA_transcr_time_end->gene_id);
  }

  if (state->replication_time_end != NULL) {
    printf("after delete_queues  replication_time_end is not NULL\n");
    printf("time=%g, gene_id=%d\n", state->replication_time_end->time, state->replication_time_end->gene_id);
  }

  /* reset the length of TF data structures to zero, we have to reconstruct them */
  state->tf_bound_num = 0;
  state->tf_hindered_num = 0;

  state->in_s_phase = 0;   /* reset S phase state to 0 */
  // TODO: check to see whether we want burn-in at the birth of each new cell ?
  state->burn_in = 0;      /* only do burn-in at beginning of runs */
  state->division_time = TIME_INFINITY;  /* reset division time for this cell */
  state->divisions = state_clone.divisions; /* copy division counter */

  state->cell_size = state_clone.cell_size*fraction;   /* reset cell size */

  /* TODO: check! growth rate, take instantaneous growth rate just before
     division, this will be updated after first new time step */
  // TODO, check that this is always inherited always from mother cell
  state->growth_rate = state_clone.growth_rate;

  /* keep thermodynamic state the same */ 
  state->RTlnKr = state_clone.RTlnKr;
  state->temperature = state_clone.temperature;

  /* initialize the rate counts */
  for (i=0; i < MAX_COPIES; i++) {
    for (j=0; j < NGENES; j++)
      state->state_change_ids[ACETYLATION][i][j] = -1;
    rates->acetylation_num[i]=0;
    rates->deacetylation_num[i]=0;
    rates->pic_assembly_num[i]=0;
    rates->transcript_init_num[i]=0;
    rates->pic_disassembly_num[i]=0;
  }
}

void initialize_new_cell_state_change_ids(CellState *state,
                                        CellState state_clone,
                                        int type,
                                        int count[MAX_COPIES],
                                        int clone_count[MAX_COPIES],
                                        int copy1,
                                        int copy2,
                                        int gene_id) 
{
  int j;
  for (j=0; j < (clone_count[copy1]); j++) {
    int gene = state_clone.state_change_ids[type][copy1][j];
    if (gene == gene_id) {
      state->state_change_ids[type][0][count[0]] = gene;
      LOG_VERBOSE(" state_change_ids[%d][0][%2d]=%2d (copy %2d)\n", type, j, state->state_change_ids[type][0][count[0]], copy1);
      (count[0])++;
    }
  }

  for (j=0; j < clone_count[copy2]; j++) {
    if (state_clone.state_change_ids[type][copy2][j] == gene_id) {
      state->state_change_ids[type][1][count[1]] = state_clone.state_change_ids[type][copy2][j];
      LOG_VERBOSE(" state_change_ids[%d][1][%2d]=%2d (copy %2d)\n", type, j, state->state_change_ids[type][1][count[1]], copy2);
      (count[1])++;
    }
  }

  /* reset the 3rd and 4th copies to zero as cell will be diploid again immediately after division */
  count[2] = 0;
  count[3] = 0;
}

void initialize_new_cell_gene(Genotype *genes, Genotype genes_clone, 
                              CellState *state, CellState *state_clone,
                              GillespieRates *rates, GillespieRates *rates_clone,
                              int gene_id, int copy1, int copy2, int *lastpos)
{
  int j, k;

  for (k=0; k < CISREG_LEN; k++) {
    genes->cisreg_seq[gene_id][0][k] = genes_clone.cisreg_seq[gene_id][copy1][k];
    genes->cisreg_seq[gene_id][1][k] = genes_clone.cisreg_seq[gene_id][copy2][k];
    genes->cisreg_seq[gene_id][2][k] = genes_clone.cisreg_seq[gene_id][copy1][k];
    genes->cisreg_seq[gene_id][3][k] = genes_clone.cisreg_seq[gene_id][copy2][k];
  }

  /* don't update TF if this gene_id is not controlling a TF */
  if (gene_id < TFGENES) {
    for (k=0; k < TF_ELEMENT_LEN; k++) {
      genes->tf_seq[gene_id][0][k] = genes_clone.tf_seq[gene_id][copy1][k];
      genes->tf_seq[gene_id][1][k] = genes_clone.tf_seq[gene_id][copy2][k];
      genes->tf_seq[gene_id][2][k] = genes_clone.tf_seq[gene_id][copy1][k];
      genes->tf_seq[gene_id][3][k] = genes_clone.tf_seq[gene_id][copy2][k];
    }
  }


  if (genes_clone.copies[gene_id] - 1 >= copy1) {
    *lastpos = move_gene_copy(copy1, 0, gene_id, &genes_clone, genes, state_clone, state, *lastpos);
  }

  if (genes_clone.copies[gene_id] - 1 >= copy2) {
    *lastpos = move_gene_copy(copy2, 1, gene_id, &genes_clone, genes, state_clone, state, *lastpos);
  }

  fflush(stdout);

  /* copy activation state */
  state->active[gene_id][0] = state_clone->active[gene_id][copy1];
  state->active[gene_id][1] = state_clone->active[gene_id][copy2];

  state->active[gene_id][2] = ON_WITH_NUCLEOSOME;
  state->active[gene_id][3] = ON_WITH_NUCLEOSOME;

  /* copy statechangeID state */
  // TODO: ultimately maybe acteylationCounts etc. should be moved to "CellState"
  // since they are not cached values

  /* ACETYLATION */
  initialize_new_cell_state_change_ids(state, *state_clone, ACETYLATION,
                                     rates->acetylation_num, rates_clone->acetylation_num,
                                     copy1, copy2, gene_id);

  /* DEACETYLATION */
  initialize_new_cell_state_change_ids(state, *state_clone, DEACETYLATION,
                                     rates->deacetylation_num, rates_clone->deacetylation_num,
                                     copy1, copy2, gene_id);
  /* PICASSEMBLY */
  initialize_new_cell_state_change_ids(state, *state_clone, PICASSEMBLY,
                                     rates->pic_assembly_num, rates_clone->pic_assembly_num,
                                     copy1, copy2, gene_id);
  /* TRANSCRIPTINIT */
  initialize_new_cell_state_change_ids(state, *state_clone, TRANSCRIPTINIT,
                                     rates->transcript_init_num, rates_clone->transcript_init_num,
                                     copy1, copy2, gene_id);
  /* PICDISASSEMBLY */
  initialize_new_cell_state_change_ids(state, *state_clone, PICDISASSEMBLY,
                                     rates->pic_disassembly_num, rates_clone->pic_disassembly_num,
                                     copy1, copy2, gene_id);

  /* initialize the per-gene (and/or copy) values of the mRNA counts
     these will be updated later */
  for (j=0; j < MAX_COPIES; j++)
    state->mRNA_transcr_num[gene_id][j] = 0;
  state->mRNA_cyto_num[gene_id] = 0;
  state->mRNA_nuclear_num[gene_id] = 0;
  state->mRNA_transl_cyto_num[gene_id] = 0;

  /* for transcribing mRNAs: split up the genetic material, as mRNA is
     actually physically attached to the DNA during the replication and
     doesn't go randomly to one or other of the cells split up
     mRNA_transcr_num, along with FixedTime events */

  FixedEvent *timeEnd = state_clone->mRNA_transcr_time_end;

  while (timeEnd != NULL) {
    for (j=0; j < MAX_COPIES; j++)
      LOG_VERBOSE("mRNA_transcr_num[%2d][%d]=%d\n", 
                  gene_id, j, state_clone->mRNA_transcr_num[gene_id][j]);
    
    int gene_idQueue = timeEnd->gene_id;
    int copy = timeEnd->copy;
    float time = timeEnd->time;
    
    if (gene_idQueue == gene_id && copy == copy1) { 
      state->mRNA_transcr_num[gene_id][0]++;
      add_fixed_event(gene_id, 0, time, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end));
      printf("move mRNATranscrTime event at time=%g on gene %d copy=%d to copy=0 of total mRNA_transcr_num=%d\n", 
             time, gene_id, copy1, state_clone->mRNA_transcr_num[gene_id][copy1]);
      /* remove from original queue */
      printf("removing mRNATranscrTime event at time=%g on gene %2d (copy %d) from clone of queue\n", time, gene_id, copy1);
      delete_fixed_event(gene_id, copy1, 0, &(state_clone->mRNA_transcr_time_end), &(state_clone->mRNA_transcr_time_end_last));
      (state_clone->mRNA_transcr_num[gene_id][copy1])--;
    } else if (gene_idQueue == gene_id && copy == copy2) {  
      state->mRNA_transcr_num[gene_id][1]++;
      add_fixed_event(gene_id, 1, time, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end));
      printf("move mRNATranscrTime event at time=%g on gene %d copy=%d to copy=1 of total mRNA_transcr_num=%d\n", 
             time, gene_id, copy2, state_clone->mRNA_transcr_num[gene_id][copy2]);
      /* remove from original queue */
      printf("removing mRNATranscrTime event at time=%g on gene %2d (copy %d) from clone of queue\n", time, gene_id, copy2);
      delete_fixed_event(gene_id, copy2, 0, &(state_clone->mRNA_transcr_time_end), &(state_clone->mRNA_transcr_time_end_last));
      (state_clone->mRNA_transcr_num[gene_id][copy2])--;
    }
    timeEnd = timeEnd->next;
    fflush(stdout);
  }
}

void realloc_cell_memory(CellState *state, float **koffvalues) 
{
  // TODO: check to see if we go back to exact computation of memory size needed
  //state->tf_bound_indexes = realloc(state->tf_bound_indexes, 2*(state->tf_bound_num+1)*sizeof(int));
  state->tf_bound_indexes = realloc(state->tf_bound_indexes, 10*maxbound*sizeof(int));

  //state->tf_hindered_indexes = realloc(state->tf_hindered_indexes, 10*state->tf_hindered_num*sizeof(int));
  state->tf_hindered_indexes = realloc(state->tf_hindered_indexes, 100*maxbound*sizeof(int));
  
  // reallocate koffvalues memory
  //*koffvalues = realloc(*koffvalues, 2*(state->tf_bound_num+1)* sizeof(float)); 
  *koffvalues = realloc(*koffvalues, 100*maxbound*sizeof(float)); 

  if (!state->tf_bound_indexes || !state->tf_hindered_indexes) {
    LOG_ERROR_NOCELLID("memory allocation error cell\n");
    exit(1);
  }
}

void do_cell_division(int motherID,
                      int daughterID,
                      Genotype *mother,
                      CellState *mother_state,
                      GillespieRates *mother_rates,
                      KonStates *mother_konStates,
                      float **mother_koffvalues,
                      float mother_transport[NGENES],
                      float mother_mRNAdecay[NGENES],
                      
                      Genotype *daughter,
                      CellState *daughter_state,
                      GillespieRates *daughter_rates,
                      KonStates *daughter_konStates,
                      float **daughter_koffvalues,
                      float daughter_mRNAdecay[NGENES],
                      float daughter_transport[NGENES],
                      float fraction,
                      float x,
                      float dt)
{
  int i, j, total;

  /* clone of cell */
  Genotype genes_clone;
  CellState state_clone;
  GillespieRates rates_clone;

  int daughter_copy1;  /* 1 */
  int daughter_copy2;  /* 2 */
  int mother_copy1;    /* 0 */
  int mother_copy2;    /* 3 */
  double r;             /* random number for gene assortment */
  int lastpos_daughter = 0;
  int lastpos_mother = 0;
  int original_bind_count = mother->binding_sites_num;
  int no_replace_mother = (motherID != daughterID) ? 1: 0;  /* set if daughter cell replaces mother 
                                                               in that case, we discard all updating of 
                                                               the mother cell */
  /* clone mother cell */
  clone_cell(mother, mother_state, mother_rates,
             &genes_clone, &state_clone, &rates_clone);
  
  /* free the existing memory for the time queue */
  delete_queues(daughter_state);
  if (no_replace_mother)
    delete_queues(mother_state);

  /* initialize the non-cis-regulatory part of the genotype 
     of the new daughter cell based on clone of mother */
  initialize_new_cell_genotype(daughter, &genes_clone);
  // print_genotype(daughter, daughterID);

  /* initialize the state of the new cell */
  /* also reset size to appropriate scale of cell size */
  initialize_new_cell_state(daughter_state, state_clone, daughter_rates, fraction);
  // TODO: fix to use scale according to the current size as per below
  //initialize_new_cell_state(daughter_state, state_clone, daughter_rates, fraction*state_clone.cell_size);

  /* set the founder_id of the daughter cell from the mother cell */
  daughter_state->founder_id = state_clone.founder_id;
  daughter_state->divisions = 0;    /* daughter cell resets divisions */ 

  printf("daughter cell %03d is founded by cell %03d with %2d divisions\n", 
         daughterID, daughter_state->founder_id, daughter_state->divisions);
  LOG_NOCELLID("daughter cell %03d is founded by cell %03d with %2d divisions\n", 
               daughterID, daughter_state->founder_id, daughter_state->divisions);

  if (no_replace_mother) {
    initialize_new_cell_genotype(mother, &genes_clone);
    //print_genotype(mother, motherID);

    initialize_new_cell_state(mother_state, state_clone, mother_rates, (1-fraction));
    // TODO (see above)
    //initialize_new_cell_state(mother_state, state_clone, mother_rates, (1-fraction)*state_clone.cell_size);
    mother_state->divisions++;      /* update divisions in mother */
    printf("mother   cell %03d is founded by cell %03d with %2d divisions\n", 
         motherID, mother_state->founder_id, mother_state->divisions);
    LOG_NOCELLID("mother   cell %03d is founded by cell %03d with %2d divisions\n", 
                 motherID, mother_state->founder_id, mother_state->divisions);
  }

  LOG_VERBOSE_NOCELLID("[cell %03d] (mother) total tf_bound_num=%d, tf_hindered_num=%d\n", motherID, 
                       state_clone.tf_bound_num, state_clone.tf_hindered_num);

  //print_all_binding_sites(mother->copies, mother->all_binding_sites, mother->binding_sites_num, 
  //                        mother->tf_seq, mother->cisreg_seq, mother->site_id_pos); 

  if (verbose) {
    LOG_NOCELLID("[cell %03d]: acetylation counts:\n", motherID);
    total = 0;
    for (j=0; j < MAX_COPIES; j++)  {
      LOG_NOCELLID(" before mother acetylation_num[%2d]=%2d\n", j, rates_clone.acetylation_num[j]);
      total += rates_clone.acetylation_num[j];
    }
    LOG_NOCELLID(" before mother total acetylation=%d\n", total);
  }

  /* now split up the genes into one or other of the new cells */
  for (i=0; i < NGENES; i++) {

    /* reset the number of copies of gene in mother and daughter after
       division to original ploidy */
    daughter->copies[i] = current_ploidy;
    if (no_replace_mother)
      mother->copies[i] = current_ploidy;

    // TODO: assume diploid for the moment, generalize for haploid
    
    r = ran1(&seed);  /* random number for independent assortment of
                         each gene */
    if (r<0.25) {
      /* take copy 0 and 1 from mother to form copy 0 and 1 in daughter */
      daughter_copy1 = 0;  /* mother    daughter */
      daughter_copy2 = 1;  /* copy 0 -> copy 0  
                              copy 1 -> copy 1  
                              copy 0 -> copy 2
                              copy 1 -> copy 3 */
      /* while mother keeps copy 2 and 3 */
      mother_copy1 = 2;    /* mother    mother */
      mother_copy2 = 3;    /* copy 2 -> copy 0
                              copy 3 -> copy 1
                              copy 2 -> copy 2
                              copy 3 -> copy 3 */
    } else if (r<0.5) {
      daughter_copy1 = 1;
      daughter_copy2 = 2;
      mother_copy1 = 0;
      mother_copy2 = 3;
    } else if (r<0.75) {
      daughter_copy1 = 2;
      daughter_copy2 = 3;
      mother_copy1 = 0;
      mother_copy2 = 1;
    } else { 
      daughter_copy1 = 0;
      daughter_copy2 = 3;
      mother_copy1 = 1;
      mother_copy2 = 2;
    }

    LOG_NOCELLID("initialize daughter=%2d, gene_id=%2d\n", daughterID, i);
    initialize_new_cell_gene(daughter, genes_clone, 
                             daughter_state, &state_clone,
                             daughter_rates, &rates_clone,
                             i, daughter_copy1, daughter_copy2, &lastpos_daughter);

    if (no_replace_mother) {
      LOG_NOCELLID("initialize mother  =%2d, gene_id=%2d\n", motherID, i);
      initialize_new_cell_gene(mother, genes_clone, 
                               mother_state, &state_clone,
                               mother_rates, &rates_clone,
                               i, mother_copy1, mother_copy2, &lastpos_mother);
    }
  }

  if (verbose) {
    total = 0;
    for (j=0; j < MAX_COPIES; j++)  {
      LOG_NOCELLID(" after daughter acetylation_num[%2d]=%2d\n", j, daughter_rates->acetylation_num[j]);
      LOG_NOCELLID(" after mother acetylation_num[%2d]=%2d\n", j, mother_rates->acetylation_num[j]);
      total += daughter_rates->acetylation_num[j];
      total += mother_rates->acetylation_num[j];
    }
    LOG_NOCELLID(" after total acetylation=%d\n", total);

    // check that we have emptied the list of transcribing mRNAs
    total = 0;
    for (i=0; i < NGENES; i++) 
      for (j=0; j < MAX_COPIES; j++) 
        total += state_clone.mRNA_transcr_num[i][j];
    
    LOG_NOCELLID("mRNA_transcr_num=%2d left in state_clone after moving to mother+daughter\n", total);
  }

  realloc_cell_memory(daughter_state, daughter_koffvalues);
  if (no_replace_mother) {
    realloc_cell_memory(mother_state, mother_koffvalues);
  }

  /* recompute *all* binding sites in daughter, then relabel sites */
  calc_all_binding_sites(daughter->copies, 
                         daughter->cisreg_seq, 
                         daughter->tf_seq, 
                         &(daughter->binding_sites_num),
                         &(daughter->all_binding_sites),
                         //daughter->hindrance_positions,
                         mother->hindrance_positions,
                         daughter->sites_per_gene,
                         daughter->site_id_pos); 

  //print_all_binding_sites(daughter->copies, daughter->all_binding_sites, daughter->binding_sites_num, 
  //                        daughter->tf_seq, daughter->cisreg_seq, daughter->site_id_pos); 

  if (no_replace_mother) {
    /* recompute *all* binding sites in mother, then relabel sites */
    calc_all_binding_sites(mother->copies, 
                           mother->cisreg_seq, 
                           mother->tf_seq, 
                           &(mother->binding_sites_num),
                           &(mother->all_binding_sites),
                           mother->hindrance_positions,
                           mother->sites_per_gene,
                           mother->site_id_pos); 
  }

  //print_all_binding_sites(mother->copies, mother->all_binding_sites, mother->binding_sites_num, 
  //                        mother->tf_seq, mother->cisreg_seq, mother->site_id_pos); 

  if (no_replace_mother) {
    LOG_NOCELLID("original number of binding sites=%d should = (mother=%d + daughter=%d) = %d\n", 
                 original_bind_count, mother->binding_sites_num, daughter->binding_sites_num, mother->binding_sites_num + daughter->binding_sites_num);
    if (original_bind_count != mother->binding_sites_num + daughter->binding_sites_num) {
      LOG_ERROR_NOCELLID("original number of binding sites=%d  != (mother=%d + daughter=%d) = %d\n", 
                         original_bind_count, mother->binding_sites_num, daughter->binding_sites_num, mother->binding_sites_num + daughter->binding_sites_num);
      exit(0);
    }
  }
  
  /* split up the volume of the cell */
  for (i=0; i < NPROTEINS; i++) {

    // first protein
    daughter_state->protein_conc[i] = fraction * state_clone.protein_conc[i];
    if (no_replace_mother) {
      mother_state->protein_conc[i] = (1-fraction) * state_clone.protein_conc[i];
      LOG_VERBOSE_NOCELLID("daughter=%g (%g), mother=%g (%g) = total protein=%g\n", 
                   daughter_state->protein_conc[i], fraction, mother_state->protein_conc[i], (1-fraction), state_clone.protein_conc[i]);
    }
  }

  /* now loop over all transcribing genes */
  for (i=0; i < NGENES; i++) {
    // mRNAs in cytoplasm (not translating)
    daughter_state->mRNA_cyto_num[i] = rint(fraction * state_clone.mRNA_cyto_num[i]);
    if (no_replace_mother) {
      mother_state->mRNA_cyto_num[i] =  state_clone.mRNA_cyto_num[i] - daughter_state->mRNA_cyto_num[i];
      LOG_VERBOSE_NOCELLID("daughter=%d, mother=%d of total mRNA_cyto_num=%d\n", 
                   daughter_state->mRNA_cyto_num[i], mother_state->mRNA_cyto_num[i], state_clone.mRNA_cyto_num[i]);
    }

    // mRNAs in nucleus
    daughter_state->mRNA_nuclear_num[i] = rint(fraction * mother_state->mRNA_nuclear_num[i]);
    if (no_replace_mother) {
      mother_state->mRNA_nuclear_num[i] =  state_clone.mRNA_nuclear_num[i] - daughter_state->mRNA_nuclear_num[i];
      LOG_VERBOSE_NOCELLID("daughter=%d, mother=%d of total mRNA_nuclear_num=%d\n", 
                   daughter_state->mRNA_nuclear_num[i], mother_state->mRNA_nuclear_num[i], state_clone.mRNA_nuclear_num[i]);
    }

    // split up mRNA_transl_cyto_num, along with FixedTime events
    split_mRNA(&(state_clone.mRNA_transl_time_end), &(state_clone.mRNA_transl_time_end_last), state_clone.mRNA_transl_cyto_num,
               &(daughter_state->mRNA_transl_time_end), &(daughter_state->mRNA_transl_time_end_last), daughter_state->mRNA_transl_cyto_num,
               &(mother_state->mRNA_transl_time_end), &(mother_state->mRNA_transl_time_end_last), mother_state->mRNA_transl_cyto_num,
               i, fraction);
  }  
  fflush(stdout);

  /* recompute rates in daughter */
  recalibrate_cell(daughter_rates, daughter_state, daughter,
                   daughter_konStates, daughter_koffvalues,
                   daughter_mRNAdecay, daughter_transport, dt);

  if (no_replace_mother) {
    /* recompute rates in mother */
    recalibrate_cell(mother_rates, mother_state, mother,
                     mother_konStates, mother_koffvalues,
                     mother_mRNAdecay, mother_transport, dt);
    LOG_VERBOSE_NOCELLID("tf_bound_num=%d (motherID=%d), tf_bound_num=%d (daughterID=%d)\n", 
                         mother_state->tf_bound_num, motherID, daughter_state->tf_bound_num, daughterID);
    LOG_VERBOSE_NOCELLID("tf_hindered_num=%d (motherID=%d), tf_hindered_num=%d (daughterID=%d)\n", 
                         mother_state->tf_hindered_num, motherID, daughter_state->tf_hindered_num, daughterID);
  }

  /* free the memory associated with temporary copy of mother */
  free_mem_CellState(&state_clone);
}

void log_snapshot(GillespieRates *rates,
                  CellState *state,
                  Genotype *genes,
                  KonStates *konStates,
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
             rates->pic_disassembly, rates->salphc+(konrate), konStates->nkon, (rates->salphc+(konrate))/(float)konStates->nkon);
  
  for (p=0; p < MAX_COPIES; p++) {
    LOG_NOFUNC(" acetylation=%g (copy %d)\n deacetylation=%g (copy %d)\n PIC assembly=%g (copy %d)\n transcriptinit=%g (copy %d)\n",
               (float)rates->acetylation_num[p]*acetylate, p, (float)rates->deacetylation_num[p]*deacetylate, p, 
               (float)rates->pic_assembly_num[p]*PICassembly, p, (float)rates->transcript_init_num[p]*transcriptinit, p);
  }
  LOG_NOFUNC(" total rates=%g=%g+%g\n", rates->total + (konrate), rates->total, konrate);
  LOG_NOFUNC(" total free=%d + total bound=%d + total hindered=%d = total sites=%d\n", 
             konStates->nkon, state->tf_bound_num, state->tf_hindered_num, genes->binding_sites_num);
  for (i = 0; i < TFGENES; i++) {
    nkon += konStates->kon_list[i]->site_count;
    LOG_NOFUNC(" unoccupied binding sites=%d available for TF=%d \n", konStates->kon_list[i]->site_count, i);
  }
  LOG_NOFUNC(" nkon recomputed=%d\n", nkon);
  LOG_NOFUNC("\n");
}

/*
 * run the model for a specified cell for a single timestep
 */
int do_single_timestep(Genotype *genes, 
                       CellState *state, 
                       KonStates *konStates, 
                       GillespieRates *rates, 
                       float *t,
                       float *koffvalues,
                       float transport[NGENES],
                       float mRNAdecay[NGENES],
                       float *x,
                       float *dt,
                       float *konrate,
                       TimeCourse *timecoursestart[NPROTEINS],
                       TimeCourse *timecourselast[NPROTEINS],
                       int maxbound2,
                       int maxbound3,
                       int no_fixed_dev_time) 
{
  int i, j;

  int event;     /* boolean to keep track of whether FixedEvent has ended */
  int total;     /* total possible translation events */
  
  float konrate2, fixed_time;

  /* if asked, check for rounding error and recompute rates */

  /* check drift from koff every 1e6 operations */
  if (recompute_koff && (rates->koff_operations >= 1e6 || rates->koff < 0.0)) {
    LOG_WARNING("(t=%g) after %d operations recompute koff rates\n", *t, rates->koff_operations);
    recompute_koff_rates(rates, state, genes, koffvalues, *t);
  }

  /* check drift from various kon rates every 1e6 operations */
  if (recompute_kon && (rates->salphc_operations >= 1e6 || rates->salphc < 0.0)) {
    LOG_WARNING("(t=%g) after %d operations recompute kon rates\n", *t, rates->salphc_operations);
    recompute_kon_rates(rates, state, genes, konStates, 0);
  } 
  
  if (*t > 0.00005 && state->burn_in) {
    printf("recalibrating cell %3d after burn-in!\n", state->cell_id);
    LOG("recalibrating cell %3d after burn-in!\n", state->cell_id);
    log_snapshot(rates,
                 state,
                 genes,
                 konStates,
                 &koffvalues,
                 mRNAdecay,
                 transport, 
                 *konrate,
                 *x,
                 *t);

    kon = kon_after_burnin;

    recalibrate_cell(rates,
                     state,
                     genes,
                     konStates,
                     &koffvalues,
                     mRNAdecay,
                     transport,
                     *dt); 

    // TODO: check!
    /* compute konrate (which is constant over the Gillespie step) */
    if (konStates->nkon==0) {
      *konrate = (-rates->salphc);    /* all binding sites are occupied, total rate of binding is zero */
      LOG_WARNING("konStates->nkon=%d, konrate=%g\n", konStates->nkon, *konrate);
    } else  {
      calc_kon_rate(*dt, konStates, konrate);           /* otherwise compute konrate */
    }

    log_snapshot(rates,
                 state,
                 genes,
                 konStates,
                 &koffvalues,
                 mRNAdecay,
                 transport,
                 *konrate,
                 *x,
                 *t);
    state->burn_in = 0;
  } 

  /* compute S-phase offsets */
  if (critical_size > 0.0 && state->cell_size >= critical_size && !state->in_s_phase)  { /* run until checkpoint reached */
    reach_s_phase(state, genes, *t);
    /* current time plus 30 mins for S phase and a further 30 mins of growth after S phase */
    state->division_time = *t + time_s_phase + time_g2_phase;   
    LOG("at t=%g add division time=%g, cell_size=%g\n", *t, state->division_time, state->cell_size);
  }
  
  // TODO: currently disable printing out the TF occupancy
#if 0
  print_tf_occupancy(state, genes->all_binding_sites, *t);
  print_rounding(state, rates, *t);
#endif
  
  *x = expdev(&seed);        /* draw random number */

  /* check for rounding error */
  // TODO 2009-03-03: will remove very soon after more testing
  if (rates->koff < 0.0){
    konrate2 = 0.0;
    for (i=0; i < state->tf_bound_num; i++) konrate2 += koffvalues[i];
    if ((verbose) || konrate2>0.0) {
      LOG_WARNING("koffvalues add up to %g rates->koff=%g < 0\n", konrate2, rates->koff);
    }
    rates->koff = konrate2;
  }
  
  /* do first Gillespie step to chose next event */
  calc_dt(x, dt, rates, konStates, mRNAdecay, genes->mRNAdecay,
          state->mRNA_cyto_num, state->mRNA_transl_cyto_num, state->cell_id);

  if (*dt < 0.0) {
    LOG_ERROR("dt=%g is negative after calc_dt, t=%g\n", *dt, *t);
  }

  LOG_VERBOSE("next stochastic event due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);
  
  // TODO: remove?  not needed anymor?
  if (!(state->mRNA_transcr_time_end_last)) {
    for (i=0; i < NGENES; i++)
      for (j=0; j < MAX_COPIES; j++) {
        LOG_VERBOSE("%d transcription events on gene %2d (copy %2d)\n", state->mRNA_transcr_num[i][j], i, j);
      }
  }
  
  /* first check to see if a fixed event occurs in current t->dt window,
     or in tdevelopment if running for a fixed development time */
  fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);

  event = does_fixed_event_end(state->mRNA_transl_time_end,
                               state->mRNA_transcr_time_end,
                               state->replication_time_end,
                               fixed_time);
  
  /* while there are either transcription or translation events
     occuring in current t->dt window */
  while (event > 0) {
    *konrate = (*x)/(*dt);
    
    switch (event) {
    case 1:   /* if a transcription event ends */
      end_transcription(dt, *t, state, transport, rates);
      
      update_protein_conc_cell_size(state->protein_conc, state, genes, *dt,
                                    rates, konStates, *t,
                                    timecoursestart, timecourselast,
                                    state->protein_conc);
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
      update_protein_conc_cell_size(state->protein_conc, state, genes, *dt,
                                    rates, konStates, *t,
                                    timecoursestart, timecourselast,
                                    state->protein_conc);
      
      /* the number of mRNAs in cytoplasm affects binding */
      change_mRNA_cytoplasm(i, genes, state, rates, konStates);
      break;
    case 3:  /* replicate gene */
      *dt = state->replication_time_end->time - *t;         /* make dt window smaller */
      
      replicate_gene(state, genes, rates, konStates, koffvalues, state->replication_time_end->gene_id, *t);
      
      /* delete the event that just happened */
      delete_fixed_event_start(&(state->replication_time_end), &(state->replication_time_end_last));
      
      update_protein_conc_cell_size(state->protein_conc, state, genes, *dt,
                                    rates, konStates, *t,
                                    timecoursestart, timecourselast,
                                    state->protein_conc);
      break;
    default:
      printf("event=%d should never get here\n", event);
      exit(-1);
      break;
    }
    
    /* advance time by the dt */
    *t += *dt;
    *x -= (*dt)*(*konrate);

    LOG_VERBOSE("dt=%g t=%g fixed event old x=%g new x=%g\n", *dt, *t, (*x)+(*dt)*(*konrate), *x);
    
    /* re-compute a new dt */
    calc_dt(x, dt, rates, konStates, mRNAdecay, 
            genes->mRNAdecay, state->mRNA_cyto_num, state->mRNA_transl_cyto_num, state->cell_id);
    
    LOG_VERBOSE("next stochastic event (2) due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

    fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);    

    /* check to see there aren't more fixed events to do */
    event = does_fixed_event_end(state->mRNA_transl_time_end, 
                                 state->mRNA_transcr_time_end, 
                                 state->replication_time_end,
                                 fixed_time);
  } 

  /* no remaining fixed events to do in dt, now do stochastic events */
  
  /* if we haven't already reached end of development with last delta-t */
  if (*t+*dt < tdevelopment || no_fixed_dev_time) {
    
    // TODO: check this logic!
    /* compute konrate (which is constant over the Gillespie step) */
    if (konStates->nkon==0) {
      *konrate = (-rates->salphc);    /* all binding sites are occupied, total rate of binding is zero */
      LOG_ERROR("konStates->nkon=%d, konrate=%g\n", konStates->nkon, *konrate);
    } else  {
      calc_kon_rate(*dt, konStates, konrate);           /* otherwise compute konrate */
    }

    // TODO: FIXME 
    // do sanity check on total rates first
    if (!(rates->total + *konrate > 0.0)) {
      log_snapshot(rates, state, genes, konStates, &koffvalues, mRNAdecay, transport, *konrate, *x, *t);
      LOG_ERROR("x should always be >0 t=%g (x=%g) rates->total=%g, konrate=%g, recalibrate cell!\n", *t, *x, rates->total, *konrate); 
      //recompute_kon_rates(rates, state, genes, konStates, 0);
      //recompute_koff_rates(rates, state, koffvalues, *t);
      recalibrate_cell(rates, state, genes, konStates, &koffvalues, mRNAdecay, transport, *dt); 
      calc_kon_rate(*dt, konStates, konrate);           /* also recompute konrate */
      log_snapshot(rates, state, genes, konStates, &koffvalues, mRNAdecay, transport, *konrate, *x, *t);

      // TODO: probably should mark cell as "dead" in this case,
      // i.e. no TFs bound, no activity etc., and remove from queue.
      if (!(rates->total + *konrate > 0.0)) {  // if still wrong, exit
        // TODO: hack time to always put the event at the very back of the queue: this cell is effectively dead and should eventually
        // be chosen by the random number to be replaced
        *t = TIME_INFINITY;
        LOG_ERROR("put event at the very end of queue, with t=%g, cell is effectively dead\n", *t); 
        return -1;
      }
    }
    
    /* 
     * choose a new uniform random number weighted by the
     * probability of all Gillespie events, note that konrate is
     * *not* included in rates->total, so it needs to be added here
     */
    *x = ran1(&seed)*(rates->total + *konrate);  
    
    if (verbose) {
      log_snapshot(rates,
                   state,
                   genes,
                   konStates,
                   &koffvalues,
                   mRNAdecay,
                   transport, 
                   *konrate,
                   *x,
                   *t);
    }
    /* JM: kon generally could be handled better, with more direct
     * references to konStates->nkonsum, probably a bit vulnerable to rounding
     * error */
    
    /* 
     * STOCHASTIC EVENT: a TF unbinds (koff) 
     */
    if (*x < rates->koff) {  
      tf_unbinding_event(rates, state, genes, konStates, koffvalues,
                         timecoursestart, timecourselast, (*konrate), *dt, *t, *x, 1);
    } else {
      *x -= rates->koff;  
      /* 
       * STOCHASTIC EVENT: a transport event
       */
      if (*x < rates->transport) {     
        transport_event(rates, state, genes, konStates, transport, 
                        timecoursestart, timecourselast, ttranslation, *dt, *t, *x);
      } else {
        
        *x -= rates->transport;
        /* 
         * STOCHASTIC EVENT: an mRNA decay event
         */
        if (*x < rates->mRNAdecay) {  
          mRNA_decay_event(rates, state, genes, konStates, mRNAdecay,
                           timecoursestart, timecourselast, *dt, *t, *x);
        } else {
          *x -= rates->mRNAdecay;
          /* 
           * STOCHASTIC EVENT: PIC disassembly
           */
          if (*x < rates->pic_disassembly) {
            disassemble_PIC_event(rates, state, genes, konStates, 
                                  timecoursestart,  timecourselast, *dt, *t, *x);
          } else {
            *x -= rates->pic_disassembly;
            /* 
             * STOCHASTIC EVENT: TF binding event
             */
            if (*x < rates->salphc + (*konrate)) {   /* add variable (salphc) and constant (konrate) */
              tf_binding_event(rates, state, genes, konStates, koffvalues,
                               timecoursestart, timecourselast, (*konrate), *dt, *t, 
                               maxbound2, maxbound3, 1);
            } else {
              *x -= (rates->salphc + (*konrate));
              /* 
               * STOCHASTIC EVENT: histone acetylation
               */
              if (*x < (float) sum_rate_counts(rates->acetylation_num) * acetylate) {
                
                histone_acteylation_event(rates, state, genes, konStates, 
                                          timecoursestart, timecourselast, *dt, *t);
              } else {
                
                *x -= (float) sum_rate_counts(rates->acetylation_num) * acetylate;
                /* 
                 * STOCHASTIC EVENT: histone deacetylation
                 */
                if (*x < (float) sum_rate_counts(rates->deacetylation_num) * deacetylate) {
                  
                  histone_deacteylation_event(rates, state, genes, konStates, 
                                              timecoursestart, timecourselast, *dt, *t);
                } else {
                  *x -= (float) sum_rate_counts(rates->deacetylation_num) * deacetylate;
                  /* 
                   * STOCHASTIC EVENT: PIC assembly
                   */
                  if (*x < (float) sum_rate_counts(rates->pic_assembly_num) * PICassembly) {
                    assemble_PIC_event(rates, state, genes, konStates, 
                                       timecoursestart, timecourselast, *dt, *t);
                  } else {
                    *x -= (float) sum_rate_counts(rates->pic_assembly_num) * PICassembly;
                    /* 
                     * STOCHASTIC EVENT: transcription initiation
                     */
                    if (*x < (float) sum_rate_counts(rates->transcript_init_num) * transcriptinit) {
                      transcription_init_event(rates, state, genes, konStates, 
                                               timecoursestart, timecourselast, *dt, *t, *x);
                    } else {
                      /*
                       * FALLBACK: shouldn't get here, previous
                       * events should be exhaustive
                       */
                      
                      LOG_ERROR("[cell %03d] t=%g no event assigned: x=%g, rates->total+konrate=%g, recalibrate cell\n", 
                                state->cell_id, *t, *x, rates->total + *konrate);

                      log_snapshot(rates, state, genes, konStates, &koffvalues, mRNAdecay,
                                   transport, *konrate, *x, *t);
                      recalibrate_cell(rates, state, genes, konStates, &koffvalues,
                                       mRNAdecay, transport, *dt); 
                      log_snapshot(rates, state, genes, konStates, &koffvalues, mRNAdecay,
                                   transport, *konrate, *x, *t);
                      if (*t > 1000.0) {
                        LOG_ERROR("cell has been running for %g which is > 1000 mins\n", *t);
                      }
                      // TODO:  2009-05-24 should probably break out of "if" statement
                      //  as we should actually ignore this timestep completely 
                      // and through out the dt we used, and not advance time until we have
                      // a successfull event
                      // break;
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
    update_protein_conc_cell_size(state->protein_conc, state, genes, *dt,
                                  rates, konStates,
                                  *t, timecoursestart, timecourselast,
                                  state->protein_conc);
    /* advance to end of development (this exits the outer while loop) */
    *t = tdevelopment;
  }
  return 0;
}

/*
 * develop: run the population of cells for a given number of
 * divisions or run cell(s) for a specific length of time
 */
//TODO: rename run_pop()?
void develop(Genotype genes[POP_SIZE],
             CellState state[POP_SIZE],
             TimeCourse *timecoursestart[POP_SIZE][NGENES],
             TimeCourse *timecourselast[POP_SIZE][NGENES],
             float temperature,   /* in Kelvin */
             float kdis[NUM_K_DISASSEMBLY],
             int hold_genotype_constant,
             int output_binding_sites,
             int no_fixed_dev_time,
             int max_divisions)
{
  /* local variables that don't require per-cell tracking */
  int i, j;
  int maxbound2, maxbound3;  
  int curr_seed;

  /* initial mRNA and protein concentrations */
  float initmRNA[NGENES], initProteinConc[NGENES];

  /* cached information about available binding sites for efficiency */
  KonStates konStates[POP_SIZE];
  GillespieRates rates[POP_SIZE];

  float t[POP_SIZE];
  float *koffvalues[POP_SIZE];   /* rates of unbinding */
  float transport[POP_SIZE][NGENES];  /* transport rates of each mRNA */
  float mRNAdecay[POP_SIZE][NGENES];  /* mRNA decay rates */
  float x[POP_SIZE];                  /* random number */
  float dt[POP_SIZE];                 /* delta-t */
  float konrate[POP_SIZE];

  /* priority queue and time step initialization */
  //FixedEvent *time_queue = NULL;
  //FixedEvent *time_queue_end = NULL;
  float t_next = 0.0;   /* initialize to zero */
  int cell = 0;         /* initialize cell number (TODO: this is only
                           necessary when running for fixed
                           development time) */
  int divisions = 0;    /* no cell divisions yet */
  int motherID;         /* mother cell at division */
  int daughterID;       /* daughter cell at division */

  int large_cell_size = 0;  // TODO: count the number of times cell_size exceeds Y=2.0

  /* keep a reaper queue from which to choose 'dead cells', if there are any, as daughter cells when dividing */
  int keep_reaper_queue = 1;  

  //long int timesteps = 0; 

  maxbound2 = maxbound;
  maxbound3 = 10*maxbound;

  /* initialize heap */
  bheap_t *queue, *empty_queue;

  queue = bh_alloc(POP_SIZE);
  empty_queue = bh_alloc(POP_SIZE);

  /* initialize protein concentrations to be used in all genes */
  // TODO: FIXME ASAP initProteinConc should be over NPROTEINS and
  // initmRNA should only loop over NGENES as they could be different
  for (i=0; i < NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }


  for (j = 0; j < POP_SIZE; j++) {

    output=1;
    /* for all cells in this replicate initialize all parts of the
       genotype, *except* for the cis-reg sequence using the initial
       genes[0]  */
    initialize_genotype(&genes[j], &genes[0], kdis, j);
    /* if genotype is held constant, start varying the seed *after*
       initialize_genotype, so we can run the same genotype with
       slightly different amounts of noise  */
    // TODO: maybe remove this?  needs to be put back at some point in more systematic way
    if (hold_genotype_constant)
      for (curr_seed=0; curr_seed<dummyrun; curr_seed++) 
        ran1(&seed);
   
    initialize_cell(&state[j], j, genes[j].copies, genes[j].mRNAdecay, initmRNA, initProteinConc, burn_in);

    /* print binding sites */
    if (output_binding_sites) 
      print_all_binding_sites(genes[j].copies, genes[j].all_binding_sites, genes[j].binding_sites_num, 
                              genes[j].tf_seq, genes[j].cisreg_seq, genes[j].site_id_pos); 

    
    /* set cell temperature and value of RTlnKr constant */
    state[j].temperature = temperature;
    state[j].RTlnKr = GasConstant * temperature * log(Kr);

    /* initialize time courses */
    for (i=0; i < NPROTEINS; i++){
      timecoursestart[j][i] = NULL;
      timecourselast[j][i] = NULL;
    } 

    add_time_points((float) 0.0, state[j].protein_conc, timecoursestart[j], timecourselast[j]);

    /* initialize konStates data structures */
    initialize_cell_cache(&(state[j]), genes[j], &(konStates[j]), &(koffvalues[j]), maxbound2, maxbound3);

    /* initialize transcriptional state of genes */
    calc_from_state(&genes[j], &state[j], &rates[j], &konStates[j], transport[j], mRNAdecay[j]);

    t[j] = 0.0;  /* time starts at zero */
    state[j].division_time = TIME_INFINITY;  /* make artificially high */

    /* do one step for the cell */
    do_single_timestep(&(genes[j]), 
                       &(state[j]), 
                       &(konStates[j]), 
                       &(rates[j]), 
                       &(t[j]),
                       koffvalues[j],
                       transport[j],
                       mRNAdecay[j],
                       &(x[j]),
                       &(dt[j]),
                       &(konrate[j]),
                       timecoursestart[j],
                       timecourselast[j],
                       maxbound2,
                       maxbound3,
                       no_fixed_dev_time);
    
    /* insert each initial time into the queue */
    int ops;
    ops = insert_with_priority_heap(queue, j, t[j]);
    LOG_NOCELLID("[cell=%03d] inserted at time=%g\n", j, t[j]); 
  }

  while ((no_fixed_dev_time && divisions < max_divisions) ||    /* no fixed dev time, run until # divisions reached */
         (!no_fixed_dev_time && t_next < tdevelopment)) {       /* or, if fixed dev time, run until tdevelopment reached */

    /* get the next cell with the smallest t to advance next */
    int ops, retval;  // TODO: remove ops references

    LOG_VERBOSE_NOCELLID("before choosing next cell (queue len reaper=%03d, main=%03d)\n", empty_queue->n, queue->n);

    if (empty_queue->n == POP_SIZE || queue->n == 0) {
      /* if all cells in population are dead, we exit */
      LOG_NOCELLID("all %03d cells in population are dead, main queue len=%03d\n", empty_queue->n, queue->n);
      break;
    } 
    /* otherwise get next event from queue */
    t_next = get_next_heap(queue, &cell, &ops);
    LOG_VERBOSE_NOCELLID("[cell=%03d] get minimum time=%g, size=%g (queue len reaper=%03d, main=%03d)\n", 
                         cell, t_next, state[cell].cell_size, empty_queue->n, queue->n);

    retval = do_single_timestep(&(genes[cell]), 
                                &(state[cell]), 
                                &(konStates[cell]), 
                                &(rates[cell]), 
                                &(t[cell]),
                                koffvalues[cell],
                                transport[cell],
                                mRNAdecay[cell],
                                &(x[cell]),
                                &(dt[cell]),
                                &(konrate[cell]),
                                timecoursestart[cell],
                                timecourselast[cell],
                                maxbound2,
                                maxbound3,
                                no_fixed_dev_time);

    // TODO: better variable name than retval
    if (retval == -1)  {
      // TODO: check, print out protein time courses when a cell is found to be dead
      print_all_protein_time_courses(timecoursestart, timecourselast);
      /* add this cell to list of empty cell locations */
      // TODO: make this configurable at run-time
      if (keep_reaper_queue)
        insert_with_priority_heap(empty_queue, cell, TIME_INFINITY);
      LOG_NOCELLID("[cell %03d] added here to empty_queue as a dead cell, t_next=%g, t=%g (queue len reaper=%03d, main=%03d)\n", 
                   cell, t_next, t[cell], empty_queue->n, queue->n);
      /* this cell is dead, so don't add back to main priority queue, skip to next event */
      continue;
    }
   
    // TODO: cleanup
    /* abort if a population of one cell, undergoing a single division
       exceeds an maximum upper limit */
    if (POP_SIZE==1 && max_divisions==1 && timemax > 0.0) {
      if (t[cell] > timemax) {
        printf("[cell %03d] at t=%g exceeds the maximum time of t=%g\n", 
               cell, t[cell], timemax);
        LOG_ERROR("[cell %03d] at t=%g exceeds the maximum time of t=%g\n", 
                  cell, t[cell], timemax);
        break;
      }
    }

    // TODO: generate new regression output based on choosing the
    // event that just happened, rather than the previous t_next event, then switch lines
    //if (t[cell] >= (state[cell]).division_time)  {  /* we have now reached cell division */
    if (t_next >= (state[cell]).division_time)  {  /* we have now reached cell division */
      float current_division_time = (state[cell]).division_time;  // store current division time
      divisions++;     /* increment number of divisions */
      motherID = cell; /* set mother cell as currently dividing cell */

      if (state[cell].cell_size > 2.0) {
        large_cell_size++;
        printf("[cell %03d] at t=%g cell size =%g, exceeding Y=2.0, %d of %d divisions (%g fraction)\n", 
               cell, t[cell], state[cell].cell_size, large_cell_size, divisions, (double)large_cell_size/(double)divisions);

        LOG_NOCELLID("[cell %03d] at t=%g cell size =%g, exceeding Y=2.0, %d of %d divisions (%g fraction)\n", 
                     cell, t[cell], state[cell].cell_size, large_cell_size, divisions, (double)large_cell_size/(double)divisions);
      }

      /* to get new daughter cell slot, first check list of empty cells */
      if (empty_queue->n > 0) {
        get_next_heap(empty_queue, &daughterID, &ops);  /* use one of the empty cells */
        LOG_NOCELLID("[cell %03d] in empty_queue (length=%03d) used as daughter cell\n", daughterID, empty_queue->n);
      } else {
        daughterID = rint((POP_SIZE-1)*ran1(&seed));  /* otherwise choose one other cell randomly */
        /* removing pending event in daughter cell from queue, if different from the mother 
           only remove if the daughter cell wasn't already removed from queue */
        if (daughterID != motherID) {
          delete_element_heap(queue, daughterID);
          LOG_NOCELLID("removing pending event from queue (length=%3d) replaced by daughter cell=%d\n", queue->n, daughterID);
        } 
      }

      printf("[cell %03d] (size=%g) dividing into mother=%03d and daughter=%03d at t=%g, division=%g, total divisions=%d (t_next=%g)\n", 
             cell, state[cell].cell_size, motherID, daughterID, t[cell], current_division_time, divisions, t_next);
      LOG_NOCELLID("[cell %03d] (size=%g) dividing into mother=%03d and daughter=%03d at t=%g, division=%g, total divisions=%d (t_next=%g)\n", 
                   cell, state[cell].cell_size, motherID, daughterID, t[cell], current_division_time, divisions, t_next);

      if (time_s_phase + time_g2_phase > 0.0) {
        do_cell_division(motherID,
                         daughterID,
                         &(genes[motherID]),
                         &(state[motherID]),
                         &(rates[motherID]),
                         &(konStates[motherID]),
                         &(koffvalues[motherID]),
                         transport[motherID],
                         mRNAdecay[motherID],
                         &(genes[daughterID]),
                         &(state[daughterID]),
                         &(rates[daughterID]),
                         &(konStates[daughterID]),
                         &(koffvalues[daughterID]),
                         transport[daughterID],
                         mRNAdecay[daughterID],
                         0.44,
                         x[motherID],
                         dt[motherID]);
      } else {
        LOG_NOCELLID("Skip division for cell=%03d S phase and G2 phase have zero length\n", motherID);
      }
        
      if (daughterID != motherID) {
        t[motherID] = current_division_time;   // reset current time in mother cell to division time
        LOG_NOCELLID("AFTER DIVISION: time at instant of division for mother cell=%03d is t=%g\n", motherID, t[motherID]);
        /* advance time for mother */
        retval = do_single_timestep(&(genes[motherID]), 
                           &(state[motherID]), 
                           &(konStates[motherID]), 
                           &(rates[motherID]), 
                           &(t[motherID]),
                           koffvalues[motherID],
                           transport[motherID],
                           mRNAdecay[motherID],
                           &(x[motherID]),
                           &(dt[motherID]),
                           &(konrate[motherID]),
                           timecoursestart[motherID],
                           timecourselast[motherID],
                           maxbound2,
                           maxbound3,
                           no_fixed_dev_time);

        if (retval == -1) {
          if (keep_reaper_queue)
            insert_with_priority_heap(empty_queue, motherID, TIME_INFINITY);
          LOG_NOCELLID("AFTER DIVISION: mother cell %03d is dead at t=%g (queue len reaper=%03d, main=%03d)\n", 
                       motherID, t[motherID], empty_queue->n, queue->n);
        } else {
          LOG_NOCELLID("AFTER DIVISION: add new timestep to queue for mother cell=%03d at t=%g (queue len reaper=%03d, main=%03d)\n", 
                       motherID, t[motherID], empty_queue->n, queue->n);
          insert_with_priority_heap(queue, motherID, t[motherID]); 
        }
        
      } else {
        printf("mother (%3d) and daughter (%3d) are the same: don't update mother it is replaced by daughter\n",
               motherID, daughterID);
        LOG_NOCELLID("mother (%3d) and daughter (%3d) are the same: don't update mother it is replaced by daughter\n",
            motherID, daughterID);
      }

      t[daughterID] = current_division_time;   // reset current time in mother cell to division time
      LOG_NOCELLID("AFTER DIVISION: time at instant of division for daughter cell=%03d is t=%g\n", daughterID, t[daughterID]);

      /* advance time for daughter */
      retval = do_single_timestep(&(genes[daughterID]), 
                                  &(state[daughterID]), 
                                  &(konStates[daughterID]), 
                                  &(rates[daughterID]), 
                                  &(t[daughterID]),
                                  koffvalues[daughterID],
                                  transport[daughterID],
                                  mRNAdecay[daughterID],
                                  &(x[daughterID]),
                                  &(dt[daughterID]),
                                  &(konrate[daughterID]),
                                  timecoursestart[daughterID],
                                  timecourselast[daughterID],
                                  maxbound2,
                                  maxbound3,
                                  no_fixed_dev_time);

      if (retval == -1) {
        if (keep_reaper_queue)
          insert_with_priority_heap(empty_queue, daughterID, TIME_INFINITY);
        LOG_NOCELLID("AFTER DIVISION: daughter cell %03d is dead at t=%g (queue len reaper=%03d, main=%03d)\n", 
                     daughterID, t[daughterID], empty_queue->n, queue->n);
      } else {
        LOG_NOCELLID("AFTER DIVISION: add new timestep in queue for daughter cell=%03d at t=%g (queue len reaper=%03d, main=%03d)\n", 
                     daughterID, t[daughterID], empty_queue->n, queue->n);
        insert_with_priority_heap(queue, daughterID, t[daughterID]);
      }
    } else {
      /* put the updated timestep back into the queue  */
      ops = insert_with_priority_heap(queue, cell, t[cell]);
    }
  }

  /* output the founder information */
  for (j = 0; j < POP_SIZE; j++) {
    printf("cell %03d derived from founder %03d had %2d divisions%s\n", 
           j, state[j].founder_id, state[j].divisions, t[j] == TIME_INFINITY ? " [dead]": "");
    LOG_NOCELLID("cell %03d derived from founder %03d had %2d divisions%s\n", 
                 j, state[j].founder_id, state[j].divisions, t[j] == TIME_INFINITY ? " [dead]": "");
  }
  // TODO: check output contents of reaper queue
  bh_dump(empty_queue);

  /* cleanup data structures */
  bh_free(queue);
  bh_free(empty_queue);

  for (j = 0; j < POP_SIZE; j++) {
    free(koffvalues[j]);
    for (i=0; i < NPROTEINS; i++) {
      free(konStates[j].kon_list[i]->available_sites);
      free(konStates[j].kon_list[i]);
    }

    // TODO: remove?
    /* free(&konStates[j]);
       free(&rates[j]); */
  }
}

void print_time_course(TimeCourse *start,
                       int i,
                       int j)
{
  FILE *fpout;
  char filename[80];
  
  /* do the normal thing on the first cell */
  if (POP_SIZE == 1)
    sprintf(filename, "%s/protein%d.dat", output_directory, i);
  else
    sprintf(filename, "%s/protein%03d-%02d.dat", output_directory, j, i);
  if ((fpout = fopen(filename,"w"))==NULL) {
    LOG_ERROR_NOCELLID("error: Can't open %s file\n", filename);
  }
  while (start) {
    fprintf(fpout,"%g %g\n", start->time, start->concentration);
    start = start->next;
  }
  fclose(fpout);  
}

void print_all_protein_time_courses(TimeCourse *timecoursestart[POP_SIZE][NPROTEINS],
                                    TimeCourse *timecourselast[POP_SIZE][NPROTEINS])
{
  int i, j;
  for (j = 0; j < POP_SIZE; j++) {
    for (i=0; i < NPROTEINS; i++) {
      if ((output)) print_time_course(timecoursestart[j][i], i, j);
    }
  }
}

// TODO: move to utility library lib.c ?

void create_output_directory(char *output_directory) {
  int directory_success;

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

}

void create_output_file(char prefix[80], char *output_directory, FILE **fp, int index) {
  char file_name[80];
  if (index != -1) 
    sprintf(file_name, "%s/%s-%03d.dat", output_directory, prefix, index);
  else    /* if index is -1, use prefix as name of file, unadorned with .dat */
    sprintf(file_name, "%s/%s", output_directory, prefix);
  if ((*fp = fopen(file_name,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", file_name);
}

void read_kdisassembly(float kdis[NUM_K_DISASSEMBLY]) {
  int j;
  FILE *fpkdis;
  /* get the kdis.txt values */
  if ((fpkdis = fopen("kdis.txt","r"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", "kdis.txt");
  for (j = 0; j < NUM_K_DISASSEMBLY; j++) {
    fscanf(fpkdis, "%f", &kdis[j]);
  }
  fclose(fpkdis);
}
