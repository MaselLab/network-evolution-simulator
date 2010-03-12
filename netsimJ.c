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
#include <assert.h>
#include <time.h>

/* local includes */
#include "randomJ.h"
#include "libJ.h"
#include "priority-queueJ.h"
#include "netsimJ.h"

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
const float NITER = 100;

const float mN = 3.3e-10 * 5e6;    /* Lynch et al. (2008), Tsai et al. (2008)  */

/* below are default options, can be changed on command line*/

float kon=1e-4;              /* lower value is so things run faster */
                             /* actual value should be kon=0.2225 is
                                based on 1 molecule taking
                                240seconds=4 minutes and 89% of the
                                proteins being in the nucleus*/
float kon_after_burnin=1e-4; /* lower value value after burn is so things run faster */

int burn_in = 0;             /* disable burn-in by default */

float tdevelopment = 30.0;/* default  development time: can be changed at runtime */
float timemax = -1.0;      /* set an upper limit to development time (default to -1.0=no limit) */
int current_ploidy = 1;    /* ploidy can be changed at run-time: 1 = haploid, 2 = diploid */
int output = 0;
long seed = 28121; //94571       /* something is wrong here: changing seed changes nothing */
int dummyrun = 7;          /* used to change seed */
int recompute_koff = 1;    /* toggle whether to recompute certain features at each time to avoid
                              compounding rounding error (off by default) */
int recompute_kon = 1;     /* likewise for kon */
float critical_size = 1.0; /* critical size at which cell divides, 
                              set to negative to prevent division  */
float growth_rate_scaling = 2.0; /* set default growth rate scaling factor */
float time_s_phase = 30.0;  /* length of S phase (default: 30 mins) */
float time_g2_phase = 30.0; /* length of G2/mitosis phase (default: 30 mins) */
int random_replication_time = 1; /* replication times in S phase are random by default */
//float *gene_active = malloc(NGENES*sizeof(float));

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

  /* for each gene determine a time point during the [0, time_s_phase min] S-phase interval */
  for (i=0; i < NGENES; i++) {
    if (random_replication_time) 
      genotype->replication_time[i] = time_s_phase*ran1(&seed); 
    else
      genotype->replication_time[i] = time_s_phase*(i/(float)NGENES);
    LOG_VERBOSE_NOCELLID("[genotype %03d] offset for replication time after S-phase starts: %g\n", genotype_id, genotype->replication_time[i]);
  }
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
  if (genotype_id == 0) {
    initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*MAX_COPIES*TFGENES, MAX_COPIES, TFGENES);
    initialize_genotype_fixed(genotype, kdis, genotype_id);
  // otherwise clone it from the first instance
  } else {
    /* copy the TF sequence data */
    memcpy(genotype->tf_seq, clone->tf_seq, sizeof(char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN]));
    initialize_new_cell_genotype(genotype, clone);
  }

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
 //DESTROY??
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
                     int cell_id,
                     int copies[NGENES],
                     float mRNAdecay[NGENES],
                     float meanmRNA[NGENES],
                     float init_protein_conc[NPROTEINS],
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
}


/*
 * initialize the cell state that are "cached", i.e. that are
 * maintained for code efficiency and are not part of the underlying
 * molecular biology.  All of these data structures can be regenerated
 * from the cell state at any point in the simulation
 */

//CHANGE JU
void calc_time (float t, 
                float x, 
                GillespieRates *rates,
                //KonStates *kon_states,
                float *f, 
                float *df)
{

 *f = 1.0e-001;

  LOG_VERBOSE_NOCELLID("x=%g,  rates->subtotal=%g,  f=%g\n", x, rates->subtotal, *f);

  // compute derivative of equation 
  //*df = x*numer/denom - 1.0;
  *df =  -1.0;

  LOG_VERBOSE_NOCELLID(" t=%g f=%g df=%g\n", t, *f, *df);
}





/* 
 * for specified gene_id and gene_copy, tests whether criterion for
 * transcription is met
 */
 //DESTROY
 //CHANGE JU
 int ready_to_transcribeJU(int gene_id, float *gene_active, CellState *state){
    //srand(time(NULL));
    float n;
    n = ran1(&seed);
    //n= rand()%10 +1;
    //n=(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
    //n= rand(0,1);
     //float m = n%1000;
     //m= m/1000;
     LOG_ERROR("n = %f, gene active = %f\n", n, gene_active[gene_id] );
    // printf("rand = %f,  gene_active = %f\n", m, gene_active[gene_id]);
    //LOG_VERBOSE("m = %f,  gene_active = %f\n", m, gene_active[gene_id]);
     if(n < gene_active[gene_id]){
        return(1);
     }else{
        return(0); 
     }
 }
 


/* 
 * calculate the rates based upon the initialization of the genotype
 * and the cell state.  Note this is only appropriate if nothing is
 * bound
 */
 //DESTROY?
void calc_from_state(Genotype *genotype,
                     CellState *state,
                     GillespieRates *rates,
                     //KonStates *kon_states,
                     float transport[],
                     float mRNAdecay[])
{
  int i, j, k;
  float salphc, protein_conc_for_tfid; 
  float protein_decay;
//QUESTION
  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to fall below cut-off, use aging term */
    protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;
    salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);
    
  }
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

  for (k=0; k < genotype->binding_sites_num; k++) {
    i = genotype->all_binding_sites[k].tf_id;
    protein_conc_for_tfid = state->protein_conc[i];
  }

  for (i=0; i < NGENES; i++) {
    // FIXME: this updates kon_states for non-TFs, even though this isn't used
    LOG_VERBOSE("after initializing kon_states for gene=%d, sites_per_gene=%d\n", 
        i, genotype->sites_per_gene[i]);

    transport[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
    rates->transport += transport[i];
    rates->transport_operations++;
    LOG_VERBOSE("mRNA_nuclear_num=%d, initializing transport[%d]=%g\n", state->mRNA_nuclear_num[i], i, transport[i]);
  }
  LOG_VERBOSE("initializing rates->transport=%g\n", rates->transport);

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

  if (1) //VERBOSE JU CHANGE
    for (j=0; j < MAX_COPIES; j++) {
      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      //LOG_ERROR("rates->deacetylation_num[%d]=%d\n", j, rates->deacetylation_num[j]);
   
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
  rates->transcript_init_num[i] = 0;

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

void revise_prob_vect(int gene_id,
                           int gene_copy,
                           Genotype *genotype,
                           CellState *state,
                           GillespieRates *rates,
                           float *gene_active)
{
 int check;
 check = ready_to_transcribeJU(gene_id, gene_active, state);
}

/*
 * update the activity state of 'gene_id' based on the new balance of
 * bound transcription factors and the current state of the gene
 */
void revise_activity_state(int gene_id,
                           int gene_copy,
                           Genotype *genotype,
                           CellState *state,
                           GillespieRates *rates,
                           float *gene_active)
{
  int transcriptrule, oldstate ; //numactive

 /* transcriptrule = ready_to_transcribe(gene_id, gene_copy, 
                                       //state->tf_bound_indexes, 
                                       //state->tf_bound_num,
                                       genotype->all_binding_sites,
                                       genotype->activating,
                                       &numactive);*/
   //srand(time(NULL));
   transcriptrule = ready_to_transcribeJU(gene_id, gene_active, state);
   /*float test = ran1(&seed);
   if(test <= .5) transcriptrule = 1;
   else transcriptrule =0;*/
   LOG_VERBOSE("transcript rule = %d\n", transcriptrule);
  /* get last state of transcription initiation */
  oldstate = state->active[gene_id][gene_copy];
  LOG_ERROR("oldstate = %d\n", oldstate);
  // printf("gene id = %d\n", gene_id);
  //if(gene_id == 0){
             //system("PAUSE");
  //}
  /*
   * first set of rules:
   * ACTIVATING TFs exceed REPRESSING TFs 
   */

  /* OFF_FULL -> ON_WITH_NUCLEOSOME  P(A) * P(AC)*/
  LOG_ERROR("transcript rule = %d, oldstate = %d,  OFF_FULL = %d\n", transcriptrule, oldstate, OFF_FULL);
  if ((transcriptrule) && oldstate==OFF_FULL){
                       // printf("A\n");
                       //LOG_ERROR("A\n");
                        LOG_ERROR("UP ACE NUM!\n");
    state->active[gene_id][gene_copy] = ON_WITH_NUCLEOSOME;
    state->state_change_ids[ACETYLATION_STATE][gene_copy][rates->acetylation_num[gene_copy]] = gene_id;
    //LOG_ERROR("ace num = %d, gene copy = %d\n", rates->acetylation_num[gene_copy], gene_copy);
    (rates->acetylation_num[gene_copy])++;
      LOG_ERROR("ace num = %d\n", rates->acetylation_num[gene_copy]);
  }
  //LOG_ERROR("after first if\n");
 // LOG_ERROR("transcript rule = %d, oldstate = %d, OFF_NO_PIC = %d\n", transcriptrule, oldstate, OFF_NO_PIC);
  /* OFF_NO_PIC -> ON_NO_PIC  P(A)*P(DEAC)*/
  if ((transcriptrule) && oldstate==OFF_NO_PIC) {
                       // printf("B\n");
                       LOG_ERROR("OFF to ON\n");
    state->active[gene_id][gene_copy] = ON_NO_PIC;
    remove_from_array(gene_id, DEACETYLATION_STATE,  state->state_change_ids[DEACETYLATION_STATE][gene_copy], 
                      &(rates->deacetylation_num[gene_copy]), (int) 1);
    //P(A)
    if (transcriptrule){// FIGURE OUT WHAT THIS IS/DOES!!
    LOG_ERROR("Increase PIC assembly num\n");
    LOG_ERROR(" Before state change num = %d\n", rates->pic_assembly_num[gene_copy]);
      //if(numactive?){
      //state->state_change_ids[PICASSEMBLY_STATE][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id;
      //(rates->pic_assembly_num[gene_copy])++;
     // LOG_ERROR(" num = %d\n", rates->pic_assembly_num[gene_copy]); 
     // }
    }
  }
  /* OFF_PIC -> ON_FULL P(A)*/
  if ((transcriptrule) && oldstate==OFF_PIC) {
                       // printf("C\n");
                       LOG_ERROR("C\n");
    state->active[gene_id][gene_copy] = ON_FULL;
  }

  /*
   * second set of rules:
   * REPRESSING TFs exceed ACTIVATING TFs 
   */

  /* ON_WITH_NUCLEOSOME -> OFF_FULL P(R)*P(A) */
  if (!(transcriptrule) && oldstate==ON_WITH_NUCLEOSOME) {
                       //  printf("D\n");
                       LOG_ERROR("D\n");
    state->active[gene_id][gene_copy] = OFF_FULL;
    LOG_VERBOSE("removing gene=%d, copy=%d from state_change_ids[ACETYLATION_STATE][%d]\n", gene_id, gene_copy, gene_copy);
    remove_from_array(gene_id, ACETYLATION_STATE, state->state_change_ids[ACETYLATION_STATE][gene_copy], 
                      &(rates->acetylation_num[gene_copy]), (int) 1);
  }
  
  /* ON_NO_PIC -> OFF_NO_PIC P(R)*P(DEAC) */
  if (!(transcriptrule) && oldstate==ON_NO_PIC){
                        // printf("E\n");   
                       LOG_ERROR("E\n");       
    state->active[gene_id][gene_copy] = OFF_NO_PIC;
    remove_from_array(gene_id, PICASSEMBLY_STATE, state->state_change_ids[PICASSEMBLY_STATE][gene_copy], 
                      &(rates->pic_assembly_num[gene_copy]), (int) 0);
    state->state_change_ids[DEACETYLATION_STATE][gene_copy][rates->deacetylation_num[gene_copy]] = gene_id;
    (rates->deacetylation_num[gene_copy])++;
    LOG_ERROR("deace num = %d\n", rates->deacetylation_num[gene_copy]);
  }
  /* ON_FULL -> OFF_PIC P(R)*/
  if (!(transcriptrule) && oldstate==ON_FULL) {
                        // printf("F\n");
                        LOG_ERROR("F\n");
    state->active[gene_id][gene_copy] = OFF_PIC;
  }
//system("PAUSE");
  /* do remaining transitions:
   * OFF_PIC -> OFF_NO_PIC
   * ON_FULL -> ON_NO_PIC 
   */
  if ((state->active[gene_id][gene_copy]==OFF_PIC || state->active[gene_id][gene_copy]==ON_FULL) ){//&& numactive==0
  LOG_ERROR("INSIDE IF julu\n");
    disassemble_PIC(state, genotype, gene_id, gene_copy, rates);
}

  if (verbose && (oldstate!=state->active[gene_id][gene_copy])) {
    LOG_VERBOSE("state change from %d to %d in gene %d, copy %d\n", oldstate, state->active[gene_id][gene_copy], gene_id, gene_copy);
  }
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
 * when S phase is reached, we set the cell state and add the
 * genotypically-specified replication times to the queue
 */
void reach_s_phase(CellState *state, 
                   Genotype *genotype, 
                   float t) 
{
  int i;

  /* set the state as being in S phase */
  state->in_s_phase = 1;
  /* only add replication events if S phase and G2 phase have non-zero length
     otherwise just jump to division */
  if (time_s_phase + time_g2_phase > 0.0) {
    for (i=0; i < NGENES; i++) {
      LOG("add replication for gene_id=%d at t=%g\n", i, t + genotype->replication_time[i]);
      /* add the replication event to the queue */
      add_fixed_event(i, -1, t + genotype->replication_time[i], 
                      &(state->replication_time_end), &(state->replication_time_end_last));
    }
  }
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
                       float ect, 
                       float ect1) 
{
  return 1.0/(pow(c,2)*Pp) * gmax * (-alpha*ect1*s_mRNA + c*(P*ect1 + alpha*deltat*s_mRNA));
}

/*
 * return the instantaneous growth rate given the current cell state,
 * also return the integrated growth rate as a pointer
 */
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
  float instantaneous_growth_rate;  /* this is returned from the function */
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
    // printf("%f  ",all_alpha[i]);
    //LOG_VERBOSE_NOCELLID("all_alpha[%d] = %f\n", i, all_alpha[i]);
    //LOG_VERBOSE_NOCELLID("all_s_mRNA[%d] = %d\n", i, all_s_mRNA[i]);
    total_alpha_s += all_alpha[i] * all_s_mRNA[i];
  }
//system("PAUSE");
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
  fprintf(fp_growthrate[cell_id], "%g %g %g %g %g %g\n", t, instantaneous_growth_rate, *integrated_growth_rate, P, s_mRNA, c);
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
 //DESTROY
void update_protein_conc_cell_size(float protein_conc[],
                                   CellState *state,
                                   Genotype *genotype,
                                   float dt,
                                   GillespieRates *rates,
                                   //KonStates *kon_states,
                                   float t,
                                   TimeCourse **timecoursestart,
                                   TimeCourse **timecourselast,
                                   float otherdata[])
{
  int i;
  float ct, ect, ect1, salph;
  float L, L_next;
  float instantaneous_growth_rate = 0.0;
  float integrated_growth_rate = 0.0;
  float adjusted_decay;
  float protein_decay;

  //printf("translation = %f\n", genotype[0].translation[0]);
  //rates->max_salphc = rates->min_salphc = 0.0;
  for (i=0; i < NPROTEINS; i++) {
    if (i == SELECTION_GENE)  /* if we are looking at the selection gene, record protein concentration before update */
      L = protein_conc[i];

    /* update protein decay rates due to dilution caused by growth */
    adjusted_decay = genotype->proteindecay[i] + state->growth_rate;
  //printf("protein_decay = %f  adjusted decay = %f,  protein ageing = %f\n", genotype->proteindecay[i] ,adjusted_decay,  protein_aging);
    /* if this results in a very small or zero decay rate, use protein aging term */
   if (adjusted_decay > protein_aging)
      protein_decay = adjusted_decay;
    else 
      protein_decay = protein_aging;

    /* print out warning if decay rates get too low */
    if (protein_decay < 1e-10) {
      LOG_WARNING("protein=%02d, protein_decay=%g, genotype->proteindecay=%g, protein_aging=%g\n", i, adjusted_decay, 
                  genotype->proteindecay[i], protein_aging);
    }

    //LOG_VERBOSE("prot decay[%d]=%g\n", i, kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX]);

    ct = protein_decay*dt;
    //ct = dt;
    ect = exp(-ct);
    if (fabs(ct)<EPSILON) ect1=ct;
    else ect1 = 1-ect;   

//printf("protein decay = %f  ect =%f \n",protein_decay, ect);
//system("PAUSE");
    /* get the new protein concentration for this gene */
    //printf("mRNA nuclear = %d  \n",state->mRNA_cyto_num[i]);
    //system("PAUSE");
    salph = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i]/protein_decay;
   //printf("\n salph = %f,  ect1 = %f,  ect = %f, protein_conc[i] = %f\n", salph, ect1, ect, protein_conc[i]);
    //rates->salphc = salph;
    protein_conc[i] = salph*ect1 + ect*protein_conc[i];
    //printf("protein_conc[i] = %f\n", protein_conc[i]);
   
    if (i == SELECTION_GENE) {  /* now find out the protein concentration at end of dt interval and compute growth rate */
      L_next = protein_conc[i];
      instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genotype->translation[SELECTION_GENE], state->mRNA_cyto_num[SELECTION_GENE], 
                                                            genotype->translation, state->mRNA_cyto_num, 
                                                            L, L_next, t, dt, 1,
                                                            //kon_states->konvalues[SELECTION_GENE][KON_PROTEIN_DECAY_INDEX], 
                                                            ect, ect1,
                                                            state->cell_id);
      /* us the integrated growth rate to compute the cell size in the next timestep */
      state->cell_size = (state->cell_size)*exp(integrated_growth_rate);
      fprintf(fp_cellsize[state->cell_id], "%g %g\n", t, state->cell_size);
    }

  }
  /* scale up the rates using kon global */
  //rates->max_salphc *= kon;
  //rates->min_salphc *= kon;
  if ((output) && (*timecourselast)->time < t+dt-0.1) 
    add_time_points(t+dt, otherdata, timecoursestart, timecourselast);

  /* update the instantaneous growth rate for the beginning of the next timestep */
  state->growth_rate = instantaneous_growth_rate;
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
    //LOG_ERROR_NOCELLID("%d  ", rate_array[i]);
  }
  //LOG_ERROR_NOCELLID("\n");
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
                     //KonStates *kon_states,
                     float transport[NGENES],
                     TimeCourse **timecoursestart, 
                     TimeCourse **timecourselast, 
                     float dt,
                     float t,
                     float x)
{
  int i,j;
  float temp_rate;
  float endtime = t + dt + TTRANSLATION; 

    for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("?rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      }

  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, rates, 
                                //kon_states, 
                                t, timecoursestart, timecourselast, 
                                state->protein_conc);
  
  i = -1;
  temp_rate = 0.0;  
  /* choose gene product (mRNA) that gets transported to cytoplasm
     based on weighting in transport[] array */
  while (i < NGENES && x > temp_rate) {
    i++;
    x -= transport[i];
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

  transport[i] -= KRNA;   /* decrease transport frequency */

  /* if below a threshold, make zero */
  if (transport[i] < 0.1*KRNA) 
    transport[i]=0.0;

  /* adjust rates */
  rates->transport -= KRNA;
  rates->transport_operations++;

  /* do similar threshold check */
  if (rates->transport < 0.1*KRNA) 
    rates->transport=0.0;
}





void mRNA_decay_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                      //KonStates *kon_states, 
                      float *mRNAdecay, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t, float x)
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
                                rates, 
                                //kon_states, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  /* decay mRNA in cytoplasm */
  if (x < (float)state->mRNA_cyto_num[i]) {
    LOG_VERBOSE("mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
                i, state->mRNA_cyto_num[i], state->mRNA_transl_cyto_num[i]);
    
    /* remove the mRNA from the cytoplasm count */
    (state->mRNA_cyto_num[i])--;  
    //change_mRNA_cytoplasm(i, genotype, state, rates); //, kon_states
    
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
                               //KonStates *kon_states, 
                               TimeCourse **timecoursestart, TimeCourse **timecourselast,
                               float dt, float t)
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
                                rates, //kon_states, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  
  /* set state: eject nucleosome, but there is no PIC yet */
  state->active[gene_id][gene_copy] = ON_NO_PIC;
  remove_from_array(gene_id, ACETYLATION_STATE, state->state_change_ids[ACETYLATION_STATE][gene_copy], &(rates->acetylation_num[gene_copy]), (int) 1);
 // state->state_change_ids[PICASSEMBLY_STATE][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id; 
    //(rates->pic_assembly_num[gene_copy])++;
  //if (is_one_activator(gene_id, gene_copy, //state->tf_bound_indexes, state->tf_bound_num, 
                      // genotype->all_binding_sites, genotype->activating)) {
    state->state_change_ids[PICASSEMBLY_STATE][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id; 
    (rates->pic_assembly_num[gene_copy])++;
  //}
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                                 //KonStates *kon_states, 
                                 TimeCourse **timecoursestart, TimeCourse **timecourselast,
                                 float dt, float t)
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
                                rates, //kon_states, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  /* set state: nucleosome returns */
  state->active[gene_id][gene_copy] = OFF_FULL;
  remove_from_array(gene_id, DEACETYLATION_STATE, state->state_change_ids[DEACETYLATION_STATE][gene_copy], &(rates->deacetylation_num[gene_copy]), (int) 1);
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genotype, 
                      //KonStates *kon_states, 
                      TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t)
{
  float x = ran1(&seed)*((float) sum_rate_counts(rates->pic_assembly_num));

  int gene_copy; 
  int gene_loc; 
 
 LOG_ERROR("PIC assembly\n");
  /* choose a particular gene and copy to change state */
  get_gene(rates->pic_assembly_num, (int)trunc(x), &gene_loc, &gene_copy);
  int gene_id = state->state_change_ids[PICASSEMBLY_STATE][gene_copy][gene_loc];

  LOG_VERBOSE("PIC assembly event gene %d copy %d\nstate change from %d to 6\n",
              gene_id, gene_copy, state->active[gene_id][gene_copy]);

  if (state->active[gene_id][gene_copy] != ON_NO_PIC) {
    LOG_ERROR("PIC assembly event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt,
                                rates, //kon_states, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);
  
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
                           //KonStates *kon_states, 
                           TimeCourse **timecoursestart, TimeCourse **timecourselast,
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
                              //KonStates *kon_states, 
                              TimeCourse **timecoursestart, TimeCourse **timecourselast,
                              float dt, float t, float x)
{
 
  int gene_id;
  int gene_copy; 
  int gene_loc; 

  x /= TRANSCRIPTINIT;

  /* choose the gene and copy that gets transcribed */
  get_gene(rates->transcript_init_num, (int)trunc(x), &gene_loc, &gene_copy);
  gene_id = state->state_change_ids[TRANSCRIPTINIT_STATE][gene_copy][gene_loc];
  LOG_VERBOSE("transcription event gene %d, copy %d\n", gene_id, gene_copy);

  if (state->active[gene_id][gene_copy] != ON_FULL && state->active[gene_id][gene_copy] != OFF_PIC) {
    LOG_ERROR("transcription event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }
  
   LOG_ERROR("gene = %d\n", gene_id);

  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt,
                                rates, 
                                //kon_states, 
                                t, timecoursestart, timecourselast,
                                state->protein_conc);

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



/* 
 * replicate the specified gene_id: eject all TFs for that gene and replicate the DNA 
 */
void replicate_gene(CellState *state,
                    Genotype *genotype,
                    GillespieRates *rates,
                    //KonStates *kon_states,
                    //float *koffvalues,
                    int gene_id,
                    float t) 
{
  int i, k, p;// l,
  int start_tfbs_pos, end_tfbs_pos, number_tfbs_pre_replication, offset  ;//tfCount
  
  LOG("[gene %2d] duplicating at t=%g\n", gene_id, t);
  
  /* double the number of copies of gene being replicating */
  genotype->copies[gene_id] = 2*current_ploidy;


  /* remove all PICs on that gene */
  for (p=0; p < current_ploidy; p++) 
    if ((state->active[gene_id][p]==OFF_PIC || state->active[gene_id][p]==ON_FULL))
      disassemble_PIC(state, genotype, gene_id, p, rates);
  
  LOG_VERBOSE("[gene %2d] number of binding sites before adding new sites=%d at t=%g\n", 
              gene_id, genotype->binding_sites_num, t);
  
  /* do mutation: mutation rate is the mN product divided by population size */
  for (p=0; p<2*current_ploidy; p++)
    mutate(genotype, gene_id, p, mN / POP_SIZE);  

  /* record number of TFBS pre-replication */
  number_tfbs_pre_replication = genotype->sites_per_gene[gene_id];
  start_tfbs_pos = 0;
  end_tfbs_pos = 0;

  /* record the beginning and end site_ids of the pre-replication list
     of binding sites  */
  for (i=0; i<=gene_id; i++) {  
    start_tfbs_pos = end_tfbs_pos;             
    end_tfbs_pos += genotype->sites_per_gene[i];
  }
  end_tfbs_pos--;  /* always one less than the end point*/

  LOG("[gene %2d] has %d TFBS before replication [run from %d to %d]\n", 
      gene_id, number_tfbs_pre_replication, start_tfbs_pos, end_tfbs_pos);
  
  
  /* recompute *all* binding sites, then relabel sites offset by
     insertion (or deletion) of new sites created by replication */
  calc_all_binding_sites(genotype->copies, 
                         genotype->cisreg_seq, 
                         genotype->tf_seq, 
                         &(genotype->binding_sites_num),
                         &(genotype->all_binding_sites),
                         genotype->hindrance_positions,
                         genotype->sites_per_gene,
                         genotype->site_id_pos); 

  /* use new sites_per_gene and pre-replication number to compute the
     offset to shift the site_ids */
  offset = genotype->sites_per_gene[gene_id] - number_tfbs_pre_replication;

  LOG_VERBOSE("[gene %2d] has %d TFBS after replication [run from %d to %d]\n", 
              gene_id, genotype->sites_per_gene[gene_id], start_tfbs_pos, end_tfbs_pos + offset);
  LOG_VERBOSE(" shift all TFBS starting at %d by %d\n", end_tfbs_pos + 1, offset);
  LOG_VERBOSE(" number of binding sites after adding new sites=%d at t=%g\n", genotype->binding_sites_num, t);

  /* starting at the original ending point, move all sites along by
     'offset'.  note this assumes that TFBS for a particular gene are
     always stored contiguously. */
 // shift_binding_site_ids(state, kon_states, end_tfbs_pos + 1, offset);
  
  /* update the kon_states data structure to make available the newly
   * created TF binding sites in the full [start, end+offset] region
   */
  LOG("[gene %2d] adding new unbound sites from=%d to=%d\n", gene_id, start_tfbs_pos, (end_tfbs_pos + offset));

  for (k=start_tfbs_pos; k <= end_tfbs_pos + offset; k++) {
   /* add_kon(state->protein_conc[genotype->all_binding_sites[k].tf_id],
            kon_states->konvalues[genotype->all_binding_sites[k].tf_id][KON_SALPHC_INDEX],
            genotype->all_binding_sites[k].tf_id,
            k,
            rates,
            kon_states);*/
  }

  /* set acetylation state in new gene copies  */  
  for (i=0; i < current_ploidy; i++) {
    p = i + current_ploidy;
    state->state_change_ids[ACETYLATION_STATE][p][rates->acetylation_num[p]] = gene_id;
    /* update the counts for the acetylation */
    rates->acetylation_num[p]++;
    LOG_VERBOSE("[gene %2d] [clone acetylation]: ploidy=%d state_change_ids[%d][%d]=%d\n", gene_id, p, p, 
                rates->acetylation_num[p], state->state_change_ids[ACETYLATION_STATE][p][rates->acetylation_num[p]]);
    LOG_VERBOSE(" rates->acetylation_num[%d]=%d\n", p, rates->acetylation_num[p]);
  }
}



 //DESTROY
 //CHANGE JU
void recalibrate_cell(GillespieRates *rates,
                      CellState *state,
                      Genotype *genotype,
                      //KonStates *kon_states,
                      //float **koffvalues,
                      float mRNAdecay[NGENES],
                      float transport[NGENES],
                      float dt) // float *active_gene
{
  int i, j; 
  
  //float protein_decay;
 // float salphc = 0.0;

  // reset the total rate for current step 
  rates->subtotal=0.0;
  
  // reset all rates and operations 
  //rates->koff=0.0;
  //rates->koff_operations=0;

  rates->transport=0.0;
  rates->transport_operations=0;

  rates->mRNAdecay=0.0;
  rates->mRNAdecay_operations=0;

  rates->pic_disassembly=0.0;
  rates->pic_disassembly_operations=0;

  

  for (i=0; i < NGENES; i++) {
    // transport rates 
    transport[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
    rates->transport += transport[i];
    rates->transport_operations++;

    // regenerate decay rates 
    mRNAdecay[i] = genotype->mRNAdecay[i] * ((float) state->mRNA_cyto_num[i] + (float) state->mRNA_transl_cyto_num[i]);
    rates->mRNAdecay += mRNAdecay[i];
    rates->mRNAdecay_operations++;

  }
  LOG_ERROR("mRNAdecay = %f\n", rates->mRNAdecay);
  /* recompute koffvalues for all sites: note that a time is required
     by the function only for diagnostic purposes, it isn't used by
     the function, and since we don't have available a current time in
     this function, we simply feed it a dummy time of 0.0  
  recompute_koff_rates(rates, state, genotype, *koffvalues, 0.0);*/

  // recompute and cache the total rate in data structure 
  rates->subtotal += rates->transport;
  rates->subtotal += rates->mRNAdecay;
  rates->subtotal += rates->pic_disassembly;

   
    //convert the counts back into rates using the constants 
  LOG_ERROR("ratesTWO = %d\n", rates->acetylation_num[0]);
  for (j=0; j < MAX_COPIES; j++) {
    rates->subtotal += (float) rates->acetylation_num[j] * ACETYLATE;
    rates->subtotal += (float) rates->deacetylation_num[j] * DEACETYLATE;
    rates->subtotal += (float) rates->pic_assembly_num[j] * PICASSEMBLY;
    rates->subtotal += (float) rates->transcript_init_num[j] * TRANSCRIPTINIT;    
  } 
}

/*
 * implement actual movement of genetic material: move 'from_copy' of
 * 'gene' in 'from_genotype' to 'to_copy' in 'to_genotype' and
 * reconstruct all the lists of bound TFs and associated hindrances
 */
/*int move_gene_copy(int from_copy,
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

  // shift all sites in tf_bound_indexes 
  for (k = 0; k < from_state->tf_bound_num; k++) {
    site_id = from_state->tf_bound_indexes[k];
    gene_id = from_genotype->all_binding_sites[site_id].cisreg_id;
    gene_copy = from_genotype->all_binding_sites[site_id].gene_copy;

    // check to see if TF is within the gene copy to be moved 
    if (site_id >= start && site_id <= end && gene_copy == from_copy && gene_id == gene) {
      // fix  
      if (start  == lastpos + 1)  // if no offset required: keep site_id 
        new_site_id = site_id;
      else                        // otherwise offset them by the difference 
        new_site_id = site_id  - (start - lastpos);
      LOG_VERBOSE_NOCELLID("gene=%d [copy=%d] orig site_id=%d, move to new [copy=%d] site_id=%d (new tf_bound_num=%d), lastpos=%d, offset=%d\n", 
                           gene, from_copy, site_id, to_copy, new_site_id, to_state->tf_bound  _num, lastpos, (start - lastpos));
      to_state->tf_bound_indexes[to_state->tf_bound_num] = new_site_id;
      to_state->tf_bound_num++;
    }
  }
//DESTROY
  // shift TF hindered indexes 
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

  lastpos += (end - start) + 1;   // update the last position in new binding site array 
  LOG_VERBOSE_NOCELLID("gene=%d [copy=%d] start=%d, end=%d, lastpos=%d\n", gene, from_copy, start, end, lastpos);
  return lastpos;
}*/

/*
 * clone a queue of pending events by adding them to the new
 * FixedEvent list and removing them from the old
 */
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

/*
 * split the mRNA and associated events between mother and daughter
 */
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
  /* regenerate initial queue(s) */
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

/*
 * make a complete 'read-only' copy of the mother's cell state.
 */
void clone_cell(Genotype *genotype_orig,
                CellState *state_orig,
                GillespieRates *rates_orig,
                Genotype *genotype_clone,
                CellState *state_clone,
                GillespieRates *rates_clone)
{
  int i, p;// k,

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
    genotype_clone->replication_time[i] =  genotype_orig->replication_time[i];
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

  /* start CellState copy */

  

  // copy activation state 
  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      state_clone->active[i][p] = state_orig->active[i][p];
      state_clone->state_change_ids[ACETYLATION_STATE][p][i] = state_orig->state_change_ids[ACETYLATION_STATE][p][i];
      state_clone->state_change_ids[DEACETYLATION_STATE][p][i] = state_orig->state_change_ids[DEACETYLATION_STATE][p][i];
      state_clone->state_change_ids[PICASSEMBLY_STATE][p][i] = state_orig->state_change_ids[PICASSEMBLY_STATE][p][i];
      state_clone->state_change_ids[TRANSCRIPTINIT_STATE][p][i] = state_orig->state_change_ids[TRANSCRIPTINIT_STATE][p][i];
      state_clone->state_change_ids[PICDISASSEMBLY_STATE][p][i] = state_orig->state_change_ids[PICDISASSEMBLY_STATE][p][i];

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

/*
 * initialize the other (non-DNA) parts of the genotype for a newly
 * divided cell based on a copy 'genotype_clone'
 */
 //DESTROY
void initialize_new_cell_genotype(Genotype *genotype, Genotype *genotype_clone)
{
  int i, p;

  /* initialize hindrance for all TFGENES */
  for (p=0; p < TFGENES; p++) {
    genotype->hindrance_positions[p] = genotype_clone->hindrance_positions[p];
  }

  /* initialize the non-cisregulatory parts of the genotype */
  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      genotype->activating[i][p]=  genotype_clone->activating[i][p];
      genotype->pic_disassembly[i][p]=  genotype_clone->pic_disassembly[i][p];
    }
    genotype->replication_time[i] =  genotype_clone->replication_time[i];
    genotype->mRNAdecay[i] =  genotype_clone->mRNAdecay[i];
    genotype->proteindecay[i] =  genotype_clone->proteindecay[i];
    genotype->translation[i] =  genotype_clone->translation[i];
  }
}

/*
 * initialize a newly divided cell's CellState based on copy 'state_clone'
 */
 //DESTROY
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
  //state->tf_bound_num = 0;
 // state->tf_hindered_num = 0;

  state->in_s_phase = 0;   /* reset S phase state to 0 */
  state->burn_in = 0;      /* only do burn-in at beginning of runs */
  state->division_time = TIME_INFINITY;  /* reset division time for this cell */
  state->divisions = state_clone.divisions; /* copy division counter */
  state->cell_size = state_clone.cell_size*fraction;   /* reset cell size */

  /* growth rate, take instantaneous growth rate of mother just before
     division, this will be updated after first new time step */
  state->growth_rate = state_clone.growth_rate;

  /* keep thermodynamic state the same */ 
  state->RTlnKr = state_clone.RTlnKr;
  state->temperature = state_clone.temperature;

  /* initialize the rate counts */
  for (i=0; i < MAX_COPIES; i++) {
    for (j=0; j < NGENES; j++)
      state->state_change_ids[ACETYLATION_STATE][i][j] = -1;
    rates->acetylation_num[i]=0;
    rates->deacetylation_num[i]=0;
    rates->pic_assembly_num[i]=0;
    rates->transcript_init_num[i]=0;
    rates->pic_disassembly_num[i]=0;
  }
}

/*
 * initialize the newly divided cells state change ID state for
 * particular gene and gene copy
 */
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

/*
 * initialize the newly divided cell's DNA genotype: cisreg and TF
 * sequences for a specified 'gene_id' and 'copy1' and 'copy2' pair
 * from the original 'genotype_clone'
 */
void initialize_new_cell_gene(Genotype *genotype, Genotype genotype_clone, 
                              CellState *state, CellState *state_clone,
                              GillespieRates *rates, GillespieRates *rates_clone,
                              int gene_id, int copy1, int copy2, int *lastpos)
{
  int j, k;

  for (k=0; k < CISREG_LEN; k++) {
    genotype->cisreg_seq[gene_id][0][k] = genotype_clone.cisreg_seq[gene_id][copy1][k];
    genotype->cisreg_seq[gene_id][1][k] = genotype_clone.cisreg_seq[gene_id][copy2][k];
    genotype->cisreg_seq[gene_id][2][k] = genotype_clone.cisreg_seq[gene_id][copy1][k];
    genotype->cisreg_seq[gene_id][3][k] = genotype_clone.cisreg_seq[gene_id][copy2][k];
  }

  /* don't update TF if this gene_id is not controlling a TF */
  if (gene_id < TFGENES) {
    for (k=0; k < TF_ELEMENT_LEN; k++) {
      genotype->tf_seq[gene_id][0][k] = genotype_clone.tf_seq[gene_id][copy1][k];
      genotype->tf_seq[gene_id][1][k] = genotype_clone.tf_seq[gene_id][copy2][k];
      genotype->tf_seq[gene_id][2][k] = genotype_clone.tf_seq[gene_id][copy1][k];
      genotype->tf_seq[gene_id][3][k] = genotype_clone.tf_seq[gene_id][copy2][k];
    }
  }


  if (genotype_clone.copies[gene_id] - 1 >= copy1) {
    //*lastpos = move_gene_copy(copy1, 0, gene_id, &genotype_clone, genotype, state_clone, state, *lastpos);
  }

  if (genotype_clone.copies[gene_id] - 1 >= copy2) {
    //*lastpos = move_gene_copy(copy2, 1, gene_id, &genotype_clone, genotype, state_clone, state, *lastpos);
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

  /* ACETYLATION_STATE */
  initialize_new_cell_state_change_ids(state, *state_clone, ACETYLATION_STATE,
                                     rates->acetylation_num, rates_clone->acetylation_num,
                                     copy1, copy2, gene_id);
  /* DEACETYLATION_STATE */
  initialize_new_cell_state_change_ids(state, *state_clone, DEACETYLATION_STATE,
                                     rates->deacetylation_num, rates_clone->deacetylation_num,
                                     copy1, copy2, gene_id);
  /* PICASSEMBLY_STATE */
  initialize_new_cell_state_change_ids(state, *state_clone, PICASSEMBLY_STATE,
                                     rates->pic_assembly_num, rates_clone->pic_assembly_num,
                                     copy1, copy2, gene_id);
  /* TRANSCRIPTINIT_STATE */
  initialize_new_cell_state_change_ids(state, *state_clone, TRANSCRIPTINIT_STATE,
                                     rates->transcript_init_num, rates_clone->transcript_init_num,
                                     copy1, copy2, gene_id);
  /* PICDISASSEMBLY_STATE */
  initialize_new_cell_state_change_ids(state, *state_clone, PICDISASSEMBLY_STATE,
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

  while (timeEnd != NULL) {   /* loop through old queue from 'state_clone' */
    for (j=0; j < MAX_COPIES; j++)
      LOG_VERBOSE("mRNA_transcr_num[%2d][%d]=%d\n", 
                  gene_id, j, state_clone->mRNA_transcr_num[gene_id][j]);
    
    int gene_id_queue = timeEnd->gene_id;
    int copy = timeEnd->copy;
    float time = timeEnd->time;
    
    if (gene_id_queue == gene_id && copy == copy1) { 
      /* move to correct new copy */
      state->mRNA_transcr_num[gene_id][0]++;
      /* add event for copy 0 to new queue */
      add_fixed_event(gene_id, 0, time, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end));
      printf("move mRNATranscrTime event at time=%g on gene %d copy=%d to copy=0 of total mRNA_transcr_num=%d\n", 
             time, gene_id, copy1, state_clone->mRNA_transcr_num[gene_id][copy1]);
      /* remove from original queue */
      LOG("removing mRNATranscrTime event at time=%g on gene %2d (copy %d) from clone of queue\n", time, gene_id, copy1);
      delete_fixed_event(gene_id, copy1, 0, &(state_clone->mRNA_transcr_time_end), &(state_clone->mRNA_transcr_time_end_last));
      (state_clone->mRNA_transcr_num[gene_id][copy1])--;
    } else if (gene_id_queue == gene_id && copy == copy2) {  
      /* move to correct new copy */
      state->mRNA_transcr_num[gene_id][1]++;
      /* add event for copy 1 to new queue */
      add_fixed_event(gene_id, 1, time, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end));
      printf("move mRNATranscrTime event at time=%g on gene %d copy=%d to copy=1 of total mRNA_transcr_num=%d\n", 
             time, gene_id, copy2, state_clone->mRNA_transcr_num[gene_id][copy2]);
      /* remove from original queue */
      LOG("removing mRNATranscrTime event at time=%g on gene %2d (copy %d) from clone of queue\n", time, gene_id, copy2);
      delete_fixed_event(gene_id, copy2, 0, &(state_clone->mRNA_transcr_time_end), &(state_clone->mRNA_transcr_time_end_last));
      (state_clone->mRNA_transcr_num[gene_id][copy2])--;
    }
    timeEnd = timeEnd->next; /* get next time in old queue */
    fflush(stdout);
  }
}


/*
 * do the actual cell division for the specified mother_id and
 * daughter_id
 */
 //DESTROY
void do_cell_division(int mother_id,
                      int daughter_id,
                      int keep_mother, 
                      int keep_daughter,
                      Genotype *mother,
                      CellState *mother_state,
                      GillespieRates *mother_rates,
                      //KonStates *mother_kon_states,
                      //float **mother_koffvalues,
                      float mother_transport[NGENES],
                      float mother_mRNAdecay[NGENES],
                      
                      Genotype *daughter,
                      CellState *daughter_state,
                      GillespieRates *daughter_rates,
                      //KonStates *daughter_kon_states,
                      //float **daughter_koffvalues,
                      float daughter_mRNAdecay[NGENES],
                      float daughter_transport[NGENES],
                      float fraction,
                      float x,
                      float dt)
{
  int i, j, total;

  /* clone of cell */
  Genotype genotype_clone;
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

  /* create a "read-only" clone of the mother cell */
  clone_cell(mother, mother_state, mother_rates,
             &genotype_clone, &state_clone, &rates_clone);
  
  /* free the existing memory for the time queue */
  if (keep_daughter) delete_queues(daughter_state);
  if (keep_mother) delete_queues(mother_state);

  if (keep_daughter) {
    /* initialize the non-cis-regulatory part of the genotype 
       of the new daughter cell based on clone of mother */
    initialize_new_cell_genotype(daughter, &genotype_clone);
    /* initialize the state of the new cell, including size */
    initialize_new_cell_state(daughter_state, state_clone, daughter_rates, fraction);
    /* set the founder_id of the daughter cell from the mother cell */
    daughter_state->founder_id = state_clone.founder_id;
    daughter_state->divisions = 0;    /* daughter cell resets divisions */ 
    printf("daughter cell %03d is founded by cell %03d with %2d divisions\n", 
           daughter_id, daughter_state->founder_id, daughter_state->divisions);
    LOG_NOCELLID("daughter cell %03d is founded by cell %03d with %2d divisions\n", 
                 daughter_id, daughter_state->founder_id, daughter_state->divisions);
  }
    
  if (keep_mother) {
    initialize_new_cell_genotype(mother, &genotype_clone);
    initialize_new_cell_state(mother_state, state_clone, mother_rates, (1.0-fraction));
    mother_state->divisions++;      /* update divisions in mother */
    printf("mother   cell %03d is founded by cell %03d with %2d divisions\n", 
           mother_id, mother_state->founder_id, mother_state->divisions);
    LOG_NOCELLID("mother   cell %03d is founded by cell %03d with %2d divisions\n", 
                 mother_id, mother_state->founder_id, mother_state->divisions);
  }

 // LOG_VERBOSE_NOCELLID("[cell %03d] (mother) total tf_bound_num=%d, tf_hindered_num=%d\n", mother_id, 
                      // state_clone.tf_bound_num, state_clone.tf_hindered_num);

  //print_all_binding_sites(mother->copies, mother->all_binding_sites, mother->binding_sites_num, 
  //                        mother->tf_seq, mother->cisreg_seq, mother->site_id_pos); 

  if (verbose) {
    LOG_NOCELLID("[cell %03d]: acetylation counts:\n", mother_id);
    total = 0;
    for (j=0; j < MAX_COPIES; j++)  {
      LOG_NOCELLID(" before mother acetylation_num[%2d]=%2d\n", j, rates_clone.acetylation_num[j]);
      total += rates_clone.acetylation_num[j];
    }
    LOG_NOCELLID(" before mother total acetylation=%d\n", total);
  }

  /* now split up the genes into one or other of the new cells */
  /* this assumes division of a diploid cell */
  for (i=0; i < NGENES; i++) {

    /* reset the number of copies of gene in mother and daughter after
       division to original ploidy */
    if (keep_daughter)
      daughter->copies[i] = current_ploidy;
    if (keep_mother)
      mother->copies[i] = current_ploidy;

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

    if (keep_daughter) {
      LOG_NOCELLID("initialize daughter=%2d, gene_id=%2d\n", daughter_id, i);
      initialize_new_cell_gene(daughter, genotype_clone, 
                               daughter_state, &state_clone,
                               daughter_rates, &rates_clone,
                               i, daughter_copy1, daughter_copy2, &lastpos_daughter);
    }

    if (keep_mother) {
      LOG_NOCELLID("initialize mother  =%2d, gene_id=%2d\n", mother_id, i);
      initialize_new_cell_gene(mother, genotype_clone, 
                               mother_state, &state_clone,
                               mother_rates, &rates_clone,
                               i, mother_copy1, mother_copy2, &lastpos_mother);
    }
  }

  if (verbose) {    /* debugging output only */
    total = 0;
    for (j=0; j < MAX_COPIES; j++)  {
      LOG_NOCELLID(" after daughter acetylation_num[%2d]=%2d\n", j, daughter_rates->acetylation_num[j]);
      LOG_NOCELLID(" after mother acetylation_num[%2d]=%2d\n", j, mother_rates->acetylation_num[j]);
      total += daughter_rates->acetylation_num[j];
      total += mother_rates->acetylation_num[j];
    }
    LOG_NOCELLID(" after total acetylation=%d\n", total);

    /* check that we have emptied the list of transcribing mRNAs */
    total = 0;
    for (i=0; i < NGENES; i++) 
      for (j=0; j < MAX_COPIES; j++) 
        total += state_clone.mRNA_transcr_num[i][j];
    
    LOG_NOCELLID("mRNA_transcr_num=%2d left in state_clone after moving to mother+daughter\n", total);
  }

  /*if (keep_daughter) {
    realloc_cell_memory(daughter_state, daughter_koffvalues);
  }
  if (keep_mother) {
    realloc_cell_memory(mother_state, mother_koffvalues);
  }*/

  if (keep_daughter) {
    // recompute *all* binding sites in daughter, then relabel sites */
    calc_all_binding_sites(daughter->copies, 
                           daughter->cisreg_seq, 
                           daughter->tf_seq, 
                           &(daughter->binding_sites_num),
                           &(daughter->all_binding_sites),
                           daughter->hindrance_positions,
                           daughter->sites_per_gene,
                           daughter->site_id_pos);
  }
    
  if (keep_mother) {
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

  if (keep_daughter && keep_mother) {
    /* check that number of binding sites in mother before division equals sum when of mother and daughter after */
    LOG_NOCELLID("original number of binding sites=%d should = (mother=%d + daughter=%d) = %d\n", 
                 original_bind_count, mother->binding_sites_num, daughter->binding_sites_num, mother->binding_sites_num + daughter->binding_sites_num);
    if (original_bind_count != mother->binding_sites_num + daughter->binding_sites_num) {
      LOG_ERROR_NOCELLID("original number of binding sites=%d  != (mother=%d + daughter=%d) = %d\n", 
                         original_bind_count, mother->binding_sites_num, daughter->binding_sites_num, mother->binding_sites_num + daughter->binding_sites_num);
      exit(0);
    }
  }
  
  /* split up the volume of the cell */
  /* do protein first */
  for (i=0; i < NPROTEINS; i++) {
    if (keep_daughter) {
      daughter_state->protein_conc[i] = fraction * state_clone.protein_conc[i];
    }
    if (keep_mother) {
      mother_state->protein_conc[i] = (1-fraction) * state_clone.protein_conc[i];
      LOG_VERBOSE_NOCELLID("daughter=%g (%g), mother=%g (%g) = total protein=%g\n", 
                   daughter_state->protein_conc[i], fraction, mother_state->protein_conc[i], (1-fraction), state_clone.protein_conc[i]);
    }
  }

  /* now loop over all transcribing genes */
  for (i=0; i < NGENES; i++) {
    /* mRNAs in cytoplasm (not translating) */
    int mRNA_cyto_to_daughter = rint(fraction * state_clone.mRNA_cyto_num[i]);
    if (keep_daughter) {
      daughter_state->mRNA_cyto_num[i] = mRNA_cyto_to_daughter;
      LOG_VERBOSE_NOCELLID("daughter=%d ", daughter_state->mRNA_cyto_num[i]);
    }
    if (keep_mother) {
      mother_state->mRNA_cyto_num[i] =  state_clone.mRNA_cyto_num[i] - mRNA_cyto_to_daughter;
      LOG_VERBOSE_NOCELLID("mother=%d ", mother_state->mRNA_cyto_num[i]);
    }
    LOG_VERBOSE_NOCELLID("total mRNA_cyto_num=%d\n", state_clone.mRNA_cyto_num[i]);

    /* mRNAs in nucleus */
    int mRNA_nuclear_to_daughter = rint(fraction * mother_state->mRNA_nuclear_num[i]);
    if (keep_daughter) {
      daughter_state->mRNA_nuclear_num[i] = mRNA_nuclear_to_daughter;
      LOG_VERBOSE_NOCELLID("daughter=%d ", daughter_state->mRNA_nuclear_num[i]);
    }
    if (keep_mother) {
      mother_state->mRNA_nuclear_num[i] =  state_clone.mRNA_nuclear_num[i] - mRNA_nuclear_to_daughter;
      LOG_VERBOSE_NOCELLID("mother=%d ", mother_state->mRNA_nuclear_num[i]);
    }
    LOG_VERBOSE_NOCELLID("of total mRNA_nuclear_num=%d\n", state_clone.mRNA_nuclear_num[i]);

    // TODO: fix to take into account mother replacing daughter or vice versa
    /* split up currently translating mRNAs (mRNA_transl_cyto_num) along with associated FixedTime events */
    split_mRNA(&(state_clone.mRNA_transl_time_end), &(state_clone.mRNA_transl_time_end_last), state_clone.mRNA_transl_cyto_num,
               &(daughter_state->mRNA_transl_time_end), &(daughter_state->mRNA_transl_time_end_last), daughter_state->mRNA_transl_cyto_num,
               &(mother_state->mRNA_transl_time_end), &(mother_state->mRNA_transl_time_end_last), mother_state->mRNA_transl_cyto_num,
               i, fraction);
  }  
  fflush(stdout);

  if (keep_daughter) {
    /* recompute rates in daughter */
    recalibrate_cell(daughter_rates, daughter_state, daughter,
                    // daughter_kon_states, daughter_koffvalues,
                     daughter_mRNAdecay, daughter_transport, dt);
   // LOG_VERBOSE_NOCELLID("tf_bound_num=%d (daughter_id=%d)\n", daughter_state->tf_bound_num, daughter_id);
    //LOG_VERBOSE_NOCELLID("tf_hindered_num=%d (daughter_id=%d)\n", daughter_state->tf_hindered_num, daughter_id);
  }

  if (keep_mother) {
    /* recompute rates in mother */
    recalibrate_cell(mother_rates, mother_state, mother,
                    // mother_kon_states, mother_koffvalues,
                     mother_mRNAdecay, mother_transport, dt);
    //LOG_VERBOSE_NOCELLID("tf_bound_num=%d (mother_id=%d)\n", mother_state->tf_bound_num, mother_id);
    //LOG_VERBOSE_NOCELLID("tf_hindered_num=%d (mother_id=%d)\n", mother_state->tf_hindered_num, mother_id);
  }

  /* free the memory associated with temporary copy of mother */
  //free_mem_CellState(&state_clone);
}

/*
 * diagnostic function to dump a copy of the GillespieRates and some
 * of the CellState to the current error file
 */
 //DESTROY?
void log_snapshot(GillespieRates *rates,
                  CellState *state,
                  Genotype *genotype,
                  //KonStates *kon_states,
                  //float **koffvalues,
                  float mRNAdecay[NGENES],
                  float transport[NGENES],
                  //float konrate,
                  float x,
                  float t)
{
  int  p, nkon = 0;//i,

  
  for (p=0; p < MAX_COPIES; p++) {
    LOG_NOFUNC(" acetylation=%g (copy %d)\n deacetylation=%g (copy %d)\n PIC assembly=%g (copy %d)\n transcriptinit=%g (copy %d)\n",
               (float)rates->acetylation_num[p]*ACETYLATE, p, (float)rates->deacetylation_num[p]*DEACETYLATE, p, 
               (float)rates->pic_assembly_num[p]*PICASSEMBLY, p, (float)rates->transcript_init_num[p]*TRANSCRIPTINIT, p);
  }
  
  LOG_NOFUNC(" nkon recomputed=%d\n", nkon);
  LOG_NOFUNC("\n");
}

//ADD JU
int intcmp(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

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
  float *prob;
   /*printf("%f\n ", initProteinConc[m]);
       for (m=0; m< NGENES; m++){
           printf("%f ", initProteinConc[m]);
       }*/
       //system("PAUSE");
       qsort((void *) &(indiv.all_binding_sites[start]), indiv.sites_per_gene[gene],                                 
            sizeof(struct AllTFBindingSites),(compfn)compare );
   
     // recordFile = fopen("recordFile.txt", "w");
     
     // if ((recordFile = fopen("recordFile.txt", "w"))) {
     
     // fprintf(recordFile, "BS num  = %d\n", indiv.binding_sites_num);
     // fprintf(recordFile,"BS num gene 1: %d\n", indiv.sites_per_gene[1]);
     /* for(n=0; n<indiv.sites_per_gene[gene]; n++){
               LOG_ERROR_NOCELLID("%d   %d \n",n,  indiv.all_binding_sites[start+n].left_edge_pos);
              // fprintf(recordFile, "%d   %d   %d\n",n,  indiv.all_binding_sites[start+n].left_edge_pos,
                                   // indiv.activating[7][indiv.all_binding_sites[0].gene_copy]);
      }*/
      
      struct Dtype *arrayD;
 arrayD = malloc(NITER*sizeof(struct Dtype));
  
 /*each struct holds ID, BS numer, start positon, hamming distance, concentration and weight*/
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
 
         //fprintf(recordFile, "%d  LEP = %d  Hd = %d   tf = %d Kon = %f weight = %.2f\n",arrayWT[lem].tfbsNum, arrayWT[lem].startPos, arrayWT[lem].hammDist, arrayWT[lem].tfIDon, arrayWT[lem].conc, arrayWT[lem].weight);
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
               /*FIX THIS! SUM IS NOT CORRECT! ROUNDING ERRORS!*/    
      }
      weightSum += arrayWT[0].weight;
      
      
      prob = malloc((TFBS+1)*sizeof(float));
    
      for(n =0; n<TFBS+1; n++){
         if(n == TFBS) { prob[n] = (arrayWT[0].weight) / weightSum;}
         else {prob[n] = (arrayWT[n].weight) / weightSum;}
        // fprintf(recordFile, "prob[%d] = %f\n", n, prob[n]);
      }
      m = rand()%1000;
      check =0.;
      check = m/1000.;
      
      //LOG_ERROR_NOCELLID("m = %f\n", m);
      //printf("check=%f\n", check);
      partition = malloc((TFBS+1)*sizeof(float));
      for(n=0; n<TFBS+1; n++){
         if(n==0){partition[n] = prob[n];}
         else if(n!=TFBS){ partition[n] = (partition[(n-1)]+prob[n]);}
         else {partition[n] = 1;}
        // fprintf(recordFile, "partition[%d] = %f\n", n, partition[n]);
      }
      b = 0;
      n =0;
     checkP =check;
    // fprintf(recordFile, "checkP = %f\n", checkP);
      
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
      if(arrayWT[b].tfIDon == 11){/*do nothing, there is nothing bound */
         unboundCount++;}
      else{
           if(ran1(&seed)<PROB_ACTIVATING){A++;}
       // if(indiv.activating[arrayWT[b].tfIDon][indiv.all_binding_sites[(arrayWT[b].tfbsNum)].gene_copy] ==1){
            // LOG_ERROR_NOCELLID("ACTIVATING!!!\n");
            // A++;
       // }
        else{R++;}
      }
      if(A+R > 9){
            // fprintf(recordFile, "PROBLEM!!!!!!!");
             system("PAUSE");
             //TO DO: fix this problem!! Too many things get bound, not sure what is wrong!
      }
     // fprintf(recordFile, "A=%d, R=%d\n", A, R);
     //LOG_ERROR_NOCELLID("A=%d, R=%d\n", A, R);  
      
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
       arrayD[size].ratio = (float)arrayD[size].active * .33 + .31;
       arrayD[size].count = 1;
       size++;
    }
    val++;
      }//closes while(val-NITER)
      
      activeCount=0;
      for(lem=0; lem<size; lem++){
                 //LOG_ERROR_NOCELLID("%d  count = %d  active = %.2f, repress = %.2f, ratio = %.3f\n", lem, arrayD[lem].count, arrayD[lem].active, arrayD[lem].repress, arrayD[lem].ratio);
                // printf("%d  count = %d  active = %.2f, repress = %.2f, ratio = %.3f\n", lem, arrayD[lem].count, arrayD[lem].active, arrayD[lem].repress, arrayD[lem].ratio);
      if((float)arrayD[lem].repress <= (float)arrayD[lem].ratio){
         activeCount += arrayD[lem].count;} 
         //printf("active_count = %d\n", activeCount);
         //LOG_ERROR_NOCELLID("active Count = %d\n", activeCount);
      }
     // printf("%d\n", activeCount);
      //system("PAUSE");
      percent = activeCount/(float)NITER;
      //LOG_ERROR_NOCELLID("active count = %d, NITER = %f\n", activeCount, NITER);
     // fprintf(recordFile, "unboundCount = %d\n", unboundCount);
     // fprintf(recordFile, "activeCount = %d\n", activeCount);
     // fprintf(recordFile, "percent = %f\n", percent);
      
      free(arrayD);
      free(partition);
      free(arrayWT); 
      free(prob); 
      free(left_edge_pos);
      //}//closes recordFile loop
      
     // fclose(recordFile);
      LOG_ERROR_NOCELLID("percent = %f\n", (float)percent);
      return percent;
 }
 
 void active_vect(Genotype indiv, float initProteinConc[NGENES], float *gene_active){
       int geneSum = 0;
       int i;
     // float prob = active_to_repress(indiv, initProteinConc, 0, 0);
      // printf("\n prob = %f\n", prob);
     
       for(i=0; i<NGENES; i++){ 
          // printf("geneSum = %d,  i=%d\n", geneSum, i);
           gene_active[i] = active_to_repress(indiv, initProteinConc, geneSum, i);
           LOG_ERROR_NOCELLID("%f\n", gene_active[i]);
          // printf("geneSum = %d, i = %d, WHAT = %f\n\n", geneSum, i, gene_active[i] );
           //fprintf(recordFile, "geneSum = %d, i = %d, WHAT = %f\n\n", geneSum, i, gene_active[i] );
           //LOG_ERROR_NOCELLID("geneSum = %d, i = %d, WHAT = %f\n", geneSum, i, gene_active[i] );
           geneSum += indiv.sites_per_gene[i];
        
  }
 }  

//END JU


/*
 * run the model for a specified cell for a single timestep:
 *  - returns 0  if cell is not "dead" (i.e. rates haven't deteroriated to zero)
 *  - returns -1 if cell is "dead"
 */
 //DESTROY
int do_single_timestep(Genotype *genotype, 
                       CellState *state, 
                       //KonStates *kon_states, 
                       GillespieRates *rates, 
                       float *t,
                       //float *koffvalues,
                       float transport[NGENES],
                       float mRNAdecay[NGENES],
                       float *x,
                       float *dt,
                      // float *konrate,
                       TimeCourse *timecoursestart[NPROTEINS],
                       TimeCourse *timecourselast[NPROTEINS],
                       int maxbound2,
                       int maxbound3,
                       int no_fixed_dev_time,
                       float *genesActive) 
{
  int i, j;
  int event;     /* boolean to keep track of whether FixedEvent has ended */
  int total;     /* total possible translation events */
  float fixed_time;
  clock_t start, stop;
  double tt = 0.0;
  
  if (1) //VERBOSE JU CHANGE
    for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
       //LOG_ERROR("rates->deacetylation_num[%d]=%d\n", j, rates->deacetylation_num[j]);
    }
    for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      //LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
       LOG_ERROR("rates->deacetylation_num[%d]=%d\n", j, rates->deacetylation_num[j]);
    }
  //LOG_ERROR("ratesTWO = %d\n", rates->acetylation_num[2]);
  
  /* if burn-in is set and we have completed the 0.003 sec (= 0.00005
     min) burn in, then recompute the rates for the new burn-in.  The
     0.003 sec length was found by multiple runs of the model with
     kon=0.225 and finding the length of time it took for the number
     of bound TFs to equilibrate  */
    // printf("INSIDE SINGLE TIME STEP!\n");
  if (*t > 0.00005 && state->burn_in) {
    printf("recalibrating cell %3d after burn-in!\n", state->cell_id);
    LOG("recalibrating cell %3d after burn-in!\n", state->cell_id);
    /*log_snapshot(rates, state, genotype, kon_states, &koffvalues,
                 mRNAdecay, transport, *konrate, *x, *t);*/

    kon = kon_after_burnin;     /* reset kon to the post-burn-in value */
    recalibrate_cell(rates,
                     state,
                     genotype,
                    // kon_states,
                     //&koffvalues,
                     mRNAdecay,
                     transport,
                     *dt); 
    
    active_vect(genotype[0], state->protein_conc, genesActive);
    //printf("STOP!");
    //system("PAUSE");
    for(j=0; j<NGENES; j++){
             //revise_prob_vect(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
    revise_activity_state(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
    }
    
    /*log_snapshot(rates, state, genotype, kon_states, &koffvalues,
                 mRNAdecay, transport, *konrate, *x, *t);*/
    state->burn_in = 0;    /* now disable burn-in */
  } 

  /* compute S-phase offsets */
  if (critical_size > 0.0 && state->cell_size >= critical_size && !state->in_s_phase)  { /* run until checkpoint reached */
    reach_s_phase(state, genotype, *t);
    /* current time plus 30 mins for S phase and a further 30 mins of growth after S phase */
    state->division_time = *t + time_s_phase + time_g2_phase;   
    LOG("at t=%g add division time=%g, cell_size=%g\n", *t, state->division_time, state->cell_size);
  }
  
#if 0 /* currently disable printing out the TF occupancy and amount of rounding */
  print_tf_occupancy(state, genotype->all_binding_sites, *t);
  print_rounding(state, rates, *t);
#endif
  
  *x = expdev(&seed);        /* draw random number */
 // LOG_ERROR("random x = %f\n", *x);
  
    
  *dt = .03 * (*x);
  /* compute the initial dt for the next event 
  calc_dt(x, dt, rates, //kon_states,
   mRNAdecay, genotype->mRNAdecay,
          state->mRNA_cyto_num, state->mRNA_transl_cyto_num, state->cell_id);*/

  

  if (*dt < 0.0) {
    LOG_ERROR("dt=%g is negative after calc_dt, t=%g\n", *dt, *t);
    exit(-1);
  }

  LOG_VERBOSE("next stochastic event due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);
 
  if((start = clock()) != -1){
  active_vect(genotype[0], state->protein_conc, genesActive);
  stop = clock();
  tt = (double)(stop-start)/CLOCKS_PER_SEC;
  //printf("Run Time: %f\n", tt);
  //system("PAUSE"); 
 }
 

  
  LOG_VERBOSE("ACTIVE VECT\n");
    //for(j=0; j<NGENES; j++){
            // printf("%f\n", genesActive[j]);
             //system("PAUSE");
             //revise_prob_vect(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
     //revise_activity_state(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
   // }
    
      
  /* first check to see if a fixed event occurs in current t->dt window,
     or in tdevelopment if running for a fixed development time */
  fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);

  event = does_fixed_event_end(state->mRNA_transl_time_end,
                               state->mRNA_transcr_time_end,
                               state->replication_time_end,
                               fixed_time);
  //printf("rates transportHere = %f\n", rates->transport);
  /* while there are either transcription or translation events
     occuring in current t->dt window */
  //LOG_ERROR("Event = %d\n", event);
  while (event > 0) {
       // printf("event > 0\n");
    //*konrate = (*x)/(*dt);
    
    switch (event) {
    case 1:   /* if a transcription event ends */
    //printf("case 1");
      end_transcription(dt, *t, state, transport, rates);
      
      update_protein_conc_cell_size(state->protein_conc, state, genotype, *dt,
                                    rates, //kon_states,
                                     *t,
                                    timecoursestart, timecourselast,
                                    state->protein_conc);
      /*active_vect(genotype[0], state->protein_conc, genesActive);
      for(j=0; j<NGENES; j++){
               revise_activity_state(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
       }*/
      break;
    case 2:            /* if a translation event ends */
    //printf("case 2");
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
                                    rates, //kon_states, 
                                    *t,
                                    timecoursestart, timecourselast,
                                    state->protein_conc);
      /*active_vect(genotype[0], state->protein_conc, genesActive);
      for(j=0; j<NGENES; j++){
               revise_activity_state(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
       }*/
      /* the number of mRNAs in cytoplasm affects binding */
      //change_mRNA_cytoplasm(i, genotype, state, rates );//kon_states
      break;
    case 3:  /* replicate gene */
    //printf("case 3");
      *dt = state->replication_time_end->time - *t;         /* make dt window smaller */
      
      replicate_gene(state, genotype, rates, state->replication_time_end->gene_id, *t);// kon_states, koffvalues, 
      
      /* delete the event that just happened */
      delete_fixed_event_start(&(state->replication_time_end), &(state->replication_time_end_last));
      
      update_protein_conc_cell_size(state->protein_conc, state, genotype, *dt,
                                    rates, //kon_states, 
                                    *t,
                                    timecoursestart, timecourselast,
                                    state->protein_conc);
      /*active_vect(genotype[0], state->protein_conc, genesActive);
      for(j=0; j<NGENES; j++){
               revise_activity_state(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
       }*/
      break;
    default:
      LOG_ERROR("event=%d should never get here\n", event);
      exit(-1);
      break;
    }
    
    *t += *dt;                  /* advance time by the dt */
    //*x -= (*dt)*(*konrate);

    //LOG_VERBOSE("dt=%g t=%g fixed event old x=%g new x=%g\n", *dt, *t, (*x)+(*dt)*(*konrate), *x);
    
    *dt = expdev(&seed);
    /* re-compute a new dt */
    /*calc_dt(x, dt, rates, //kon_states,
     mRNAdecay, 
            genotype->mRNAdecay, state->mRNA_cyto_num, state->mRNA_transl_cyto_num, state->cell_id);*/
    
    LOG_VERBOSE("next stochastic event (2) due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

    fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);    

    /* check to see there aren't more fixed events to do */
    event = does_fixed_event_end(state->mRNA_transl_time_end, 
                                 state->mRNA_transcr_time_end, 
                                 state->replication_time_end,
                                 fixed_time);
  } 
  recalibrate_cell(rates, state, genotype, mRNAdecay, transport, *dt);
  LOG_ERROR("RATES SUBTOTAL = %f\n", rates->subtotal);
  // printf("no more events");
  /* no remaining fixed events to do in dt, now do stochastic events */
  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */
  if (*t+*dt < tdevelopment || no_fixed_dev_time) {
    //printf("inside if");
    /* compute konrate */
    /*if (kon_states->nkon==0) {
      *konrate = (-rates->salphc);    // all binding sites are occupied, total rate of binding is zero 
      LOG_WARNING("kon_states->nkon=%d, konrate=%g\n", kon_states->nkon, *konrate);
    } else  {
      calc_kon_rate(*dt, kon_states, konrate);  // otherwise compute konrate 
    }*/

    /* if the total rates falls below zero, we do an emergency recalibration of cell */
    if (!(rates->subtotal > 0.0)) {//+ *konrate 
     //log_snapshot(rates, state, genotype, kon_states, &koffvalues, mRNAdecay, transport, *konrate, *x, *t);
      LOG_ERROR("x should always be >0 t=%g (x=%g) rates->subtotal=%g, recalibrate cell!\n", *t, *x, rates->subtotal); 
      recalibrate_cell(rates, state, genotype, mRNAdecay, transport, *dt); // kon_states,&koffvalues
     // log_snapshot(rates, state, genotype, kon_states, &koffvalues, mRNAdecay, transport, *konrate, *x, *t);*/
     }

      /* if this still results in either zero or negative total rates,
         this most likely due the cell being "dead" no TFs bound, no
         activity etc.  We mark cell as "dead" in this case, and
         remove from queue. */
         //recalibrate_cell(rates, state, genotype, kon_states, &koffvalues, mRNAdecay, transport, *dt);
      if (!(rates->subtotal > 0.0)) {  //+ *konrate 
        LOG_ERROR("cell is effectively dead\n"); 
        return -1;        /* return cell status as "dead" */
      }
//    }
    
    /* 
     * choose a new uniform random number weighted by the rate of all
     * Gillespie events, note that konrate is *not* included in
     * rates->subtotal, so it needs to be added here
     */
    // rates->subtotal = 104833.04;
     LOG_ERROR("\n HERE!!\n");
    for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
       //LOG_ERROR("rates->deacetylation_num[%d]=%d\n", j, rates->deacetylation_num[j]);
    }
    LOG_ERROR("\n HERE!!\n");
    
    LOG_ERROR("rates subtotal = %f\n", rates->subtotal);
    LOG_ERROR("koff = %f, transport = %f, mrnadecay = %f, pic dis = %f, salphc = %f, max s = %f, min s = %f, \n", 
    rates->koff, rates->transport, rates->mRNAdecay, rates->pic_disassembly, rates->salphc, rates->max_salphc, rates->min_salphc);
    for(i=0; i<MAX_COPIES; i++){
             LOG_ERROR("ace_num[%d] = %d, deace_num[%d] = %d, pic[%d] = %d, trans_num[%d] = %d, pic_dis[%d] = %d\n",
             i,rates->acetylation_num[i],i,rates->deacetylation_num[i],i, rates->pic_assembly_num[i],i, rates->transcript_init_num[i],i, rates->pic_disassembly_num[i]);
    }
 
 
    *x = ran1(&seed)*rates->subtotal;//*(876.76);//*(rates->subtotal);  // + *konrate
    LOG_ERROR("NEW x = %f\n", *x);
    float test2, *testProb, *deaceProb;
    testProb = malloc(NGENES*sizeof(float));
    deaceProb = malloc(NGENES*sizeof(float));
    test2 = *x; /// rates->subtotal;
    LOG_ERROR("test2 = %f\n", test2);
    for(i=0;i<NGENES; i++){
       testProb[i] = genesActive[i]*ACETYLATE*(float)sum_rate_counts(rates->acetylation_num)*10;
       //deaceProb[i] = (1-genesActive[i])*DEACETYLATE*(float)sum_rate_counts(rates->deacetylation_num);
       LOG_ERROR("testProb[%d] = %f\n",i, testProb[i]);
    }
    for(i=0;i<NGENES; i++){
         LOG_ERROR("deaceProb[%d] = %f\n",i, testProb[i]);  
    }
    for(i=0; i<NGENES; i++){    
    if(test2 < testProb[i]){
             revise_activity_state(i,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
    }
    }
    LOG_ERROR("\n HERE!!\n");
    
    for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
       //LOG_ERROR("rates->deacetylation_num[%d]=%d\n", j, rates->deacetylation_num[j]);
    }
    LOG_ERROR("\n HERE!!\n");
    /*if(test2 <deaceProb[i]){
             LOG_ERROR("DEACE %d!!\n", i);
             revise_activity_state(i,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
    }
    }*/
    
    //testProb = genesActive[0]*ACETYLATE*(float)sum_rate_counts(rates->acetylation_num);
    //LOG_ERROR("testProb = %f\n", testProb);
   // printf("x* = %f \n rates_transport = %f \n rates subtotal =%f \n rates->mRNAdecay= %f \n rates->pic_disassembly= %f \n sumrate = %f\n", *x, rates->transport, rates->subtotal, rates->mRNAdecay,  rates->pic_disassembly,  (float) sum_rate_counts(rates->acetylation_num) * ACETYLATE);
    //LOG_ERROR("*x = %f, rates subtotal = %f\n", *x, rates->subtotal);
    LOG_ERROR("rates transport = %f\n rates mRNAdecay = %f\n", rates->transport, rates->mRNAdecay);
    LOG_ERROR("rates pic = %f\n sum rate ACE = %f\n sum rate DEACE = %f\n sum rate PIC = %f\n sum rate trans = %f\n", rates->pic_disassembly, sum_rate_counts(rates->acetylation_num) * ACETYLATE, sum_rate_counts(rates->deacetylation_num) * DEACETYLATE, sum_rate_counts(rates->pic_assembly_num) * PICASSEMBLY, sum_rate_counts(rates->transcript_init_num) * TRANSCRIPTINIT);
    
    
    if (verbose) {
      log_snapshot(rates,
                   state,
                   genotype,
                   
                   mRNAdecay,
                   transport, 
                  
                   *x,
                   *t);
    }//kon_states, //&koffvalues,  //*konrate,
    
     float test;
      test = active_to_repress( *genotype, state->protein_conc, 0,  0);

     LOG_ERROR("a to r = %f, ace = %f\n", genesActive[0], (float)sum_rate_counts(rates->acetylation_num)*ACETYLATE);
     //*x = ran
     LOG_ERROR("x = %f\n", *x);
     
     float part1, part2, part3, part4, part5, part6, part7;
     part1 = .3;
     part2= .5;
     part3= .53;
     part4= .7;
     part5= .72;
     part6= .98;
     part7 =1;
     // if (*x < rates->transport) { 
           if (*x < part1) { 
             LOG_ERROR("transport event\n");  
             //LOG_ERROR("Inside if statement\n"); 
               
        transport_event(rates, state, genotype,  transport, 
                        timecoursestart, timecourselast, *dt, *t, *x);//kon_states,
      
      } else {
        //*x -= part1;
        //*x -= rates->transport;
        /* 
         * STOCHASTIC EVENT: an mRNA decay event
         */
       // if (*x < rates->mRNAdecay) {  
             if ( *x < part2) { 
               LOG_ERROR("decay event\n");
          mRNA_decay_event(rates, state, genotype,  mRNAdecay,
                           timecoursestart, timecourselast, *dt, *t, *x);//kon_states,
           
        } else {
          //*x -= part2;
          //*x -= rates->mRNAdecay;
          /* 
           * STOCHASTIC EVENT: PIC disassembly
           */
          //if (*x < rates->pic_disassembly) {
               if (*x < part3) {
                  LOG_ERROR("pic disassembly event\n");
            disassemble_PIC_event(rates, state, genotype,  
                                  timecoursestart,  timecourselast, *dt, *t, *x);//kon_states,
          } else {
           // *x -= part3;
            //*x -= rates->pic_disassembly;
            
               //LOG_ERROR("rates salphc = %f\n", rates->salphc);
               LOG_ERROR("transcript init num = %d, pic assemb = %d, pic assem *PIC = %f\n", sum_rate_counts(rates->transcript_init_num), sum_rate_counts(rates->pic_assembly_num), sum_rate_counts(rates->pic_assembly_num)*PICASSEMBLY);
               LOG_ERROR("CHECK! x = %f, sum_rate_counts = %d, acetylate = %f\n", *x, sum_rate_counts(rates->acetylation_num), ACETYLATE );
            
            
            
             // if (*x < (float) sum_rate_counts(rates->acetylation_num) * ACETYLATE) {
                   if (*x < part4){
                LOG_ERROR("hist act event\n");
                histone_acteylation_event(rates, state, genotype, //kon_states, 
                                          timecoursestart, timecourselast, *dt, *t);
                //active_vect(genotype[0], state->protein_conc, genesActive);
              } else {
                LOG_ERROR("x = %f\n", *x);
               
               
              // *x -= part4;
               // *x -= (float) sum_rate_counts(rates->acetylation_num) * ACETYLATE;
                LOG_ERROR("x = %f\n", *x);
                /* 
                 * STOCHASTIC EVENT: histone deacetylation
                 */
                  LOG_ERROR("CHECK! x = %f, sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
               
               
               
               if (*x < part5){
               // if (*x < (float) sum_rate_counts(rates->deacetylation_num) * DEACETYLATE) {
                   LOG_ERROR("deact event\n");
                  histone_deacteylation_event(rates, state, genotype, //kon_states, 
                                              timecoursestart, timecourselast, *dt, *t);
                 
                  //active_vect(genotype[0], state->protein_conc, genesActive);
                } else {
                 
                 //*x -= part5;
                 // *x -= (float) sum_rate_counts(rates->deacetylation_num) * DEACETYLATE;
                  /* 
                   * STOCHASTIC EVENT: PIC assembly
                   */
                 
                 if (*x < part6){
                 // if (*x < (float) sum_rate_counts(rates->pic_assembly_num) * PICASSEMBLY) {
                    LOG_ERROR("pic assembly event\n"); 
                    assemble_PIC_event(rates, state, genotype, //kon_states, 
                                       timecoursestart, timecourselast, *dt, *t);
                    //active_vect(genotype[0], state->protein_conc, genesActive);
                  } else {
                   
                   
                   
                  // *x -= part6;
                   // *x -= (float) sum_rate_counts(rates->pic_assembly_num) * PICASSEMBLY;
                    /* 
                     * STOCHASTIC EVENT: transcription initiation
                     */
                     LOG_ERROR("transcript init num = %d\n", sum_rate_counts(rates->transcript_init_num));
                   
                   
                   if (*x < part7){
                   // if (*x < (float) sum_rate_counts(rates->transcript_init_num) * TRANSCRIPTINIT) {
                       LOG_ERROR("transcript init event time = %f\n", *t);
                      transcription_init_event(rates, state, genotype, //kon_states, 
                                               timecoursestart, timecourselast, *dt, *t, *x);
                      //active_vect(genotype[0], state->protein_conc, genesActive);
                    } else {
                      /*
                       * FALLBACK: shouldn't get here, previous
                       * events should be exhaustive
                       */
                      
                      LOG_ERROR("[cell %03d] t=%g no event assigned: x=%g, rates->subtotal+konrate=%g, recalibrate cell\n", 
                                state->cell_id, *t, *x, rates->subtotal);

                      log_snapshot(rates, state, genotype, 
                      mRNAdecay,transport,  *x, *t);//kon_states, &koffvalues, *konrate,
                      recalibrate_cell(rates, state, genotype, mRNAdecay, transport, *dt); //kon_states, &koffvalues,
                      log_snapshot(rates, state, genotype, 
                      mRNAdecay, transport, *x, *t);//kon_states, &koffvalues, *konrate, 
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
   // printf("rates transportHere = %f\n", rates->transport);
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
                                  rates, //kon_states,
                                  *t, timecoursestart, timecourselast,
                                  state->protein_conc);
    //printf("last active_vect\n");
    //active_vect(genotype[0], state->protein_conc, genesActive);
    
    /* advance to end of development (this exits the outer while loop) */
    *t = tdevelopment;
  }
   /*active_vect(genotype[0], state->protein_conc, genesActive);
     for(j=0; j<NGENES; j++){
    revise_activity_state(j,genotype->all_binding_sites->gene_copy, genotype, state,rates, genesActive);
    }*/
  for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("4rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      }
  return 0;
}

/*
 * initialize the population of cells, then run for a given
 * number of divisions or run cell(s) for a specific length of time
 */
 //DESTROY
void init_run_pop(Genotype genotype[POP_SIZE],
                  CellState state[POP_SIZE],
                  TimeCourse *timecoursestart[POP_SIZE][NGENES],
                  TimeCourse *timecourselast[POP_SIZE][NGENES],
                  float temperature,   /* in Kelvin */
                  float kdis[NUM_K_DISASSEMBLY],
                  int output_binding_sites,
                  int no_fixed_dev_time,
                  int max_divisions)
{
  /* local variables that don't require per-cell tracking */
  int i, j;
  int maxbound2, maxbound3;  

  /* initial mRNA and protein concentrations */
  float init_mRNA[NGENES], init_protein_conc[NGENES];

  /* cached information about available binding sites for efficiency */
  //KonStates kon_states[POP_SIZE];
  GillespieRates rates[POP_SIZE];

  float t[POP_SIZE];             /* time of last event */
  //float *koffvalues[POP_SIZE];   /* rates of unbinding */
  float transport[POP_SIZE][NGENES];  /* transport rates of each mRNA */
  float mRNAdecay[POP_SIZE][NGENES];  /* mRNA decay rates */
  float x[POP_SIZE];                  /* random number */
  float dt[POP_SIZE];                 /* delta-t */
 // float konrate[POP_SIZE];

  /* priority queue and time step initialization */
  bheap_t *priority_queue, *reaper_queue; /*living cells and empty slots*/
  float t_next = 0.0;  
  int cell = 0;         /* next cell to have an event */
  int divisions = 0;    /* no cell divisions yet */
  int mother_id;         /* mother cell at division */
  int daughter_id;       /* daughter cell at division */
  int large_cell_size = 0;  /* count the number of times cell_size exceeds Y=2.0*/
  
  float *genesActive = malloc(NGENES*sizeof(float));

  maxbound2 = MAXBOUND;
  maxbound3 = 10*MAXBOUND;

  priority_queue = bh_alloc(POP_SIZE);
  reaper_queue = bh_alloc(POP_SIZE);

  clock_t start, stop;
  double tt = 0.0;

  /* initialize protein concentrations to be used in all genes */
  // TODO: FIXME ASAP init_protein_conc should be over NPROTEINS and
  // init_mRNA should only loop over NGENES as they could be different
  for (i=0; i < NGENES; i++) {
    init_protein_conc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    init_mRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }
  
  for(i=0; i<NGENES; i++){
    printf("initProtein Conc %f\n", init_protein_conc[i]);
    }
 

  /* initialize the population of cells */
  for (j = 0; j < POP_SIZE; j++) {

    output=1;
    /* for all cells in this replicate initialize all parts of the
       genotype, *except* for the cis-reg sequence using the initial
       genes[0]  */
    initialize_genotype(&genotype[j], &genotype[0], kdis, j);
    initialize_cell(&state[j], j, genotype[j].copies, genotype[j].mRNAdecay, init_mRNA, init_protein_conc, burn_in);
    //printf("proteindecay[0] = %f\n", genotype[j].proteindecay[]);

    /* print binding sites */
    if (output_binding_sites) 
      print_all_binding_sites(genotype[j].copies, genotype[j].all_binding_sites, genotype[j].binding_sites_num, 
                              genotype[j].tf_seq, genotype[j].cisreg_seq, genotype[j].site_id_pos); 

    
    /* set cell temperature and value of RTlnKr constant */
    state[j].temperature = temperature;
    state[j].RTlnKr = GASCONSTANT * temperature * log(KR);

    /* initialize time courses */
    for (i=0; i < NPROTEINS; i++){
      timecoursestart[j][i] = NULL;
      timecourselast[j][i] = NULL;
    } 

    add_time_points((float) 0.0, state[j].protein_conc, timecoursestart[j], timecourselast[j]);

    for(i=0; i<NGENES; i++){
             LOG_VERBOSE("mRNAdecay[%d] = %g\n", i, genotype->mRNAdecay[i]);
    }

    /* initialize KonStates data structures */
    //CHANGE JU
    //initialize_cell_cache(&(state[j]), genotype[j], &(kon_states[j]), &(koffvalues[j]), maxbound2, maxbound3);

    /* initialize transcriptional state of genes and computes initial rates */
    //QUESTION
 
    /*init_protein_conc[0] = 1134.99;
    init_protein_conc[1] = 3223.03;
    init_protein_conc[2] = 713.07;
    init_protein_conc[3] = 1809.57;
    init_protein_conc[4] = 4540.77;
    init_protein_conc[5] = 6446.1;
    init_protein_conc[6] = 9715.71;
    init_protein_conc[7] = 257.74;
    init_protein_conc[8] = 1512.38;
    init_protein_conc[9] = 303.28;*/
    // = [1134.99, 3223.03, 713.07,1809.57,4540.77,6446.1,9715.71,257.74,1512.38,5113.58,303.28];
    /*for(i=0; i<NGENES; i++){
    printf("initProtein Conc %f\n", init_protein_conc[i]);
    }*/
    //active_vect(genotype[0], init_protein_conc, genesActive);
    /*
    for(i=0; i<NGENES; i++){
           printf( " WHAT = %f\n", genesActive[i] );
    }*/
    //system("PAUSE");
    calc_from_state(&genotype[j], &state[j], &rates[j],  transport[j], mRNAdecay[j]);//&kon_states[j],

    //active_vect(genotype[0], init_protein_conc, genesActive);
    
    t[j] = 0.0;  /* time starts at zero */
    state[j].division_time = TIME_INFINITY;  /* make artificially high */
   //printf("rates->transport = %f\n", rates->transport);
    /* do one step for the cell */
    if((start = clock()) != -1){
    do_single_timestep(&(genotype[j]), 
                       &(state[j]), 
                       //&(kon_states[j]), 
                       &(rates[j]), 
                       &(t[j]),
                       //koffvalues[j],
                       transport[j],
                       mRNAdecay[j],
                       &(x[j]),
                       &(dt[j]),
                       //&(konrate[j]),
                       timecoursestart[j],
                       timecourselast[j],
                       maxbound2,
                       maxbound3,
                       no_fixed_dev_time,
                       genesActive);
    stop = clock();
    tt = (double)(stop-start)/CLOCKS_PER_SEC;
    //printf("Run time: %f\n", tt);
    //system("PAUSE");
    }
    //RUN MY CODE HERE 
    /* insert each initial time into the priority_queue */
    insert_with_priority_heap(priority_queue, j, t[j]);
    LOG_NOCELLID("[cell=%03d] inserted at time=%g\n", j, t[j]); 
  }
    
    for(i=0; i<NGENES; i++){
             LOG_VERBOSE("mRNAdecay[%d] = %g\n", i, genotype->mRNAdecay[i]);
    }
  //printf("rates->transport = %f\n", rates->transport);
  //active_vect(genotype[0], init_protein_conc, genesActive);
  

  while ((no_fixed_dev_time && divisions < max_divisions) ||    /* no fixed dev time, run until # divisions reached */
         (!no_fixed_dev_time && t_next < tdevelopment)) {       /* or, if fixed dev time, run until tdevelopment reached */
                                                                /* TODO: when evolutionary goal is achieved! */
    int cell_status; 
    
    LOG_VERBOSE_NOCELLID("before choosing next cell (queue len reaper=%03d, main=%03d)\n", reaper_queue->n, priority_queue->n);
    /* get the next cell with the smallest t to advance next */
    if (reaper_queue->n == POP_SIZE || priority_queue->n == 0) {
      /* if all cells in population are dead, we exit */
      LOG_NOCELLID("all %03d cells in population are dead, main queue len=%03d\n", reaper_queue->n, priority_queue->n);
      break;
    } 
    /* otherwise get next event from priority_queue */
    t_next = get_next_heap(priority_queue, &cell);
    LOG_VERBOSE_NOCELLID("[cell=%03d] get minimum time=%g, size=%g (queue len reaper=%03d, main=%03d)\n", 
                         cell, t_next, state[cell].cell_size, reaper_queue->n, priority_queue->n);

     if (1) //VERBOSE JU CHANGE
    for (j=0; j < MAX_COPIES; j++) {
      //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
      LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
    }
                        
    //printf("cellStatus\n");
    /* update a single timestep on this cell, returns new timestep in t[cell] */
    cell_status = do_single_timestep(&(genotype[cell]), 
                                     &(state[cell]), 
                                     //&(kon_states[cell]), 
                                     &(rates[cell]), 
                                     &(t[cell]),
                                     //koffvalues[cell],
                                     transport[cell],
                                     mRNAdecay[cell],
                                     &(x[cell]),
                                     &(dt[cell]),
                                     //&(konrate[cell]),
                                     timecoursestart[cell],
                                     timecourselast[cell],
                                     maxbound2,
                                     maxbound3,
                                     no_fixed_dev_time,
                                     genesActive);
   // printf("rates->transport = %f\n", rates->transport);
    if (cell_status == -1)  {   /* if cell is dead */
      /* for debugging purposes, print out protein time courses when a cell is found to be dead */
      print_all_protein_time_courses(timecoursestart, timecourselast);
      /* add this cell to list of empty cell locations */
      insert_with_priority_heap(reaper_queue, cell, TIME_INFINITY);
      LOG_NOCELLID("[cell %03d] added here to reaper_queue as a dead cell, t_next=%g, t=%g (queue len reaper=%03d, main=%03d)\n", 
                   cell, t_next, t[cell], reaper_queue->n, priority_queue->n);
      /* this cell is dead, so don't add back to main priority queue, skip to next event */
      continue;
    }
   
    if (t[cell] >= (state[cell]).division_time)  {  /* we have now reached cell division */
      int keep_mother, keep_daughter;
      float current_division_time = (state[cell]).division_time;  // store current division time
      divisions++;     /* increment number of divisions */
      mother_id = cell; /* set mother cell as currently dividing cell */

      if (state[cell].cell_size > 2.0) {
        large_cell_size++;
        printf("[cell %03d] at t=%g cell size =%g, exceeding Y=2.0, %d of %d divisions (%g fraction)\n", 
               cell, t[cell], state[cell].cell_size, large_cell_size, divisions, (double)large_cell_size/(double)divisions);

        LOG_NOCELLID("[cell %03d] at t=%g cell size =%g, exceeding Y=2.0, %d of %d divisions (%g fraction)\n", 
                     cell, t[cell], state[cell].cell_size, large_cell_size, divisions, (double)large_cell_size/(double)divisions);
        // TODO: currently this just keeps track of the number of divisions that the cell exceeds 2.0
        // to implement:
        //  0. let this number go down after a fat cell divides
        //  1. a check that a threshold has been reached
        //  2. choose a new, lower growth_rate_scaling
        //  3. re-run initialize_growth_rate_parameters() to recompute gpeak etc.
        //  4. for each cell in the population, run recalibrate_cell() to compute new rates etc.
        //  5. now have machinery for simulated annealing, get evolution to work.
      }

      /* to get new daughter cell slot, first check list of empty cells */
      if (reaper_queue->n > 0) {
        get_next_heap(reaper_queue, &daughter_id);  /* use one of the empty cells */
        keep_mother = 1;                            /* keep both */
        keep_daughter = 1;                   
        LOG_NOCELLID("[cell %03d] in reaper_queue (length=%03d) used as daughter cell\n", daughter_id, reaper_queue->n);
      } else {
        /* otherwise choose one of the now POP_SIZE + 1 (since the mother has reproduced) other cell randomly to die */
        daughter_id = trunc((POP_SIZE+1)*ran1(&seed));  /*for Alex: regression change here*/

        if (daughter_id == mother_id) {          /* daughter replaces mother */
          LOG_NOCELLID("daughter=%03d replaces mother=%03d\n", daughter_id, mother_id);
          keep_mother = 0;                     /* discard mother */
          keep_daughter = 1;                   /* keep daughter */

        } else if (daughter_id == POP_SIZE) {    /* arbitrarily define POP_SIZE as being mother replaces daughter */
          LOG_NOCELLID("mother=%03d replaces daughter=%03d\n", mother_id, daughter_id);
          keep_mother = 1;                     /* keep mother */
          keep_daughter = 0;                   /* discard daughter */

        } else {  /* removing pending event in original daughter cell
                     from queue, since we are replacing it only remove
                     if the daughter cell wasn't already removed from queue */
          keep_mother = 1;                     /* keep both mother and daughter */
          keep_daughter = 1;                   
          delete_element_heap(priority_queue, daughter_id);
          LOG_NOCELLID("removing pending event from queue (length=%3d) replaced by daughter cell=%d\n", priority_queue->n, daughter_id);
        } 
      }

      printf("[cell %03d] (size=%g) dividing into mother=%03d and daughter=%03d at t=%g, division=%g, total divisions=%d (t_next=%g)\n", 
             cell, state[cell].cell_size, mother_id, daughter_id, t[cell], current_division_time, divisions, t_next);
      LOG_NOCELLID("[cell %03d] (size=%g) dividing into mother=%03d and daughter=%03d at t=%g, division=%g, total divisions=%d (t_next=%g)\n", 
                   cell, state[cell].cell_size, mother_id, daughter_id, t[cell], current_division_time, divisions, t_next);

      do_cell_division(mother_id, daughter_id,
                         keep_mother, keep_daughter,
                         &(genotype[mother_id]),
                         &(state[mother_id]),
                         &(rates[mother_id]),
                         //&(kon_states[mother_id]),
                         //&(koffvalues[mother_id]),
                         transport[mother_id],
                         mRNAdecay[mother_id],
                         &(genotype[daughter_id]),
                         &(state[daughter_id]),
                         &(rates[daughter_id]),
                         //&(kon_states[daughter_id]),
                         //&(koffvalues[daughter_id]),
                         transport[daughter_id],
                         mRNAdecay[daughter_id],
                         0.44,
                         x[mother_id],
                         dt[mother_id]);
      if (keep_mother) {
        t[mother_id] = current_division_time;   /* reset current time in mother cell to division time */
        LOG_NOCELLID("AFTER DIVISION: time at instant of division for mother cell=%03d is t=%g\n", mother_id, t[mother_id]);
        /* advance time for mother */
        cell_status = do_single_timestep(&(genotype[mother_id]), 
                                         &(state[mother_id]), 
                                         //&(kon_states[mother_id]), 
                                         &(rates[mother_id]), 
                                         &(t[mother_id]),
                                         //koffvalues[mother_id],
                                         transport[mother_id],
                                         mRNAdecay[mother_id],
                                         &(x[mother_id]),
                                         &(dt[mother_id]),
                                         //&(konrate[mother_id]),
                                         timecoursestart[mother_id],
                                         timecourselast[mother_id],
                                         maxbound2,
                                         maxbound3,
                                         no_fixed_dev_time,
                                         genesActive);

        if (cell_status == -1) {  /* if cell dies, add it to the reaper queue */
          insert_with_priority_heap(reaper_queue, mother_id, TIME_INFINITY);
          LOG_NOCELLID("AFTER DIVISION: mother cell %03d is dead at t=%g (queue len reaper=%03d, main=%03d)\n", 
                       mother_id, t[mother_id], reaper_queue->n, priority_queue->n);
        } else {                  /* otherwise we add the new mother timestep back to the priority queue */
          LOG_NOCELLID("AFTER DIVISION: add new timestep to queue for mother cell=%03d at t=%g (queue len reaper=%03d, main=%03d)\n", 
                       mother_id, t[mother_id], reaper_queue->n, priority_queue->n);
          insert_with_priority_heap(priority_queue, mother_id, t[mother_id]); 
        }
        
      } else {
        printf("don't update mother [%03d] it replaced by daughter [%03d]\n", mother_id, daughter_id);
        LOG_NOCELLID("don't update mother [%03d] it replaced by daughter [%03d]\n", mother_id, daughter_id);
      }

      if (keep_daughter) {   /* only update if daughter cell isn't replaced by mother */
        t[daughter_id] = current_division_time;   /* reset current time in daughter cell to division time */
        LOG_NOCELLID("AFTER DIVISION: time at instant of division for daughter cell=%03d is t=%g\n", daughter_id, t[daughter_id]);

        /* advance time for daughter */
        cell_status = do_single_timestep(&(genotype[daughter_id]), 
                                         &(state[daughter_id]), 
                                         //&(kon_states[daughter_id]), 
                                         &(rates[daughter_id]), 
                                         &(t[daughter_id]),
                                         //koffvalues[daughter_id],
                                         transport[daughter_id],
                                         mRNAdecay[daughter_id],
                                         &(x[daughter_id]),
                                         &(dt[daughter_id]),
                                         //&(konrate[daughter_id]),
                                         timecoursestart[daughter_id],
                                         timecourselast[daughter_id],
                                         maxbound2,
                                         maxbound3,
                                         no_fixed_dev_time,
                                         genesActive);

        if (cell_status == -1) {   /* if cell dies, add to reaper queue */
          insert_with_priority_heap(reaper_queue, daughter_id, TIME_INFINITY);
          LOG_NOCELLID("AFTER DIVISION: daughter cell %03d is dead at t=%g (queue len reaper=%03d, main=%03d)\n", 
                       daughter_id, t[daughter_id], reaper_queue->n, priority_queue->n);
        } else {                   /* otherwise we add the new daughter timestep back to the priority queue */
          LOG_NOCELLID("AFTER DIVISION: add new timestep in queue for daughter cell=%03d at t=%g (queue len reaper=%03d, main=%03d)\n", 
                       daughter_id, t[daughter_id], reaper_queue->n, priority_queue->n);
          insert_with_priority_heap(priority_queue, daughter_id, t[daughter_id]);
        }
      } else {
        printf("don't update daughter [%03d] it is replaced by mother [%03d]\n", daughter_id, mother_id);
        LOG_NOCELLID("don't update daughter [%03d] it is replaced by mother [%03d]\n", daughter_id, mother_id);
      }

    } else {   /* no division happens at all */
      /* put the updated timestep back into the priority_queue  */
      insert_with_priority_heap(priority_queue, cell, t[cell]);
    }
  }

  /* output the founder information */
  for (j = 0; j < POP_SIZE; j++) {
    printf("cell %03d derived from founder %03d had %2d divisions%s\n", 
           j, state[j].founder_id, state[j].divisions, t[j] == TIME_INFINITY ? " [dead]": "");
    LOG_NOCELLID("cell %03d derived from founder %03d had %2d divisions%s\n", 
                 j, state[j].founder_id, state[j].divisions, t[j] == TIME_INFINITY ? " [dead]": "");
  }

  /* cleanup data structures */
  bh_free(priority_queue);
  bh_free(reaper_queue);

 // for (j = 0; j < POP_SIZE; j++) {
    //free(koffvalues[j]);
   //for (i=0; i < NPROTEINS; i++) {
     // free(kon_states[j].kon_list[i]->available_sites);
      //free(kon_states[j].kon_list[i]);
    //}
  //}
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




