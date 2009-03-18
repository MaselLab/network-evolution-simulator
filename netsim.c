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

const float mN = 0.1;
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

  printf("[genotype %03d] hindPos: ", genotypeID);
  for (i=0; i < TFGENES; i++) {
    printf("%d ", genotype->hindrancePositions[i]);
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
      printf("PICd[%d]=%g ", p, genotype->PICdisassembly[i][p]);
    }
    printf("\n");
  }
}


void print_all_binding_sites(int copies[NGENES],
                             AllTFBindingSites *allBindingSites, 
                             int numElements,
                             char transcriptionFactorSeq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                             char cisRegSeq[NGENES][MAX_COPIES][CISREG_LEN],
                             int tfsStart[NGENES][MAX_COPIES][2])
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
        printf("TF sequence gene %2d (copy %d): %.*s\n", i, j, TF_ELEMENT_LEN, transcriptionFactorSeq[i][j]);
      else
        printf("            gene %2d (copy %d) does not encode a TF\n", i, j);
      if (i < NGENES) {
        printf("cis-reg     gene %2d (copy %d): %.*s\n", i, j, CISREG_LEN, cisRegSeq[i][j]);
        printf("ID range    gene %2d (copy %d): [%3d, %3d]\n", i, j, tfsStart[i][j][0], tfsStart[i][j][1]);
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
    printf("       cis-reg region: %3d", allBindingSites[i].cisregID);
    printf("         cis-reg copy: %3d", allBindingSites[i].geneCopy);
    printf(" (sequence %.*s)\n", CISREG_LEN, cisRegSeq[allBindingSites[i].cisregID][allBindingSites[i].geneCopy]);
    printf(" transcription-factor: %3d", allBindingSites[i].tfID);
    printf(" (sequence: %.*s)\n", TF_ELEMENT_LEN, transcriptionFactorSeq[allBindingSites[i].tfID][allBindingSites[i].geneCopy]); 
    printf("  L-edge of %2dbp hind: %3d\n", HIND_LENGTH, allBindingSites[i].leftEdgePos);        
    printf("  Hind offset position: %3d\n", allBindingSites[i].hindPos); 
    printf("               strand: %3d\n", allBindingSites[i].strand);
    printf("         Hamming dist: %3d\n", allBindingSites[i].hammingDist); 
  }
}

void print_tf_occupancy(CellState *state,
                        AllTFBindingSites *allBindingSites,
                        float t)
{
  int bound_count[NGENES][MAX_COPIES];
  int i, j, geneID, geneCopy;

  for (i = 0; i < NGENES; i++) 
    for (j = 0; j < MAX_COPIES; j++) 
      bound_count[i][j] = 0;
  
  /* for all the currently bound sites */
  for (j = 0; j < state->tfBoundCount; j++) {
    geneID = allBindingSites[state->tfBoundIndexes[j]].cisregID;
    geneCopy = allBindingSites[state->tfBoundIndexes[j]].geneCopy;
    bound_count[geneID][geneCopy]++;
  }

  // TODO: currently comment-out generation of TF bound information
  fprintf(fp_tfsbound[state->cellID], "%g %d ", t, state->tfBoundCount);
  for (i = 0; i < NGENES; i++) 
    for (j = 0; j < MAX_COPIES; j++) 
      fprintf(fp_tfsbound[state->cellID], "%d ", bound_count[i][j]);

  fprintf(fp_tfsbound[state->cellID], "\n");
}

void print_rounding(CellState *state, GillespieRates *rates, float t)
{
  fprintf(fp_rounding[state->cellID], "%g %d %d %d %d %d %d %d %d %d %d %d %d\n", 
          t, rates->koff_operations, rates->transport_operations, rates->mRNAdecay_operations, 
          rates->picDisassembly_operations, rates->salphc_operations, rates->maxSalphc_operations, rates->minSalphc_operations, 
          rates->acetylationCount_operations, rates->deacetylationCount_operations, rates->picAssemblyCount_operations, 
          rates->transcriptInitCount_operations, rates->picDisassemblyCount_operations);
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
      indiv->hindrancePositions[p]=0;
    } else  {
      // TODO: only keep until new regression output is generated
#ifdef USE_RAND
      indiv->hindrancePositions[p]=rand()%10;
#else
      indiv->hindrancePositions[p]=rint(ran1(&seed)*(HIND_LENGTH - TF_ELEMENT_LEN));
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
      indiv->PICdisassembly[i][p] = kdis[j];
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

void initialize_genotype(Genotype *indiv, 
                         Genotype clone,
                         float kdis[],
                         int genotypeID)
{
  int i, k, p;
  
  initialize_sequence((char *)indiv->cisRegSeq, CISREG_LEN*MAX_COPIES*NGENES, MAX_COPIES, NGENES);

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
    initialize_sequence((char *)indiv->transcriptionFactorSeq, TF_ELEMENT_LEN*MAX_COPIES*TFGENES, MAX_COPIES, TFGENES);
    initialize_genotype_fixed(indiv, kdis, genotypeID);
  // otherwise clone it from the first instance
  } else {
    // TODO: FIXME: should refactor this to have common code to clone TF sequences
    for (i=0; i < TFGENES; i++) {
      for (k=0; k < TF_ELEMENT_LEN; k++) {
        for (p=0; p < MAX_COPIES; p++) {
          indiv->transcriptionFactorSeq[i][p][k] = clone.transcriptionFactorSeq[i][p][k];
        }
      }
    }
    initialize_new_cell_genotype(indiv, clone);
  }

  /* start number of copies of gene at current_ploidy */
  for (p=0; p < NGENES; p++) {
    indiv->copies[p] = current_ploidy;
  }

  calc_all_binding_sites(indiv->copies, indiv->cisRegSeq, indiv->transcriptionFactorSeq, 
                         &(indiv->bindSiteCount), &(indiv->allBindingSites), indiv->hindrancePositions,
                         indiv->tfsPerGene, indiv->tfsStart);
}

void mutate(Genotype *gene,
            int geneID,
            int geneCopy,
            float m)
{
  int k, p;
  char x;
  float q;
  
  for (k=0; k<CISREG_LEN; k++) {

    p = geneCopy;
    if (m > ran1(&seed)) {
      x = gene->cisRegSeq[geneID][p][k]; /* because sometimes gene and new are the same */

      LOG_VERBOSE_NOCELLID("in geneID=%2d, mutating pos=%3d, copy=%1d:", geneID, k, p);
      while (gene->cisRegSeq[geneID][p][k] == x) {
        q = ran1(&seed);
        gene->cisRegSeq[geneID][p][k] = set_base_pair(q);
      }
      if (verbose)
        LOG_NOFUNC(" '%c' -> '%c'\n", x, gene->cisRegSeq[geneID][p][k]);
    }
  }
}

int calc_all_binding_sites_copy(char cisRegSeq[NGENES][MAX_COPIES][CISREG_LEN],
                                char transcriptionFactorSeq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                int bindSiteCount,
                                AllTFBindingSites **allBindingSites,
                                int *maxAlloc,
                                int geneID,
                                int geneCopy,
                                int hindPos[TFGENES])
{
  int i, j, tfind, match, maxBindingSiteAlloc;

  maxBindingSiteAlloc = *maxAlloc;

  for (i=0; i < CISREG_LEN-TF_ELEMENT_LEN; i++) {  /* scan forwards */
    for (tfind=0; tfind < TFGENES; tfind++) {      /* only loop through TF genes */
      match=0;
      for (j=i; j < i+TF_ELEMENT_LEN; j++) {
        if (cisRegSeq[geneID][geneCopy][j] == transcriptionFactorSeq[tfind][geneCopy][j-i])
          match++;
      }
      if (match >= nmin){
        if (bindSiteCount + 1 >= maxBindingSiteAlloc) {
          maxBindingSiteAlloc = 2*maxBindingSiteAlloc;
          *allBindingSites = realloc(*allBindingSites, maxBindingSiteAlloc*sizeof(AllTFBindingSites));
          if (!(*allBindingSites)) {
            LOG_ERROR_NOCELLID("realloc of allBindingSites to bindSiteCount = %d failed.\n", maxBindingSiteAlloc);
            exit(1);
          }
          else LOG_VERBOSE_NOCELLID("realloc of allBindingSites to bindSiteCount = %d succeeded\n", maxBindingSiteAlloc);
        }
        if(((i - hindPos[tfind]) >=0) && (((i - hindPos[tfind]) + (HIND_LENGTH -1)) < CISREG_LEN)) {
          (*allBindingSites)[bindSiteCount].cisregID = geneID;
          (*allBindingSites)[bindSiteCount].geneCopy = geneCopy; 
          (*allBindingSites)[bindSiteCount].tfID = tfind;
          (*allBindingSites)[bindSiteCount].strand = 0;
          (*allBindingSites)[bindSiteCount].hammingDist = TF_ELEMENT_LEN-match;
          (*allBindingSites)[bindSiteCount].hindPos = hindPos[tfind];
          (*allBindingSites)[bindSiteCount].leftEdgePos = i - hindPos[tfind];
          bindSiteCount++;
        }
      }
    }
  }
  for (i=CISREG_LEN-1; i>=TF_ELEMENT_LEN-1; i--) {  /* scan backwards */
    for (tfind=0; tfind < TFGENES; tfind++) {       /* only loop through TF genes */
      match=0;
      for (j=i; j>i-TF_ELEMENT_LEN; j--)
        if (
            (cisRegSeq[geneID][geneCopy][j]=='a' && transcriptionFactorSeq[tfind][geneCopy][i-j]=='t')
            || (cisRegSeq[geneID][geneCopy][j]=='t' && transcriptionFactorSeq[tfind][geneCopy][i-j]=='a')
             || (cisRegSeq[geneID][geneCopy][j]=='c' && transcriptionFactorSeq[tfind][geneCopy][i-j]=='g')
            || (cisRegSeq[geneID][geneCopy][j]=='g' && transcriptionFactorSeq[tfind][geneCopy][i-j]=='c')            
            ) match++;
      if (match >= nmin){
        if (bindSiteCount + 1 >= maxBindingSiteAlloc){
          maxBindingSiteAlloc = 2*maxBindingSiteAlloc;
          *allBindingSites = realloc(*allBindingSites, maxBindingSiteAlloc*sizeof(AllTFBindingSites));
          if (!(*allBindingSites)){
            LOG_ERROR_NOCELLID("realloc of allBindingSites to bindSiteCount = %d failed.\n", maxBindingSiteAlloc);
            exit(1);
          }
          else LOG_VERBOSE_NOCELLID("realloc of allBindingSites to bindSiteCount = %d succeeded\n", maxBindingSiteAlloc);
        }
        if (((i-TF_ELEMENT_LEN+1 - hindPos[tfind]) >=0) && (((i-TF_ELEMENT_LEN+1 - hindPos[tfind])+(HIND_LENGTH-1))< CISREG_LEN)) {
          (*allBindingSites)[bindSiteCount].cisregID = geneID;
          (*allBindingSites)[bindSiteCount].geneCopy = geneCopy; 
          (*allBindingSites)[bindSiteCount].tfID = tfind;
          (*allBindingSites)[bindSiteCount].strand = 1;
          (*allBindingSites)[bindSiteCount].hammingDist = TF_ELEMENT_LEN-match;
          (*allBindingSites)[bindSiteCount].hindPos = hindPos[tfind];
          (*allBindingSites)[bindSiteCount].leftEdgePos = i-TF_ELEMENT_LEN+1 - hindPos[tfind];
          bindSiteCount++;
        }
      }
    }
  }
  *maxAlloc = maxBindingSiteAlloc;
  return bindSiteCount;
}

void calc_all_binding_sites(int copies[NGENES],
                            char cisRegSeq[NGENES][MAX_COPIES][CISREG_LEN],
                            char transcriptionFactorSeq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                            int *newBindSiteCount,
                            AllTFBindingSites **allBindingSites,
                            int hindPos[TFGENES],
                            int tfsPerGene[NGENES],
                            int tfsStart[NGENES][MAX_COPIES][2])
{
  int p, maxBindingSiteAlloc, bindSiteCount;
  int geneID;

  maxBindingSiteAlloc = maxelements;
  *allBindingSites = malloc(maxBindingSiteAlloc*sizeof(AllTFBindingSites));
  if (!(*allBindingSites)) {
    LOG_ERROR_NOCELLID("initial setting of allBindingSites failed.\n");
    exit(1);
  }
  bindSiteCount = 0;

  /* initialize per gene # of sites */
  for (geneID=0; geneID < NGENES; geneID++) {
    tfsPerGene[geneID] = 0;
    for (p=0; p < MAX_COPIES; p++) {
      tfsStart[geneID][p][0] = -1;
      tfsStart[geneID][p][1] = -1;
    }
  }

  for (geneID=0; geneID < NGENES; geneID++) {  /* now which cis-reg region */
    for (p=0; p < MAX_COPIES; p++) {     /* loop through the maximum copies possible */
      if (p < copies[geneID]) {  
        int before = bindSiteCount;

        /* record initial siteID for this copy */
        tfsStart[geneID][p][0] = bindSiteCount;

        /* if this particular gene has this copy number then generate binding sites
           for all the relevant gene copies (assume no gene divergence) */
        bindSiteCount = calc_all_binding_sites_copy(cisRegSeq, 
                                                    transcriptionFactorSeq, 
                                                    bindSiteCount,
                                                    allBindingSites,
                                                    &maxBindingSiteAlloc,
                                                    geneID,
                                                    p, 
                                                    hindPos);

        /* record end siteID for this copy */
        tfsStart[geneID][p][1] = bindSiteCount - 1;

        /* add the new number of sites */
        tfsPerGene[geneID] += (bindSiteCount-before);
      }
    }
  }

  *allBindingSites = realloc(*allBindingSites, bindSiteCount*sizeof(AllTFBindingSites));
  if (!(*allBindingSites)) {
    LOG_ERROR_NOCELLID("realloc of allBindingSites down to bindSiteCount = %d failed.\n", bindSiteCount);
    exit(1);
  }
  *newBindSiteCount = bindSiteCount;
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
  newtime->geneID = i;
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
  newtime->geneID = i;
  newtime->copy = p;
  newtime->time = t;
  LOG_VERBOSE_NOCELLID("adding event end at time=%f for gene=%d (copy=%d)\n", t, i, p);
  sls_store_end(newtime, start, last);
}

void delete_fixed_event(int geneID,
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
    if ((info->geneID==geneID && info->copy==p)) {
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
                       j+1, i, geneID, p);
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

void initialize_cell(CellState *state,
                     int cellID,
                     int copies[NGENES],
                     float mRNAdecay[NGENES],
                     float meanmRNA[NGENES],
                     float initProteinConc[NPROTEINS],
                     int burn_in)
{
  int i, j, k, totalmRNA;
  float t;

  /* initialize the ID of the cell */
  state->cellID = cellID;

  /* initialize the founder cell using the current ID */
  state->founderID = cellID;

  /* don't start in S phase */
  state->in_s_phase = 0;

  /* initialize whether to do kon burn-in or not */
  state->burn_in = burn_in;

  /* start cell size at 0.5 */
  state->cellSize = 0.5;

  /* start cell with no divisions */
  state->divisions = 0;

  /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
  state->growthRate = 0.0;

  state->mRNATranscrTimeEnd = NULL;
  state->mRNATranscrTimeEndLast = NULL;
  state->mRNATranslTimeEnd = NULL;
  state->mRNATranslTimeEndLast = NULL;
  state->replicationTimeEnd = NULL;
  state->replicationTimeEndLast = NULL;

  state->tfBoundCount = 0;  /* initialize with nothing bound */
  state->tfHinderedCount = 0;
  state->tfBoundIndexes = NULL;
  state->tfHinderedIndexes = NULL;

  for (i=0; i < NGENES; i++) {

    for (j=0; j < MAX_COPIES; j++) {
      state->active[i][j] = ON_WITH_NUCLEOSOME;
    }

    totalmRNA = (int) poidev(meanmRNA[i],&seed);
    state->mRNANuclearCount[i] = (int) bnldev(startnucleus, totalmRNA, &seed);

    //LOG_VERBOSE("mRNANuclearCount=%02d for gene=%02d\n", state->mRNANuclearCount[i], i);

    state->mRNACytoCount[i] = totalmRNA - state->mRNANuclearCount[i];
    state->mRNATranslCytoCount[i] = 0;

    for (k=0; k<state->mRNACytoCount[i]; k++) {
      t = expdev(&seed) / mRNAdecay[i];
      if (t < ttranslation) {
        (state->mRNACytoCount[i])--;
        (state->mRNATranslCytoCount[i])++;
        LOG_VERBOSE("add translation event time=%g for gene=%d\n", (ttranslation-t), i);
        add_fixed_event(i, -1, ttranslation-t, &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));
      }
    } 

    int total_mRNA_transcribing = (int) poidev(meanmRNA[i]*ttranscription*mRNAdecay[i], &seed);
    
    /* split it up evenly between the copies */
    int mRNA_copy1 = trunc(total_mRNA_transcribing/current_ploidy);
    int mRNA_copy2 = total_mRNA_transcribing - mRNA_copy1;

    for (j=0; j < MAX_COPIES; j++) {
      if (j < current_ploidy)  {
        state->mRNATranscrCount[i][j] = (j==0) ? mRNA_copy1 : mRNA_copy2;
        LOG_VERBOSE_NOCELLID("initializing state->mRNATranscrCount[%2d][%2d]=%d\n", i, j, state->mRNATranscrCount[i][j]);
        for (k=0; k < state->mRNATranscrCount[i][j]; k++)
          add_fixed_event(i, j, ran1(&seed)*ttranscription, &(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));
      } else {
        state->mRNATranscrCount[i][j] = 0;
      }
    }
  }
  for (i=0; i < NPROTEINS; i++) {
    state->proteinConc[i] = initProteinConc[i];
  }
}

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
    konStates->konList[i] = malloc(sizeof(KonList));
    konStates->konList[i]->available_sites = malloc(genes.bindSiteCount*sizeof(int));
  }
  
  state->tfBoundIndexes = realloc(state->tfBoundIndexes, maxbound2*sizeof(int));
  *koffvalues = malloc(maxbound2*sizeof(float)); 
  state->tfHinderedIndexes = realloc(state->tfHinderedIndexes, 2*maxbound3*sizeof(int));
  
  if (!konStates->konvalues || !state->tfBoundIndexes || !koffvalues ||
      !state->tfHinderedIndexes || !konStates->konList) {
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
    /* if this particular TF is bound somewhere */
    if (konStates->nkonsum[i] > 0) {
      ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] * t;
      if (fabs(ct)<EPSILON) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) konStates->nkonsum[i]) * konStates->konvalues[i][KON_DIFF_INDEX] * ect;
      numer += ((float) konStates->nkonsum[i]) * konStates->konvalues[i][KON_DIFF_INDEX] * (ect-ct*exp(-ct));
    }
  }
  numer *= kon;
  r *= kon;
  denom = r + t*(rates->total + rates->salphc);
  denom = denom * denom;
  r /= t;

  /* compute delta-t */
  *f = x/(r + rates->total + rates->salphc) - t;

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
    if (konStates->nkonsum[i] > 0) {
      ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]*t;
      if (fabs(ct)<EPSILON) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) konStates->nkonsum[i])*konStates->konvalues[i][KON_DIFF_INDEX]*ect;
    }
  }
  *konrate = kon*r/t;
  LOG_VERBOSE_NOCELLID("r=%g t=%g konrate=%g\n", r, t, *konrate);
}

/* must have already updated proteinConc first */
void change_mRNA_cytoplasm(int i,
                           Genotype *genes,
                           CellState *state,
                           GillespieRates *rates,
                           KonStates *konStates)
{
  float salphc; 
  
  /* number of mRNAs in cytoplasm affects kon rates */
  salphc = (float) (state->mRNACytoCount[i]) * genes->translation[i] / konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  
  LOG_VERBOSE("change_mRNA_cytoplasm[%d]: mRNA=%d, transl rate=%g, protein decay=%g, salphc=%g\n", 
              i, state->mRNACytoCount[i], genes->translation[i], konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX], salphc);
  
  rates->salphc += konStates->nkonsum[i]*kon*(salphc - konStates->konvalues[i][KON_SALPHC_INDEX]);
  rates->salphc_operations++;

  rates->maxSalphc += konStates->nkonsum[i]*kon*(fmaxf(state->proteinConc[i], salphc) - fmaxf(state->proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]));
  rates->maxSalphc_operations++;
  rates->minSalphc += konStates->nkonsum[i]*kon*(fminf(state->proteinConc[i], salphc) - fminf(state->proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]));    
  rates->minSalphc_operations++;

  konStates->konvalues[i][KON_DIFF_INDEX] = (state->proteinConc[i] - salphc) / konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
}

void calc_koff(int k,
               AllTFBindingSites *allBindingSites,
               CellState *state,
               float *koff,
               float t)
{
  float Gibbs;  /*free energy in kJ/mol*/
  int posdiff, front, back, j;
  
  front = back = 0;
  Gibbs = (((float) allBindingSites[k].hammingDist)/3.0 - 1.0) * state->RTlnKr; /* subject to revision of TF_ELEMENT_LEN */
  for (j=0; j < state->tfBoundCount; j++) {
    if (allBindingSites[k].cisregID==allBindingSites[state->tfBoundIndexes[j]].cisregID &&
        allBindingSites[k].geneCopy==allBindingSites[state->tfBoundIndexes[j]].geneCopy &&
        !(k==state->tfBoundIndexes[j])) {
      posdiff = allBindingSites[k].leftEdgePos - allBindingSites[state->tfBoundIndexes[j]].leftEdgePos;
      //printf("diff=%d\n", posdiff);
      if (abs(posdiff) < HIND_LENGTH) {/*Phey*/
        LOG_ERROR("t=%g steric hindrance breached with site %d (at pos %d, strand %d copy %d gene %d), %d away from site %d \
                   (at pos %d, strand %d copy %d of gene %d)\n",
                  t, k, allBindingSites[k].leftEdgePos, allBindingSites[k].strand, allBindingSites[k].geneCopy, allBindingSites[k].cisregID,
                  posdiff, state->tfBoundIndexes[j], allBindingSites[state->tfBoundIndexes[j]].leftEdgePos, 
                  allBindingSites[state->tfBoundIndexes[j]].strand, allBindingSites[state->tfBoundIndexes[j]].geneCopy, 
                  allBindingSites[state->tfBoundIndexes[j]].cisregID);
        exit(-1);
      }
      if (abs(posdiff) < cooperative_distance) {
        if (posdiff>0) front++; else back++;
      }
    }
    if ((front) && (back)) 
      j=state->tfBoundCount;
  }
  if (front>0) Gibbs -= cooperativity*state->RTlnKr/3;
  if (back>0) Gibbs -= cooperativity*state->RTlnKr/3;
  *koff = NumSitesInGenome*kon*0.25/exp(-Gibbs/(GasConstant*state->temperature));
  /*  fprintf(fperrors,"state->RTlnKr=%g front=%d back=%d H=%d Gibbs=%g koff=%g\n",state->RTlnKr,front,back,allBindingSites[k][4],Gibbs,*koff); */
  /* 25% protein in nucleus is implicit in formula above */
}

/* when TF binding changes, adjust cooperativity at neighbouring sites */
void scan_nearby_sites(int indexChanged,
                       AllTFBindingSites *allBindingSites,
                       CellState *state,
                       GillespieRates *rates,
                       float *koffvalues,
                       float t)
{
  int posdiff, j;
  float diff;
  
  /* for all the currently bound sites */
  for (j = 0; j < state->tfBoundCount; j++) {
    
    /* we are on the same cisreg sequence and we aren't the same TF binding  */
    if (allBindingSites[indexChanged].cisregID == allBindingSites[state->tfBoundIndexes[j]].cisregID && 
        allBindingSites[indexChanged].geneCopy == allBindingSites[state->tfBoundIndexes[j]].geneCopy && 
        !(indexChanged==state->tfBoundIndexes[j])) {

      /* how close are we on the sequence */
      posdiff = allBindingSites[indexChanged].leftEdgePos - allBindingSites[state->tfBoundIndexes[j]].leftEdgePos;
      //printf("diff3=%d\n", posdiff);
      if (abs(posdiff) < HIND_LENGTH) { /* within HIND_LENGTH: bad: shouldn't happen Phey*/
        LOG_ERROR("t=%g steric hindrance 2 has been breached with site %d %d away from site %d\n",
                  t, indexChanged, posdiff, state->tfBoundIndexes[j]);
      }
      if (abs(posdiff) < cooperative_distance) {  /* within 20, adjust koff */

        /* save old value */
        diff = -koffvalues[j];

        /* recompute koffvalues */
        calc_koff(state->tfBoundIndexes[j], allBindingSites, state, &(koffvalues[j]), t);

        /* calculating how koff changes  */
        diff += koffvalues[j];

        /* adjust rates by difference */
        rates->koff += diff;
        rates->koff_operations++;
      }
    }
  }
}

void remove_kon(int siteID,
                int TFID,
                GillespieRates *rates,
                float salphc,
                KonStates *konStates,
                float proteinConcTFID)
{
  int k;
  
  k = 0;
  
  while (!(konStates->konList[TFID]->available_sites[k] == siteID) && k < konStates->konList[TFID]->site_count) {
    k++;
  }

  LOG_VERBOSE_NOCELLID(">>> remove site %d konList (k=%d of %d total sites for TF %d, grandtotal=%d)\n", 
                       siteID, k, konStates->konList[TFID]->site_count, TFID, konStates->nkon);

  /* make sure that we have enough unoccupied sites left */
  if (k < konStates->konList[TFID]->site_count && k < konStates->nkon) { 
    /* adjust rates */
    rates->salphc -= kon*salphc;
    rates->salphc_operations++;

    rates->maxSalphc -= kon*fmaxf(proteinConcTFID, salphc);
    rates->maxSalphc_operations++;

    rates->minSalphc -= kon*fminf(proteinConcTFID, salphc);
    rates->minSalphc_operations++;

    /* one less site available for binding of total */
    (konStates->nkon)--;

    /* also per gene */
    (konStates->nkonsum[TFID])--;

    (konStates->konList[TFID]->site_count)--;

    /* move the last element end of array into space vacated by site k */
    konStates->konList[TFID]->available_sites[k] = konStates->konList[TFID]->available_sites[konStates->konList[TFID]->site_count];
  } else {
    LOG_VERBOSE_NOCELLID("||| couldn't remove site %d from TF %d in konList (k=%d)\n", siteID, TFID, k);
  }
  /* else do nothing: there is likely a redundancy in steric
     hindrance, hence no site to remove */
}

void add_kon(float proteinConcTFID,
             float salphc,
             int TFID,
             int siteID,
             GillespieRates *rates,
             KonStates *konStates)
{

  /* update rates because new site is now available */
  rates->salphc += kon*salphc;
  rates->salphc_operations++;

  rates->maxSalphc += fmaxf(proteinConcTFID, salphc);
  rates->maxSalphc_operations++;

  rates->minSalphc += fminf(proteinConcTFID, salphc);
  rates->minSalphc_operations++;

  /* add back siteID to pool of available sites */
  konStates->konList[TFID]->available_sites[konStates->konList[TFID]->site_count] = siteID;

  /* one more site available */
  (konStates->konList[TFID]->site_count)++;
  (konStates->nkonsum[TFID])++;
  (konStates->nkon)++;
}

/* tests whether criterion for transcription is met */
int ready_to_transcribe(int geneID,
                        int geneCopy,
                        int *tfBoundIndexes,
                        int tfBoundCount,
                        AllTFBindingSites *allBindingSites,
                        int activating[NGENES][MAX_COPIES],
                        int *on)
{
  int i,off;
  
  *on=off=0;
  for (i=0; i < tfBoundCount; i++) {
    if (geneID==allBindingSites[tfBoundIndexes[i]].cisregID &&
        geneCopy==allBindingSites[tfBoundIndexes[i]].geneCopy)
      {
        if (activating[allBindingSites[tfBoundIndexes[i]].tfID][geneCopy]) (*on)++;
        else off++;
    }
  }
  if ((float)off <= 0.33442*(float)(*on) + 0.31303) 
    return (1);
  else 
    return (0);
}

int is_one_activator(int geneID,
                     int geneCopy,
                     int *tfBoundIndexes,
                     int tfBoundCount,
                     AllTFBindingSites *allBindingSites,
                     int activating[NGENES][MAX_COPIES])
{
  int i;
  
  for (i=0; i < tfBoundCount; i++)
    if (geneID==allBindingSites[tfBoundIndexes[i]].cisregID && 
        geneCopy==allBindingSites[tfBoundIndexes[i]].geneCopy &&
        (activating[allBindingSites[tfBoundIndexes[i]].tfID][geneCopy])) 
      return (1);
  return (0);
}

/* only appropriate if nothing is bound */
void calc_from_state(Genotype *genes,
                     CellState *state,
                     GillespieRates *rates,
                     KonStates *konStates,
                     float transport[],
                     float mRNAdecay[])
/* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit */
{
  int i, j, k;
  float salphc, proteinConcTFID; 
  float protein_decay;

  // FIXME: this updates konStates for non-TFs, even though this isn't used
  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to be zero, use aging term */
    protein_decay = genes->proteindecay[i] > 0 ? genes->proteindecay[i] : protein_aging;
    salphc = (float) (state->mRNACytoCount[i]) * genes->translation[i] / (protein_decay);
    konStates->konvalues[i][KON_DIFF_INDEX] = (state->proteinConc[i] - salphc) / (protein_decay);
    konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
    konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
    konStates->nkonsum[i] = 0;
    konStates->konList[i]->site_count = 0;
    //printf("calc_from_state: prot decay[%d]=%g\n", i, konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]);
  }
  state->tfBoundCount=0;

  rates->koff=0.0;
  rates->koff_operations = 0;

  rates->transport=0.0;
  rates->transport_operations = 0;

  rates->mRNAdecay=0.0;
  rates->mRNAdecay_operations = 0;

  rates->picDisassembly=0.0;
  rates->picDisassembly_operations = 0;

  rates->salphc=0.0;
  rates->salphc_operations=0;

  rates->maxSalphc=0.0;
  rates->maxSalphc_operations=0;

  rates->minSalphc=0.0;
  rates->minSalphc_operations=0;

  for (k=0; k < genes->bindSiteCount; k++) {
    i = genes->allBindingSites[k].tfID;
    proteinConcTFID = state->proteinConc[i];
    salphc = konStates->konvalues[i][KON_SALPHC_INDEX];

    rates->salphc += salphc;
    rates->salphc_operations++;

    rates->maxSalphc += fmaxf(proteinConcTFID, salphc);
    rates->maxSalphc_operations++;

    rates->minSalphc += fminf(proteinConcTFID, salphc);
    rates->minSalphc_operations++;

    /* update the list of sites that bind for a particular TF, i */
    konStates->konList[i]->available_sites[konStates->konList[i]->site_count] = k;
    (konStates->konList[i]->site_count)++;
    (konStates->nkonsum[i])++;
  }

  /* initialize konStates->nkon as the total number of binding sites */
  konStates->nkon = genes->bindSiteCount;

  for (i=0; i < NGENES; i++) {

    // FIXME: this updates konStates for non-TFs, even though this isn't used
    LOG("after initializing konStates for gene=%d, nkonsum=%d, site_count=%d, tfsPerGene=%d, nkon=%d\n", 
        i, konStates->nkonsum[i], konStates->konList[i]->site_count, genes->tfsPerGene[i], konStates->nkon);

    transport[i] = kRNA * (float) (state->mRNANuclearCount[i]);
    rates->transport += transport[i];
    rates->transport_operations++;
    LOG_VERBOSE("initializing transport[%d]=%g\n", i, transport[i]);
  }
  LOG_VERBOSE("initializing rates->transport=%g\n", rates->transport);

  /* start all genes in acteylated state */
  for (j=0; j < MAX_COPIES; j++) {  
    int pos = 0;
    for (i=0; i < NGENES; i++) {
      if (genes->copies[i] > j) {  
        state->statechangeIDs[ACETYLATION][j][pos] = i;
        LOG_VERBOSE("Initializing statechange gene=%d, ploidy=%d statechangeIDs[%d][%d]=%d\n", i, j, j, 
                    pos, state->statechangeIDs[ACETYLATION][j][pos]);
        pos++;
      }
    }
  }

  rates->salphc *= kon;
  rates->salphc_operations++;

  rates->maxSalphc *= kon;
  rates->maxSalphc_operations++;

  rates->minSalphc *= kon;
  rates->minSalphc_operations++;

  /* first initialize everything at zero */
  for (j=0; j < MAX_COPIES; j++) {
    rates->acetylationCount[j]=0;
    rates->acetylationCount_operations=0;

    rates->deacetylationCount[j]=0;
    rates->deacetylationCount_operations=0;

    rates->picAssemblyCount[j]=0;
    rates->picAssemblyCount_operations=0;

    rates->transcriptInitCount[j]=0;
    rates->transcriptInitCount_operations=0;

    rates->picDisassemblyCount[j]=0;
    rates->picDisassemblyCount_operations=0;
  }

  /* now set the per copy acetylation rate  */
  for (i=0; i < NGENES; i++) {
    for (j=0; j < genes->copies[i]; j++) {
      rates->acetylationCount[j]++;
      rates->acetylationCount_operations++;
    }
    
  }

  if (verbose) 
    for (j=0; j < MAX_COPIES; j++) {
      LOG_VERBOSE("rates->acetylationCount[%d]=%d\n", j, rates->acetylationCount[j]);
    }
}

/* returns:
 *  0 if there is no fixed event occuring before time t
 *  1 if a transcription event happens before time t
 *  2 if a translation event happens before time t
 *  3 if a gene replication event happens before time t
 */
int does_fixed_event_end(FixedEvent *mRNATranslTimeEnd,
                         FixedEvent *mRNATranscrTimeEnd,
                         FixedEvent *replicationTimeEnd,
                         float t)
{
  int retval;
  float t1;
  float t2;
  float t3;

  if (mRNATranscrTimeEnd == NULL && mRNATranslTimeEnd==NULL && replicationTimeEnd==NULL) {
    retval = 0;
  } else {
    // TODO: rewrite this to avoid use of magic number
    t1 = mRNATranscrTimeEnd ? mRNATranscrTimeEnd->time : TIME_INFINITY;
    t2 = mRNATranslTimeEnd ? mRNATranslTimeEnd->time : TIME_INFINITY;
    t3 = replicationTimeEnd ? replicationTimeEnd->time : TIME_INFINITY;

    LOG_VERBOSE_NOCELLID("check fixed event: t1=%g, t2=%g, t3=%f [t=%g] ", t1, t2, t3, t);

    if ((t1 < t2) && (t1 < t3) && (t1 < t)) { 
      if (mRNATranscrTimeEnd == NULL) retval = 0;
      else retval = 1;
    } else
      if ((t2 < t3) && (t2 < t1) && (t2 < t)) {
        if (mRNATranslTimeEnd == NULL)  retval = 0;
        else retval = 2;
      }
      else 
      if ((t3 < t1) && (t3 < t2) && (t3 < t)) {
        if (replicationTimeEnd == NULL) retval = 0;
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
             int mRNACytoCount[],
             int mRNATranslCytoCount[],
             int cellID)
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
     cytoplasm (mRNACytoCount) and ones that have only just recently arrived
     (mRNATranslCytoCount) */
  for (i=0; i < NGENES; i++) {
    mRNAdecay[i] = mRNAdecayrates[i] * ((float) mRNACytoCount[i] + (float) mRNATranslCytoCount[i]);
    rates->mRNAdecay += mRNAdecay[i];
    rates->mRNAdecay_operations++;
  }

  /* recompute and cache the total rate in data structure */
  rates->total += rates->koff;
  rates->total += rates->transport;
  rates->total += rates->mRNAdecay;
  rates->total += rates->picDisassembly;
  rates->total += rates->salphc;

  /* 
   * convert the counts back into rates using the constants 
   */
  for (j=0; j < MAX_COPIES; j++) {
    rates->total += (float) rates->acetylationCount[j] * acetylate;
    rates->total += (float) rates->deacetylationCount[j] * deacetylate;
    rates->total += (float) rates->picAssemblyCount[j] * PICassembly;
    rates->total += (float) rates->transcriptInitCount[j] * transcriptinit;    
  } 

  tbound1 = *x/(rates->total + rates->maxSalphc);
  tbound2 = *x/(rates->total + rates->minSalphc);
  LOG_VERBOSE_NOCELLID("[cell %03d] bounds %g %g\n", cellID, tbound1, tbound2);

  /* if bounds are the same, simply choose tbound1 */
  if (tbound1==tbound2){
    if (konStates->nkon!=0) {
      LOG_ERROR_NOCELLID("[cell %03d] nkon=%d when it should be zero x=%f rates->maxSalphc=%g rates->minSalphc=%g rates->total=%g\n",
                         cellID, konStates->nkon, *x, rates->maxSalphc, rates->minSalphc, rates->total);
    }
    *dt = tbound1;
  } else {
    /* otherwise get delta t by solving the equation using Newton-Raphson method */
    *dt = rtsafe(&calc_time, *x, rates, konStates, tbound1, tbound2, (float) RT_SAFE_EPSILON); 
  }
}

void end_transcription(float *dt,
                       float t,
                       CellState *state,
                       float transport[NGENES],
                       GillespieRates *rates)
{
  int i, j, total;
  
  /* recompute the delta-t based on difference between now and the
     time of transcription end */
  *dt = state->mRNATranscrTimeEnd->time - t;

  if (verbose) {
    total = 0;
    for (i=0; i < NGENES; i++) 
      for (j=0; j < MAX_COPIES; j++) 
        total += state->mRNATranscrCount[i][j];
    
    LOG_VERBOSE("\ntranscription event finishes out of %d possible t=%g dt=%g\n", total, t, *dt);
  }

  /* get the gene which is ending transcription */
  i = state->mRNATranscrTimeEnd->geneID;
  j = state->mRNATranscrTimeEnd->copy;

  /* increase number of mRNAs in nucleus */
  (state->mRNANuclearCount[i])++;

  /* decrease the number of mRNAs undergoing transcription */
  (state->mRNATranscrCount[i][j])--;

  /* delete the fixed even which has just occurred */
  delete_fixed_event_start(&(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));

  /* add rate kRNA to transport and Gillespie rates */
  transport[i] += kRNA;
  rates->transport += kRNA;
  rates->transport_operations++;

  LOG_VERBOSE("add one new mRNA in nucleus, updating transport[%d]=%g, rates->transport=%g\n", i, transport[i], rates->transport);

}


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
                     int geneID,
                     int geneCopy,
                     GillespieRates *rates)
{
  float disassembly = genes->PICdisassembly[geneID][geneCopy];
  remove_from_array(geneID, TRANSCRIPTINIT, state->statechangeIDs[TRANSCRIPTINIT][geneCopy], &(rates->transcriptInitCount[geneCopy]), (int) 1);
  remove_from_array(geneID, PICDISASSEMBLY, state->statechangeIDs[PICDISASSEMBLY][geneCopy], &(rates->picDisassemblyCount[geneCopy]), (int) 1);
  rates->picDisassembly -= disassembly;
  rates->picDisassembly_operations++;
  
  /* disassemble PIC in OFF state */
  if (state->active[geneID][geneCopy] == OFF_PIC) {
    (state->active[geneID][geneCopy]) = OFF_NO_PIC;
    state->statechangeIDs[DEACETYLATION][geneCopy][rates->deacetylationCount[geneCopy]] = geneID;
    (rates->deacetylationCount[geneCopy])++;
    rates->deacetylationCount_operations++;
  }
  /* disassemble PIC in ON state */
  if (state->active[geneID][geneCopy] == ON_FULL) {
    (state->active[geneID][geneCopy]) = ON_NO_PIC;
  }
}

void revise_activity_state(int geneID,
                           int geneCopy,
                           Genotype *genes,
                           CellState *state,
                           GillespieRates *rates)
{
  int transcriptrule, oldstate, numactive;

  transcriptrule = ready_to_transcribe(geneID, geneCopy, 
                                       state->tfBoundIndexes, 
                                       state->tfBoundCount,
                                       genes->allBindingSites,
                                       genes->activating,
                                       &numactive);

  /* get last state of transcription initiation */
  oldstate = state->active[geneID][geneCopy];

  /* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit, 4-PIC disassembly*/

  /*
   * first set of rules:
   * ACTIVATING TFs exceed REPRESSING TFs 
   */

  /* OFF_FULL -> ON_WITH_NUCLEOSOME */
  if ((transcriptrule) && oldstate==OFF_FULL){
    state->active[geneID][geneCopy] = ON_WITH_NUCLEOSOME;
    state->statechangeIDs[ACETYLATION][geneCopy][rates->acetylationCount[geneCopy]] = geneID;
    (rates->acetylationCount[geneCopy])++;
    rates->acetylationCount_operations++;
  }
  /* OFF_NO_PIC -> ON_NO_PIC */
  if ((transcriptrule) && oldstate==OFF_NO_PIC) {
    state->active[geneID][geneCopy] = ON_NO_PIC;
    remove_from_array(geneID, DEACETYLATION,  state->statechangeIDs[DEACETYLATION][geneCopy], &(rates->deacetylationCount[geneCopy]), (int) 1);
    if (numactive){
      state->statechangeIDs[PICASSEMBLY][geneCopy][rates->picAssemblyCount[geneCopy]] = geneID;
      (rates->picAssemblyCount[geneCopy])++;
      rates->picAssemblyCount_operations++;
    }
  }
  /* OFF_PIC -> ON_FULL */
  if ((transcriptrule) && oldstate==OFF_PIC) {
    state->active[geneID][geneCopy] = ON_FULL;
  }

  /*
   * second set of rules:
   * REPRESSING TFs exceed ACTIVATING TFs 
   */

  /* ON_WITH_NUCLEOSOME -> OFF_FULL */
  if (!(transcriptrule) && oldstate==ON_WITH_NUCLEOSOME) {
    state->active[geneID][geneCopy] = OFF_FULL;
    LOG_VERBOSE("removing gene=%d, copy=%d from statechangeIDs[ACETYLATION][%d]\n", geneID, geneCopy, geneCopy);
    remove_from_array(geneID, ACETYLATION, state->statechangeIDs[ACETYLATION][geneCopy], &(rates->acetylationCount[geneCopy]), (int) 1);
  }
  
  /* ON_NO_PIC -> OFF_NO_PIC */
  if (!(transcriptrule) && oldstate==ON_NO_PIC){          
    state->active[geneID][geneCopy] = OFF_NO_PIC;
    remove_from_array(geneID, PICASSEMBLY, state->statechangeIDs[PICASSEMBLY][geneCopy], &(rates->picAssemblyCount[geneCopy]), (int) 0);
    state->statechangeIDs[DEACETYLATION][geneCopy][rates->deacetylationCount[geneCopy]] = geneID;
    (rates->deacetylationCount[geneCopy])++;
    rates->deacetylationCount_operations++;
  }
  /* ON_FULL -> OFF_PIC  */
  if (!(transcriptrule) && oldstate==ON_FULL) {
    state->active[geneID][geneCopy] = OFF_PIC;
  }

  /* do remaining transitions:
   * OFF_PIC -> OFF_NO_PIC
   * ON_FULL -> ON_NO_PIC 
   */
  if ((state->active[geneID][geneCopy]==OFF_PIC || state->active[geneID][geneCopy]==ON_FULL) && numactive==0)
    disassemble_PIC(state, genes, geneID, geneCopy, rates);

  if (verbose && (oldstate!=state->active[geneID][geneCopy])) {
    LOG_VERBOSE("state change from %d to %d in gene %d, copy %d\n", oldstate, state->active[geneID][geneCopy], geneID, geneCopy);
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
  int i, j, k, bound, siteID, geneID, geneCopy;

  i = 0;

  /* given site 'site', look for the index in the list of bound sites */
  while ((state->tfBoundIndexes[i] != site) && (i < state->tfBoundCount)) 
    i++;
  if (i == state->tfBoundCount) {  /* couldn't find the site */
    LOG_ERROR("t=%g could not find site %d with %d possibilities\n Bound sites are\n",
              t, site, state->tfBoundCount);
    for (j = 0; j < state->tfBoundCount; j++)  {
      LOG_NOFUNC("%d\n", state->tfBoundIndexes[j]);
    }
  }
  else {
    j = 0;
    /* loop through the sterically hindered sites */
    while (j < state->tfHinderedCount) {

      /* check all sites hindered by binding to location 'site' */
      if (state->tfHinderedIndexes[j][1] == site) {
        k = bound = 0;

        /* is anything else hindering the same site? */
        while (bound == 0 && k < state->tfHinderedCount) {
          if (state->tfHinderedIndexes[j][0] == state->tfHinderedIndexes[k][0] && j != k) 
            bound=1;
          k++;
        }

        /* if nothing else is hindering this site then allow site to be (re-)bound */
        if (bound==0) {
          siteID = state->tfHinderedIndexes[j][0];
          LOG_VERBOSE("Site %d leftEdgePos %d on gene %d freed from steric hindrance\n",
                      siteID, genes->allBindingSites[siteID].leftEdgePos, genes->allBindingSites[siteID].cisregID);

          /* adjust rates by returning kon to pool */
          add_kon(state->proteinConc[genes->allBindingSites[siteID].tfID],
                  konStates->konvalues[genes->allBindingSites[siteID].tfID][KON_SALPHC_INDEX],
                  genes->allBindingSites[siteID].tfID,
                  siteID,
                  rates,
                  konStates);
        }

        /* now we have one less sterically hindered site */
        (state->tfHinderedCount)--;

        if (j < state->tfHinderedCount) {
          /* shorten array by moving the end of the array to the hole opened up by removed site */
          state->tfHinderedIndexes[j][0] = state->tfHinderedIndexes[state->tfHinderedCount][0];
          state->tfHinderedIndexes[j][1] = state->tfHinderedIndexes[state->tfHinderedCount][1];
        }
      } else { /* only increment if we haven't shortened array */
        j++;
      }
    }    

    /* reduce the koff rate by the amount */
    rates->koff -= koffvalues[i];
    rates->koff_operations++;

    /* one less bound site */
    (state->tfBoundCount)--;

    /* shift end of array to hole opened up */
    state->tfBoundIndexes[i] = state->tfBoundIndexes[state->tfBoundCount];

    /* likewise with koffvalues */
    koffvalues[i] = koffvalues[state->tfBoundCount];

    /* find the gene and copy whose cisreg region has an unbinding event */
    geneID = genes->allBindingSites[site].cisregID;
    geneCopy = genes->allBindingSites[site].geneCopy;
    LOG_VERBOSE("Add site %d at leftEdgePos %d on gene %d copy %d freed by unbinding\n",
                site, genes->allBindingSites[site].leftEdgePos, geneID, geneCopy);

    /* adjust kon */
    add_kon(state->proteinConc[genes->allBindingSites[site].tfID],
            konStates->konvalues[genes->allBindingSites[site].tfID][KON_SALPHC_INDEX],
            genes->allBindingSites[site].tfID,
            site,
            rates,
            konStates);

    /* adjust the state of the gene */
    revise_activity_state(geneID, geneCopy, genes, state, rates);

    /* when TF unbinds adjust the co-operativity at close sites */
    scan_nearby_sites(site, genes->allBindingSites, state, rates, koffvalues, t);
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
  int geneID, geneCopy, k, posdiff;

  LOG_VERBOSE("kon1 event at site %d out of %d possible, %d TFs previously bound bindSiteCount=%d\n",
              site, konStates->nkon, state->tfBoundCount, genes->bindSiteCount);

  /* if we have run out of space, double memory  */
  if (state->tfBoundCount >= *maxbound2){
    (*maxbound2) *= 2;

    state->tfBoundIndexes = realloc(state->tfBoundIndexes, (*maxbound2)*sizeof(int));

    /* do the copy */
    *koffvalues = realloc(*koffvalues, (*maxbound2)*sizeof(float));

    /* check return value */
    if (!state->tfBoundIndexes || !(*koffvalues)) {
      LOG_ERROR("memory allocation error resetting maxbound2=%d\n", *maxbound2);
      exit(1);
    }
  }

  /* append the site to end of indexes */
  state->tfBoundIndexes[state->tfBoundCount] = site;
  LOG_VERBOSE("remove site %3d on gene %2d (copy %d)\n", 
              site, genes->allBindingSites[site].cisregID, genes->allBindingSites[site].geneCopy);

  /* remove the site from the kon pool */
  remove_kon(site,
             genes->allBindingSites[site].tfID,
             rates, 
             konStates->konvalues[genes->allBindingSites[site].tfID][KON_SALPHC_INDEX],
             konStates,
             state->proteinConc[genes->allBindingSites[site].tfID]);

  /* recompute the koffvalues */
  calc_koff(site, genes->allBindingSites, state, &((*koffvalues)[state->tfBoundCount]), t);

  LOG_VERBOSE("new koff = %g is number %d\n",
              (*koffvalues)[state->tfBoundCount], (state->tfBoundCount+1));

  /* adjust rates by adding the new koffvalue to rates->koff */
  rates->koff += (*koffvalues)[state->tfBoundCount];
  rates->koff_operations++;

  state->tfBoundIndexes[state->tfBoundCount] = site;

  /* increment number of bound TFs */
  (state->tfBoundCount)++;

  /* adjust co-operative binding in context of new TF */
  scan_nearby_sites(site, genes->allBindingSites, state, rates, *koffvalues, t);

  /* get the gene that the TF is binding to */
  geneID = genes->allBindingSites[site].cisregID;

  /* get the copy that the TF is binding to */
  geneCopy = genes->allBindingSites[site].geneCopy;
  
  /* update steric hindrance data structures */
  /* JM: this cycles over all sites, not just bound ones, in order to
     record redundancy in steric hindrance*/
  for (k = 0; k < genes->bindSiteCount; k++) {

    /* if we are on the same gene and not the same binding site */
    if (geneID == genes->allBindingSites[k].cisregID &&
        geneCopy == genes->allBindingSites[k].geneCopy &&
        !(k==site)) {

      /* check distance from current binding site (k) to the original (site) */
      //fprintf(fperrors, "site=%d k=%d\n", genes->allBindingSites[site].leftEdgePos, genes->allBindingSites[k].leftEdgePos);
      posdiff = genes->allBindingSites[site].leftEdgePos - genes->allBindingSites[k].leftEdgePos;

      /* if within HIND_LENGTH, we prevent binding by adding to steric hindrance */
      if (abs(posdiff) < HIND_LENGTH) {/* Phey */
        /* if not enough memory, reallocate */
        if (state->tfHinderedCount > *maxbound3 - 1) {
          (*maxbound3) *= 2;
          state->tfHinderedIndexes = realloc(state->tfHinderedIndexes,2*(*maxbound3)*sizeof(int));
        }
        
        /* record hindrance: 'site' blocks 'k' */
        state->tfHinderedIndexes[state->tfHinderedCount][0] = k;
        state->tfHinderedIndexes[state->tfHinderedCount][1] = site;

        /* update list of hindered count */
        (state->tfHinderedCount)++;

        LOG_VERBOSE("%d steric hindrance sites after %d blocks site %d\n", state->tfHinderedCount, site, k);

        /* remove the kon from pool */
        remove_kon(k,
                   genes->allBindingSites[k].tfID,
                   rates,
                   konStates->konvalues[genes->allBindingSites[k].tfID][KON_SALPHC_INDEX],
                   konStates,
                   state->proteinConc[genes->allBindingSites[k].tfID]);
      }
    }
  }
  LOG_VERBOSE("tfBoundCount=%d tfHinderedCount=%d maxbound2=%d maxbound3=%d\n",
              state->tfBoundCount, state->tfHinderedCount, *maxbound2, *maxbound3);

  /* gene activity may change as a result of binding */
  revise_activity_state(geneID, geneCopy, genes, state, rates);
}

/*
 * time course of [TF]s represented as array of TimeCourse lists.
 */
void add_time_points(float time,
                     float proteinConc[NPROTEINS],
                     TimeCourse **timecoursestart,
                     TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i < NPROTEINS; i++)
    add_time_point(time, proteinConc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

void add_integer_time_points(float time,
                             int proteinConc[NPROTEINS],
                             TimeCourse **timecoursestart,
                             TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i < NPROTEINS; i++)
    add_time_point(time, (float) proteinConc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

void reach_s_phase(CellState *state, Genotype *genes, float t) {
  int i;

  /* set the state as being in S phase */
  state->in_s_phase = 1;

  // only add replication events if S phase and G2 phase have non-zero length
  // otherwise just jump to division
  if (time_s_phase + time_g2_phase > 0.0) {
    for (i=0; i < NGENES; i++) {
      LOG("add replication for geneID=%d at t=%g\n", i, t + genes->replication_time[i]);

      /* add the replication event to the queue */
      add_fixed_event(i, -1, t + genes->replication_time[i], 
                      &(state->replicationTimeEnd), &(state->replicationTimeEndLast));
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
                                int cellID) 
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
                       cellID, P, P_next, c, t, t+deltat, s_mRNA);
    LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
                       cellID);
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
  fprintf(fp_growthrate[cellID], "%g %g %g %g %g %g\n", t, instantaneous_growth_rate, *integrated_growth_rate, P, s_mRNA, c);
#endif

  return (instantaneous_growth_rate);
}

/* 
 * need some sort of control in case it declines to essentially zero.
 * Add in discrete, stochastic and/or zero values, but this may create
 * false attractor if time steps are small and rising tide keeps
 * getting rounded down
 */
void update_protein_conc_cell_size(float proteinConc[],
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

  rates->maxSalphc = rates->minSalphc = 0.0;
  for (i=0; i < NPROTEINS; i++) {
    if (i == SELECTION_GENE)  /* if we are looking at the selection gene, record protein concentration before update */
      L = proteinConc[i];

    /* update protein decay rates due to dilution caused by growth */
    adjusted_decay = genes->proteindecay[i] + state->growthRate;

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
    proteinConc[i] = konStates->konvalues[i][KON_SALPHC_INDEX]*ect1 + ect*proteinConc[i];

    konStates->konvalues[i][KON_DIFF_INDEX] = (proteinConc[i] - konStates->konvalues[i][KON_SALPHC_INDEX]) / konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX];
    rates->maxSalphc += ((float) konStates->nkonsum[i]) * fmaxf(proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]);
    rates->minSalphc += ((float) konStates->nkonsum[i]) * fminf(proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]);
    
    if (i == SELECTION_GENE) {  /* now find out the protein concentration at end of dt interval and compute growth rate */
      L_next = proteinConc[i];
      instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genes->translation[SELECTION_GENE], state->mRNACytoCount[SELECTION_GENE], 
                                                            genes->translation, state->mRNACytoCount, 
                                                            L, L_next, t, dt, 
                                                            konStates->konvalues[SELECTION_GENE][KON_PROTEIN_DECAY_INDEX], ect, ect1,
                                                            state->cellID);
      state->cellSize = (state->cellSize)*exp(integrated_growth_rate);
      fprintf(fp_cellsize[state->cellID], "%g %g\n", t, state->cellSize);
    }

  }
  rates->maxSalphc *= kon;
  rates->minSalphc *= kon;
  if ((output) && (*timecourselast)->time < t+dt-0.1) 
    add_time_points(t+dt, otherdata, timecoursestart, timecourselast);

  /* update the growth rate for next timestep */
  state->growthRate = instantaneous_growth_rate;
}

void calc_num_bound(float proteinConc[],
                    int tfBoundCount)
{
  float sum;
  int i;
  
  sum = 0.0;
  for (i=0; i < TFGENES; i++) 
    sum += proteinConc[i];
  /* if this is wrong for random sequences adjust Kr accordingly */
  /* fprintf(fperrors, "%d bound %g expected\n", tfBoundCount, 0.0003*sum);*/

  LOG_VERBOSE_NOCELLID("%d bound %g expected\n", tfBoundCount, (CISREG_LEN*TFGENES*sum)/NumSitesInGenome);
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

void get_gene(int rate_array[MAX_COPIES], int pos, int *geneLoc, int *geneCopy)
{
  int i = 0;
  int total_rate = 0;
  *geneCopy = -1;   /* haven't found the copy yet */

  while (i < MAX_COPIES && *geneCopy < 0) {
    if (pos < (total_rate + rate_array[i])) {
      *geneCopy = i;
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

  update_protein_conc_cell_size(state->proteinConc, state, genes, dt, 
                                rates, konStates, 
                                t, timecoursestart, timecourselast, 
                                state->proteinConc);
  
  i = -1;
  konrate2 = 0.0;  

  /* choose gene product (mRNA) that gets transported to cytoplasm
     based on weighting in transport[] array */
  while (i < NGENES && x > konrate2) {
    i++;
    x -= transport[i];
  }

  if (i >= NGENES) {
    LOG_ERROR("[cell %03d] attempted to choose mRNA for gene=%d which doesn't exist\n", state->cellID, i);
    exit(0);
  } 
  
  LOG_VERBOSE("do transport event mRNA from gene=%d from %d copies (x=%g)\n", i, state->mRNANuclearCount[i], x);

  /* one less mRNAs in nucleus */
  (state->mRNANuclearCount[i])--;

  /* it has just arrived in cytoplasm, ready to be translated */
  (state->mRNATranslCytoCount[i])++;
  
  /* add the endtime for translation */
  LOG_VERBOSE("add translation event endtime=%f for mRNA encoded by gene=%d \n", endtime, i);
  //add_fixed_event_end(i, endtime, &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));
  add_fixed_event_end(i, -1, endtime, &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));

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
  int siteID = -1;

  /* loop through all TFs, then choose a particular binding site */
  for (k=0; k < TFGENES; k++) {

    /* if no sites available for this TF, skip to next gene */
    if (konStates->konList[k]->site_count == 0) {
      LOG_VERBOSE("looking at TF: %d, has %d sites available, skipping\n", k, konStates->konList[k]->site_count);
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
    total_konrate2 = ((konStates->konList[k]->site_count)) * konrate2_for_TF;

    LOG_VERBOSE("looking at TF: %d, has %d sites available [konrate2: %g, total_konrate2: %g, x: %g]\n", 
                k, konStates->konList[k]->site_count, konrate2_for_TF, total_konrate2, x); 

    /* if we are already in the appropriate TF, now choose a binding site */
    if (!(x > total_konrate2) || (k == TFGENES - 1)) {
      float konrate2 = 0.0;
      
      LOG_VERBOSE("selecting TF: %d, konrate2: %g, total_konrate2: %g, x: %g\n", k, konrate2_for_TF, total_konrate2, x); 
      
      while (l < (konStates->konList[k]->site_count - 1) && x > konrate2) {
        /* this will record the last binding site before we
           reach the appropriate binding site  */
        l++;
        
        /* get ID of site */
        siteID = konStates->konList[k]->available_sites[l];
        
        LOG_VERBOSE("l: %d, site: %d, binds to TF: %d, x = %g (site_count=%d)\n", l, siteID, k, x, konStates->konList[k]->site_count); 
        
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
  if (siteID == -1) {
    LOG_ERROR("no binding site could be found  TF: total_konrate2: %g, x: %g\n", total_konrate2, x);
  }
  else {
    LOG_VERBOSE("found a binding site l: %d, site: %d, binds to TF: %d, konrate2: %g, x: %g\n", l, siteID, k, konrate2_for_TF, x);  
  }

  if (update_protein_flag) 
    /* update protein concentration before doing the binding */
    update_protein_conc_cell_size(state->proteinConc, state, genes, dt, 
                                  rates, konStates, 
                                  t, timecoursestart, timecourselast,
                                  state->proteinConc);
  
  /* bind siteID, only if found */
  if (siteID != -1)
    attempt_tf_binding(genes, state, rates, &koffvalues, konStates, &maxbound2, &maxbound3, siteID, t);

  /* calculate the number of TFs bound */
  calc_num_bound(state->proteinConc, state->tfBoundCount);
}

void tf_unbinding_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                        KonStates *konStates, float *koffvalues, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                        float konrate, float dt, float t, float x, int update_protein_flag)
{
  int i, j = -1;
  int site;

  if (state->tfBoundCount <= 0) {
    LOG_ERROR("tfBoundCount=%d, t=%g, dt=%g, x=%g, rates->koff=%g\n", state->tfBoundCount, t, dt, x, rates->koff);
    // TODO: disable jump to critical size to stop sim
    // state->cellSize = critical_size; 
    
    for (i=0; i < NGENES; i++) {
      int k;
      LOG_ERROR("[gene %2d]: mRNA nuclear=%d, cyto=%d, transl=%d, transcribing ", 
                i, state->mRNANuclearCount[i], state->mRNACytoCount[i], state->mRNATranslCytoCount[i]);
      for (k=0; k < MAX_COPIES; k++)
        LOG_NOFUNC("copy%d=%d ", k, state->mRNATranscrCount[i][k]);
      LOG_NOFUNC("\n");
    }
    LOG_ERROR("attempting to unbind when nothing bound, recomputing koff and kon rates\n");
    recompute_koff_rates(rates, state, genes, koffvalues, t);    
    recompute_kon_rates(rates, state, genes, konStates, 0);

    // TODO: check that this is correct
    return;
  }

  while (j < state->tfBoundCount && x > 0) {
    j++;
    x -= koffvalues[j];
  }
  if (j==state->tfBoundCount) {
    float konrate2 = 0.0;
    for (i = 0; i < state->tfBoundCount; i++) 
      konrate2 += koffvalues[i];
    LOG_WARNING("t=%g koffvalues add up to %g instead of rates->koff=%g\n",
                t, konrate2, rates->koff);
    rates->koff = konrate2;
    j--; /* a bit of a fudge for rounding error, really should move on to rates->transport, but too complicated for something so minor */
  } 
  site = state->tfBoundIndexes[j];
  LOG_VERBOSE("t=%g koff event %d of %d at site %d\n", t, j, state->tfBoundCount,site);
  if (j < 0) { LOG_ERROR("t=%g (x=%g) koff event %d of %d at site %d\n", t, x, j, state->tfBoundCount,site); exit(-1); }
  
  if (update_protein_flag) 
    /* update protein concentration before removing TF */
    update_protein_conc_cell_size(state->proteinConc, state, genes, dt, 
                                  rates, konStates, 
                                  t, timecoursestart, timecourselast, 
                                  state->proteinConc);
  
  /* remove TF binding from 'site' */
  remove_tf_binding(genes, state, rates, konStates, site, koffvalues, t);
  calc_num_bound(state->proteinConc, state->tfBoundCount);
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
  x = ran1(&seed)*((float) (state->mRNACytoCount[i] + state->mRNATranslCytoCount[i]));

  update_protein_conc_cell_size(state->proteinConc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->proteinConc);
  /* 
   * decay mRNA in cytoplasm 
   */
  if (x < (float)state->mRNACytoCount[i]) {
    LOG_VERBOSE("mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
                i, state->mRNACytoCount[i], state->mRNATranslCytoCount[i]);
    
    /* remove the mRNA from the cytoplasm count */
    (state->mRNACytoCount[i])--;  
    change_mRNA_cytoplasm(i, genes, state, rates, konStates); 
    
  } else {
    /* 
     * decay mRNA in process of translating
     */
    x = ran1(&seed)*((float) state->mRNATranslCytoCount[i]);
    LOG_VERBOSE("mRNA decay event gene %d not from %d copies in cytoplasm but %f from %d copies translating\n",
                i, state->mRNACytoCount[i], trunc(x), state->mRNATranslCytoCount[i]);
    
    /* delete this fixed event: this mRNA will never be translated */
    LOG_VERBOSE("delete fixed TRANSLATION EVENT at time =%d for gene=%d\n", (int) trunc(x), i);
    //delete_fixed_event(i, (int) trunc(x), &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));
    delete_fixed_event(i, -1, (int) trunc(x), &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));
    
    /* remove the mRNA from the count */
    (state->mRNATranslCytoCount[i])--;
    if (verbose) 
      for (j=0; j < NGENES; j++) {
        LOG_VERBOSE("%d copies of gene %d translating\n", state->mRNATranslCytoCount[j], j);
      }
  }
}

void histone_acteylation_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                               KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                               float dt, float t)
{
  int geneLoc, geneCopy;
  float x = ran1(&seed)*((float) sum_rate_counts(rates->acetylationCount));
  
  get_gene(rates->acetylationCount, (int)trunc(x), &geneLoc, &geneCopy);

  /* choose a particular gene to change state */
  int geneID = state->statechangeIDs[ACETYLATION][geneCopy][geneLoc];

  LOG_VERBOSE("acetylation event gene %d (copy %d)\nstate change from %d to 4\n",
              geneID, geneCopy, state->active[geneID][geneCopy]);
  if (state->active[geneID][geneCopy] != ON_WITH_NUCLEOSOME) {
    LOG_ERROR("acetylation event on gene %d (copy %d) attempted from state %d\n", geneLoc, geneCopy, state->active[geneID][geneCopy]);
  }

  /* update protein concentration and cell size */
  update_protein_conc_cell_size(state->proteinConc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->proteinConc);
  
  /* set state: eject nucleosome, but there is no PIC yet */
  state->active[geneID][geneCopy] = ON_NO_PIC;
  remove_from_array(geneID, ACETYLATION, state->statechangeIDs[ACETYLATION][geneCopy], &(rates->acetylationCount[geneCopy]), (int) 1);
  if (is_one_activator(geneID, geneCopy, state->tfBoundIndexes, state->tfBoundCount, 
                       genes->allBindingSites, genes->activating)) {
    state->statechangeIDs[PICASSEMBLY][geneCopy][rates->picAssemblyCount[geneCopy]] = geneID; 
    (rates->picAssemblyCount[geneCopy])++;
    rates->picAssemblyCount_operations++;
  }
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                                 KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                                 float dt, float t)
{
  float x = ran1(&seed)*((float) sum_rate_counts(rates->deacetylationCount));

  int geneCopy; 
  int geneLoc; 

  get_gene(rates->deacetylationCount, (int)trunc(x), &geneLoc, &geneCopy);

  /* choose a particular gene and copy to change state */
  int geneID = state->statechangeIDs[DEACETYLATION][geneCopy][geneLoc];

  LOG_VERBOSE("deacetylation event gene %d (copy %d)\nstate change from %d to 1\n",
              geneID, geneCopy, state->active[geneID][geneCopy]);
  if (state->active[geneID][geneCopy] != OFF_NO_PIC) {
    LOG_ERROR("deacetylation event attempted from state %d\n", state->active[geneID][geneCopy]);
  }

  update_protein_conc_cell_size(state->proteinConc, state, genes, dt, 
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->proteinConc);
  /* set state: nucleosome returns */
  state->active[geneID][geneCopy] = OFF_FULL;
  remove_from_array(geneID, DEACETYLATION, state->statechangeIDs[DEACETYLATION][geneCopy], &(rates->deacetylationCount[geneCopy]), (int) 1);
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                      KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t)
{

  float x = ran1(&seed)*((float) sum_rate_counts(rates->picAssemblyCount));

  int geneCopy; 
  int geneLoc; 

  get_gene(rates->picAssemblyCount, (int)trunc(x), &geneLoc, &geneCopy);

  /* choose a particular gene and copy to change state */
  int geneID = state->statechangeIDs[PICASSEMBLY][geneCopy][geneLoc];

  LOG_VERBOSE("PIC assembly event gene %d copy %d\nstate change from %d to 6\n",
              geneID, geneCopy, state->active[geneID][geneCopy]);

  if (state->active[geneID][geneCopy] != ON_NO_PIC) {
    LOG_ERROR("PIC assembly event attempted from state %d\n", state->active[geneID][geneCopy]);
  }

  update_protein_conc_cell_size(state->proteinConc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->proteinConc);
  
  /* turn gene fully on: ready for transcription */
  state->active[geneID][geneCopy] = ON_FULL;
  remove_from_array(geneID, PICASSEMBLY, state->statechangeIDs[PICASSEMBLY][geneCopy], &(rates->picAssemblyCount[geneCopy]), (int) 1);
  state->statechangeIDs[TRANSCRIPTINIT][geneCopy][rates->transcriptInitCount[geneCopy]] = geneID;

  (rates->transcriptInitCount[geneCopy])++;
  rates->transcriptInitCount_operations++;

  state->statechangeIDs[PICDISASSEMBLY][geneCopy][rates->picDisassemblyCount[geneCopy]] = geneID;

  (rates->picDisassemblyCount[geneCopy])++;
  rates->picDisassemblyCount_operations++;

  rates->picDisassembly += genes->PICdisassembly[geneID][geneCopy];
  rates->picDisassembly_operations++;
}

void disassemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                           KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                           float dt, float t, float x)
{
  int geneCopy, geneLoc, geneID;
  int j=-1;
  while (j < NGENES*current_ploidy && x>0) {
    j++;

    get_gene(rates->picDisassemblyCount, j, &geneLoc, &geneCopy);

    x -= genes->PICdisassembly[state->statechangeIDs[PICDISASSEMBLY][geneCopy][geneLoc]][geneCopy];
  }
  if (j==NGENES*current_ploidy) { LOG_ERROR("error in PIC disassembly\n"); }
  geneID = state->statechangeIDs[PICDISASSEMBLY][geneCopy][geneLoc];
  LOG_VERBOSE("PIC disassembly event in copy %d of gene %d\n", geneCopy, geneID);
  disassemble_PIC(state, genes, geneID, geneCopy, rates);
}

void transcription_init_event(GillespieRates *rates, CellState *state, Genotype *genes,
                              KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                              float dt, float t, float x)
{
  int geneID;

  x /= transcriptinit;

  int geneCopy; 
  int geneLoc; 

  /* choose the gene and copy that gets transcribed */
  get_gene(rates->transcriptInitCount, (int)trunc(x), &geneLoc, &geneCopy);
  geneID = state->statechangeIDs[TRANSCRIPTINIT][geneCopy][geneLoc];
  LOG_VERBOSE("transcription event gene %d, copy %d\n", geneID, geneCopy);

  if (state->active[geneID][geneCopy] != ON_FULL && state->active[geneID][geneCopy] != OFF_PIC) {
    LOG_ERROR("transcription event attempted from state %d\n", state->active[geneID][geneCopy]);
  }

  update_protein_conc_cell_size(state->proteinConc, state, genes, dt,
                                rates, konStates, 
                                t, timecoursestart, timecourselast,
                                state->proteinConc);

  /* now that transcription of gene has been initiated, 
   * we add the time it will end transcription, 
   * which is dt+time of transcription from now */
  //add_fixed_event_end(geneID, t+dt+ttranscription, 
  //                    &(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));
  add_fixed_event_end(geneID, geneCopy, t+dt+ttranscription, 
                      &(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));

  /* increase the number mRNAs being transcribed */
  (state->mRNATranscrCount[geneID][geneCopy])++;                      
}
/* -----------------------------------------------------
 * END
 * Functions that handle each possible Gillespie event 
 * ----------------------------------------------------- */

/* Helper function that shifts siteIDs in KonStates, tfHinderedIndexes
   and tfBoundIndexes after a position by the specified offset */
void shift_binding_site_ids(CellState *state, 
                            KonStates *konStates,
                            int end,
                            int offset)
{
  int i, j, k, siteID;

  /* shift all sites in konList */
  for (i=0; i < TFGENES; i++) {
    k = 0;
    while (k < konStates->konList[i]->site_count) {
      siteID = konStates->konList[i]->available_sites[k];
      if (siteID >= end)
        konStates->konList[i]->available_sites[k] += offset;
      k++;
    }
  }

  /* shift all sites in tfHinderedIndexes */
  j = 0;
  while (j < state->tfHinderedCount) {
    siteID = state->tfHinderedIndexes[j][0];
    if (siteID >= end)
      state->tfHinderedIndexes[j][0] += offset;
    siteID = state->tfHinderedIndexes[j][1];
    if (siteID >= end)
      state->tfHinderedIndexes[j][1] += offset;
    j++;
  }

  /* shift all sites in tfBoundIndexes */
  for (j = 0; j < state->tfBoundCount; j++) {
    siteID = state->tfBoundIndexes[j];
    if (siteID >= end)
      state->tfBoundIndexes[j] += offset;
  }
}

/* eject TFs and replicate DNA */
void replicate_gene(CellState *state,
                    Genotype *genes,
                    GillespieRates *rates,
                    KonStates *konStates,
                    float *koffvalues,
                    int geneID,
                    float t) 
{
  int i, k, l, p;
  int start_tfbs_pos, end_tfbs_pos, number_tfbs_pre_replication, offset;
  
  LOG("[gene %2d] duplicating at t=%g\n", geneID, t);
  
  /* double the number of copies of gene being replicating */
  genes->copies[geneID] = 2*current_ploidy;

  LOG("[gene %2d] before removing all TFs: nkonsum=%d, site_count=%d, tfsPerGene=%d, nkon=%d\n", 
      geneID, konStates->nkonsum[geneID], konStates->konList[geneID]->site_count, genes->tfsPerGene[geneID],
      konStates->nkon);
  LOG("[gene %2d] before removing all TFs: tfBoundCount=%d, tfHinderedCount=%d\n", geneID, state->tfBoundCount, state->tfHinderedCount);
  
  /* first eject all TFs from this particular gene */
  
  int tfCount = state->tfBoundCount;
  for (k=0; k < tfCount; k++) {
    i = state->tfBoundIndexes[k];            /* first get the binding site ID */
    l = genes->allBindingSites[i].cisregID;  /* now get the gene it belongs to */
    p = genes->allBindingSites[i].geneCopy;
    
    // if we looking at the gene in question
    if (l == geneID) {
      /* remove TF binding and update rates */
      LOG_VERBOSE("[gene %2d] ejecting TF on binding siteID=%d on copy=%d\n", i, l, p);
      remove_tf_binding(genes, state, rates, konStates, i, koffvalues, t);
    }
  }

  /* remove all PICs on that gene */
  for (p=0; p < current_ploidy; p++) 
    if ((state->active[geneID][p]==OFF_PIC || state->active[geneID][p]==ON_FULL))
      disassemble_PIC(state, genes, geneID, p, rates);
  
  LOG_VERBOSE("[gene %2d] number of binding sites before adding new sites=%d at t=%g\n", 
              geneID, genes->bindSiteCount, t);
  
  /* do mutation */
  // TODO: make mutation rate a parameter 
  for (p=0; p<2*current_ploidy; p++)
    mutate(genes, geneID, p, 0.01);  

  /* record number of TFBS pre-replication */
  number_tfbs_pre_replication = genes->tfsPerGene[geneID];
  start_tfbs_pos = 0;
  end_tfbs_pos = 0;

  /* record the beginning and end siteIDs of the pre-replication list
     of binding sites  */
  for (i=0; i<=geneID; i++) {  
    start_tfbs_pos = end_tfbs_pos;             
    end_tfbs_pos += genes->tfsPerGene[i];
  }
  end_tfbs_pos--;  /* always one less than the end point*/

  LOG("[gene %2d] has %d TFBS before replication [run from %d to %d]\n", 
      geneID, number_tfbs_pre_replication, start_tfbs_pos, end_tfbs_pos);

  LOG("[gene %2d] after removing all TFs: nkonsum=%d, site_count=%d, tfsPerGene=%d, nkon=%d\n", 
      geneID, konStates->nkonsum[geneID], konStates->konList[geneID]->site_count, genes->tfsPerGene[geneID],
      konStates->nkon);

  LOG("[gene %2d] after removing all TFs: tfBoundCount=%d, tfHinderedCount=%d\n", geneID, state->tfBoundCount, state->tfHinderedCount);  


  /* remove all of these old TFBS from konStates, some of them may no
     longer exist after mutation, we re-add them with add_kon() call
     after new binding sites computed */
  for (k=start_tfbs_pos; k < end_tfbs_pos + 1; k++) {
    remove_kon(k,
               genes->allBindingSites[k].tfID,
               rates, 
               konStates->konvalues[genes->allBindingSites[k].tfID][KON_SALPHC_INDEX],
               konStates,
               state->proteinConc[genes->allBindingSites[k].tfID]);
  }
  
  /* recompute *all* binding sites, then relabel sites offset by
     insertion (or deletion) of new sites created by replication */
  calc_all_binding_sites(genes->copies, 
                         genes->cisRegSeq, 
                         genes->transcriptionFactorSeq, 
                         &(genes->bindSiteCount),
                         &(genes->allBindingSites),
                         genes->hindrancePositions,
                         genes->tfsPerGene,
                         genes->tfsStart); 

  /* print_all_binding_sites(genes->copies, genes->allBindingSites, genes->bindSiteCount, 
     genes->transcriptionFactorSeq, genes->cisRegSeq);  */

  /* use new tfsPerGene and pre-replication number to compute the
     offset to shift the siteIDs */
  offset = genes->tfsPerGene[geneID] - number_tfbs_pre_replication;

  LOG_VERBOSE("[gene %2d] has %d TFBS after replication [run from %d to %d]\n", 
              geneID, genes->tfsPerGene[geneID], start_tfbs_pos, end_tfbs_pos + offset);
  LOG_VERBOSE(" shift all TFBS starting at %d by %d\n", end_tfbs_pos + 1, offset);
  LOG_VERBOSE(" number of binding sites after adding new sites=%d at t=%g\n", genes->bindSiteCount, t);

  /* starting at the original ending point, move all sites along by
     'offset'.  Note this assumes that TFBS for a particular gene are
     always stored contiguously. */
  shift_binding_site_ids(state, konStates, end_tfbs_pos + 1, offset);
  
  /* update the konStates data structure to make available the newly
   * created TF binding sites in the full [start, end+offset] region
   */
  LOG("[gene %2d] adding new unbound sites from=%d to=%d\n", geneID, start_tfbs_pos, (end_tfbs_pos + offset));

  for (k=start_tfbs_pos; k <= end_tfbs_pos + offset; k++) {
    add_kon(state->proteinConc[genes->allBindingSites[k].tfID],
            konStates->konvalues[genes->allBindingSites[k].tfID][KON_SALPHC_INDEX],
            genes->allBindingSites[k].tfID,
            k,
            rates,
            konStates);
  }
  
  for (i=0; i < current_ploidy; i++) {
    p = i + current_ploidy;
    
    /* set acetylation state in new gene copy  */
    state->statechangeIDs[ACETYLATION][p][rates->acetylationCount[p]] = geneID;
    
    /* update the counts for the acetylation */
    rates->acetylationCount[p]++;
    rates->acetylationCount_operations++;
    
    LOG_VERBOSE("[gene %2d] [clone acetylation]: ploidy=%d statechangeIDs[%d][%d]=%d\n", geneID, p, p, 
                rates->acetylationCount[p], state->statechangeIDs[ACETYLATION][p][rates->acetylationCount[p]]);
    LOG_VERBOSE(" rates->acetylationCount[%d]=%d\n", p, rates->acetylationCount[p]);
  }
}

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

  rates->maxSalphc=0.0;
  rates->maxSalphc_operations=0;

  rates->minSalphc=0.0;
  rates->minSalphc_operations=0;

  // go over all the currently unbound and increment salphc
  for (k=0; k < genes->bindSiteCount; k++) {
    int siteID;
    int tfID = genes->allBindingSites[k].tfID;
    
    int notbound = 1;
    j = 0;
    while (j < state->tfBoundCount && notbound) {
      siteID = state->tfBoundIndexes[j];
      if (siteID == k)  // site is not available, don't add to kon 
        notbound = 0;
      j++;
    }
    int nothindered = 1;
    j = 0;
    while (j < state->tfHinderedCount && nothindered) {
      siteID = state->tfHinderedIndexes[j][0];
      if (siteID == k)  // site is not available, don't add to kon 
        nothindered = 0;
      j++;
    }

    if (notbound && nothindered) {

      salphc = konStates->konvalues[tfID][KON_SALPHC_INDEX];

      LOG_VERBOSE("for unbound site %03d [tfID=%02d] adding salphc=%g to rates->salphc=%g\n", k, tfID, salphc, rates->salphc);

      rates->salphc += salphc;
      rates->salphc_operations++;

      rates->maxSalphc += fmaxf(state->proteinConc[tfID], salphc);
      rates->maxSalphc_operations++;

      rates->minSalphc += fminf(state->proteinConc[tfID], salphc);
      rates->minSalphc_operations++;

      // TODO: check
      // only do if recalibrating state of cell from scratch
      if (recalibrate) {
        // update the list of sites that bind for a particular TF, i
        konStates->konList[tfID]->available_sites[konStates->konList[tfID]->site_count] = k;
        (konStates->konList[tfID]->site_count)++;
        (konStates->nkonsum[tfID])++;
        konStates->nkon++;
      }
    } 
  }
  rates->salphc = kon*rates->salphc;
  rates->maxSalphc = kon*rates->maxSalphc;
  rates->minSalphc = kon*rates->minSalphc;

  /* now that it is recomputed, add salphc back to total */
  rates->total += rates->salphc;
  LOG("AFTER rates->total=%g, rates->salphc=%g, percent difference=%g\n",   
      rates->total, rates->salphc, 100.0*(fabs(rates->salphc-orig_rates)/rates->salphc));
}

// TODO: refactor to use this functions in  recalibrate_cell() to avoid code duplication.
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
  
  for (i = 0; i < state->tfBoundCount; i++) {
    int site = state->tfBoundIndexes[i];
    calc_koff(site, genes->allBindingSites, state, &(koffvalues[i]), t);
    temprate += koffvalues[i];
  }

#if 0
  fprintf(fp_koff[state->cellID], "%g %g\n", t, rates->koff - temprate);
#endif

  rates->koff = temprate;
  rates->koff_operations = 0;

  /* add back new koff */
  rates->total += rates->koff;
  LOG("AFTER rates->total=%g, rates->koff=%g, percent difference=%g\n",   
      rates->total, rates->koff, 100.0*(fabs(rates->koff-orig_rates)/rates->koff));
}

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
  //int siteID, tfID;

  /* reset the total rate for current step */
  rates->total=0.0;
  
  /* reset all rates and operations */
  rates->koff=0.0;
  rates->koff_operations=0;

  rates->transport=0.0;
  rates->transport_operations=0;

  rates->mRNAdecay=0.0;
  rates->mRNAdecay_operations=0;

  rates->picDisassembly=0.0;
  rates->picDisassembly_operations=0;

  rates->salphc=0.0;
  rates->salphc_operations=0;

  rates->maxSalphc=0.0;
  rates->maxSalphc_operations=0;

  rates->minSalphc=0.0;
  rates->minSalphc_operations=0;

  /* regenerate konStates */
  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to be zero, use aging term */
    protein_decay = genes->proteindecay[i] > 0.0 ? genes->proteindecay[i] : protein_aging;
    salphc = (float) (state->mRNACytoCount[i]) * genes->translation[i] / (protein_decay);
    konStates->konvalues[i][KON_DIFF_INDEX] = (state->proteinConc[i] - salphc) / (protein_decay);
    konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
    konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
    konStates->nkonsum[i] = 0;  
    konStates->konList[i]->site_count = 0;
  }
  konStates->nkon = 0;   /* initialize to zero */

  /* regenerate konStates and rates->{salphc,maxSalphc,minSalphc} */
  recompute_kon_rates(rates, state, genes, konStates, 1);

  for (i=0; i < NGENES; i++) {
    /* transport rates */
    transport[i] = kRNA * (float) (state->mRNANuclearCount[i]);
    rates->transport += transport[i];
    rates->transport_operations++;

    /* regenerate decay rates */
    mRNAdecay[i] = genes->mRNAdecay[i] * ((float) state->mRNACytoCount[i] + (float) state->mRNATranslCytoCount[i]);
    rates->mRNAdecay += mRNAdecay[i];
    rates->mRNAdecay_operations++;

  }

  /* recompute koffvalues for all sites */
  recompute_koff_rates(rates, state, genes, *koffvalues, TIME_INFINITY);

  /* recompute and cache the total rate in data structure */
  rates->total += rates->transport;
  rates->total += rates->mRNAdecay;
  rates->total += rates->picDisassembly;

  /* 
   * convert the counts back into rates using the constants 
   */
  for (j=0; j < MAX_COPIES; j++) {
    rates->total += (float) rates->acetylationCount[j] * acetylate;
    rates->total += (float) rates->deacetylationCount[j] * deacetylate;
    rates->total += (float) rates->picAssemblyCount[j] * PICassembly;
    rates->total += (float) rates->transcriptInitCount[j] * transcriptinit;    
  } 
}



/*
 * recompute rates from scratch to avoid compounding rounding error
 */
// TODO: cleanup to call recompute_*_rates calls, but without requiring
// a complete recalibrate_cell() call
void recompute_rates(GillespieRates *rates,
                     CellState *state,
                     Genotype *genes,
                     KonStates *konStates,
                     float *koffvalues) 
{
  int i, j;
  float protein_decay;
  float salphc, new_salphc = 0.0;
  float new_maxSalphc = 0.0;
  float new_minSalphc = 0.0;
  float new_transport = 0.0;

  float konrate2;

  /* loop through all proteins */
  for (i=0; i < NPROTEINS; i++) {
    /* look at all unoccupied sites to regenerate salphc */
    for (j=0; j < konStates->konList[i]->site_count; j++) {
      //geneID = genes->allBindingSites[state->tfBoundIndexes[j]].tfID;
      /* if protein decay is otherwise going to be zero, use aging term */
      protein_decay = genes->proteindecay[i] > 0 ? genes->proteindecay[i] : protein_aging;
      salphc = ((float) (state->mRNACytoCount[i]) * genes->translation[i] / protein_decay);;
      new_salphc += salphc;
      new_maxSalphc += fmaxf(state->proteinConc[i], salphc);
      new_minSalphc += fminf(state->proteinConc[i], salphc);
    }

    /* transport rates */
    new_transport += kRNA * (float) (state->mRNANuclearCount[i]);
  }

  /* check for rounding error */
  //if (rates->koff < 0.0){
  konrate2 = 0.0;
  for (i=0; i < state->tfBoundCount; i++) konrate2 += koffvalues[i];


  /*if (rates->salphc < 0.0 || 
      rates->maxSalphc < 0.0 ||
      rates->minSalphc < 0.0 ||
      rates->transport < 0.0 ||
      rates->koff < 0.0) { */
    printf("salphc:     old=%g new=%g\n", rates->salphc, new_salphc*kon);
    printf("maxSalphc:  old=%g new=%g\n", rates->maxSalphc, new_maxSalphc*kon);
    printf("minSalphc:  old=%g new=%g\n", rates->minSalphc, new_minSalphc*kon);
    printf("transport:  old=%g new=%g\n", rates->transport, new_transport);
    printf("koffvalues: old=%g new=%g\n", rates->koff, konrate2);
    printf("\n");
    /*}*/

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
  int siteID, newSiteID, geneID, geneCopy;
  int start = from_genotype->tfsStart[gene][from_copy][0];
  int end = from_genotype->tfsStart[gene][from_copy][1];

  /* shift all sites in tfBoundIndexes */
  for (k = 0; k < from_state->tfBoundCount; k++) {
    siteID = from_state->tfBoundIndexes[k];
    geneID = from_genotype->allBindingSites[siteID].cisregID;
    geneCopy = from_genotype->allBindingSites[siteID].geneCopy;

    /* check to see if TF is within the gene copy to be moved */
    if (siteID >= start && siteID <= end && geneCopy == from_copy && geneID == gene) {
      /* fix  */
      if (start  == lastpos + 1)  /* if no offset required: keep siteID */
        newSiteID = siteID;
      else                        /* otherwise offset them by the difference */
        newSiteID = siteID  - (start - lastpos);
      LOG_VERBOSE_NOCELLID("gene=%d [copy=%d] orig siteID=%d, move to new [copy=%d] siteID=%d (new tfBoundCount=%d), lastpos=%d, offset=%d\n", 
                           gene, from_copy, siteID, to_copy, newSiteID, to_state->tfBoundCount, lastpos, (start - lastpos));
      to_state->tfBoundIndexes[to_state->tfBoundCount] = newSiteID;
      to_state->tfBoundCount++;
    }
  }

  /* shift TF hindered indexes */
  for (k = 0; k < from_state->tfHinderedCount; k++) {
    siteID = from_state->tfHinderedIndexes[k][0];
    geneID = from_genotype->allBindingSites[siteID].cisregID;
    geneCopy = from_genotype->allBindingSites[siteID].geneCopy;
    if (siteID >= start && siteID <= end && geneCopy == from_copy && geneID == gene) {
      newSiteID = siteID - (start - lastpos);
      to_state->tfHinderedIndexes[to_state->tfHinderedCount][0] = newSiteID;
      LOG_VERBOSE_NOCELLID("k=%d, gene=%d [copy=%d] orig hindering siteID=%d, move to new [copy=%d] siteID=%d (new tfHinderedCount=%d), lastpos= %d, offset=%d\n", 
                         k, gene, from_copy, siteID, to_copy, newSiteID, to_state->tfHinderedCount, lastpos, (start - lastpos));
      siteID = from_state->tfHinderedIndexes[k][1];
      newSiteID = siteID - (start - lastpos);
      to_state->tfHinderedIndexes[to_state->tfHinderedCount][1] = newSiteID;
      LOG_VERBOSE_NOCELLID("k=%d, gene=%d [copy=%d] orig hindered siteID=%d, move to new [copy=%d] siteID=%d (new tfHinderedCount=%d), lastpos= %d, offset=%d\n", 
                         k, gene, from_copy, siteID, to_copy, newSiteID, to_state->tfHinderedCount, lastpos, (start - lastpos));
      to_state->tfHinderedCount++;
    }
  }

  lastpos += (end - start) + 1;   /* update the last position in new binding site array */
  LOG_VERBOSE_NOCELLID("gene=%d [copy=%d] start=%d, end=%d, lastpos=%d\n", gene, from_copy, start, end, lastpos);
  return lastpos;
}

void clone_queue(FixedEvent **start_orig,
                 FixedEvent **last_orig,
                 FixedEvent **start_clone,
                 FixedEvent **last_clone) 
{
  while (*start_orig != NULL) {
    LOG_NOCELLID("adding geneID=%d time=%g to new clone, removing from orig\n", (*start_orig)->geneID, (*start_orig)->time);
    //add_fixed_event((*start_orig)->geneID, (*start_orig)->time, start_clone, last_clone);
    add_fixed_event((*start_orig)->geneID, (*start_orig)->copy, (*start_orig)->time, start_clone, last_clone);
    delete_fixed_event_start(start_orig, last_orig);
  }

  if (*start_clone != NULL) {
    LOG_NOCELLID("clone queue is not empty: head is: geneID=%d time=%g\n", (*start_clone)->geneID, (*start_clone)->time);
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
      int geneID = (*start_clone)->geneID;
      int copy = (*start_clone)->copy;
      float time = (*start_clone)->time;
      
      if (ran1(&seed) <= fraction) { /* fraction of time move to daughter */
        count_daughter[geneID]++;
        add_fixed_event(geneID, copy, time, start_daughter, last_daughter);
        LOG_NOCELLID("move event at time=%g on gene %d to daughter=%d of total count=%d\n", 
                     time, geneID, count_daughter[i], count_clone[i]);
      } else {  /* (1-fraction of time move to mother */
        count_mother[geneID]++;
        add_fixed_event(geneID, copy, time, start_mother, last_mother);
        LOG_NOCELLID("move event at time=%g on gene %d to mother=%d of total count=%d\n", 
                     time, geneID, count_daughter[i], count_clone[i]);
      }
      /* remove from original queue */
      LOG_NOCELLID("removing event at time=%g on gene %d from clone of queue\n", time, geneID);
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

  state_clone->founderID = state_orig->founderID;
  state_clone->in_s_phase = state_orig->in_s_phase;
  state_clone->burn_in = state_orig->burn_in;
  state_clone->divisions = state_orig->divisions;

  state_clone->cellSize = state_orig->cellSize;
  state_clone->growthRate =  state_orig->growthRate;

  state_clone->mRNATranscrTimeEnd = NULL;
  state_clone->mRNATranscrTimeEndLast = NULL;
  state_clone->mRNATranslTimeEnd = NULL;
  state_clone->mRNATranslTimeEndLast = NULL;
  
  state_clone->RTlnKr = state_orig->RTlnKr;
  state_clone->temperature = state_orig->temperature;

  /* copy Genotype */
  // TODO: check to see how much of these states should be set for all 
  // genotypes or copied from all genotypes, either way it should be 
  // copied from the parent cell
  for (i=0; i < NGENES; i++) {
    for (k=0; k < CISREG_LEN; k++) {
      for (p=0; p < MAX_COPIES; p++) {
        genes_clone->cisRegSeq[i][p][k] = genes_orig->cisRegSeq[i][p][k];
      }
    }
    for (p=0; p < MAX_COPIES; p++) {
      genes_clone->tfsStart[i][p][0] =  genes_orig->tfsStart[i][p][0];
      genes_clone->tfsStart[i][p][1] =  genes_orig->tfsStart[i][p][1];
      genes_clone->activating[i][p]=  genes_orig->activating[i][p];
      genes_clone->PICdisassembly[i][p]=  genes_orig->PICdisassembly[i][p];
    }
    genes_clone->tfsPerGene[i] = genes_orig->tfsPerGene[i];
    genes_clone->copies[i] = genes_orig->copies[i];
    genes_clone->replication_time[i] =  genes_orig->replication_time[i];
    genes_clone->mRNAdecay[i] =  genes_orig->mRNAdecay[i];
    genes_clone->proteindecay[i] =  genes_orig->proteindecay[i];
    genes_clone->translation[i] =  genes_orig->translation[i];
  }
  // FIXME: split out from NGENES loop above
  for (i=0; i < TFGENES; i++) {
    genes_clone->hindrancePositions[i] = genes_orig->hindrancePositions[i];
    for (k=0; k < TF_ELEMENT_LEN; k++) {
      for (p=0; p < MAX_COPIES; p++) {
        genes_clone->transcriptionFactorSeq[i][p][k] = genes_orig->transcriptionFactorSeq[i][p][k];
      }
    }
  }


  // make a pointer to original all_binding_sites, we don't modify it yet
  genes_clone->allBindingSites =  genes_orig->allBindingSites;

  state_clone->tfBoundCount  = 0;  //state_orig->tfBoundCount;
  state_clone->tfHinderedCount = 0; // state_orig->tfHinderedCount;

  state_clone->tfBoundIndexes = NULL;
  state_clone->tfHinderedIndexes = NULL;

  // copy bound indexes if there are any
  if (state_orig->tfBoundCount > 0) {
    state_clone->tfBoundIndexes = malloc(state_orig->tfBoundCount*sizeof(int));
    for (k=0; k < state_orig->tfBoundCount; k++) {
      state_clone->tfBoundIndexes[k] = state_orig->tfBoundIndexes[k];
      state_clone->tfBoundCount++;
    }
  }

  // copy bound hindrance positions if there are any
  if (state_orig->tfHinderedCount > 0) {
    state_clone->tfHinderedIndexes = malloc(2*state_orig->tfHinderedCount*sizeof(int));
    for (k=0; k < state_orig->tfHinderedCount; k++) {
      state_clone->tfHinderedIndexes[k][0] = state_orig->tfHinderedIndexes[k][0];
      state_clone->tfHinderedIndexes[k][1] = state_orig->tfHinderedIndexes[k][1];
      state_clone->tfHinderedCount++;
    }
  }

  /* more state.... */

  /* copy activation state */

  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      state_clone->active[i][p] = state_orig->active[i][p];
      state_clone->statechangeIDs[ACETYLATION][p][i] = state_orig->statechangeIDs[ACETYLATION][p][i];
      state_clone->statechangeIDs[DEACETYLATION][p][i] = state_orig->statechangeIDs[DEACETYLATION][p][i];
      state_clone->statechangeIDs[PICASSEMBLY][p][i] = state_orig->statechangeIDs[PICASSEMBLY][p][i];
      state_clone->statechangeIDs[TRANSCRIPTINIT][p][i] = state_orig->statechangeIDs[TRANSCRIPTINIT][p][i];
      state_clone->statechangeIDs[PICDISASSEMBLY][p][i] = state_orig->statechangeIDs[PICDISASSEMBLY][p][i];

      rates_clone->acetylationCount[p] = rates_orig->acetylationCount[p];
      rates_clone->deacetylationCount[p] =  rates_orig->deacetylationCount[p];
      rates_clone->picAssemblyCount[p] = rates_orig->picAssemblyCount[p];
      rates_clone->transcriptInitCount[p] = rates_orig->transcriptInitCount[p];
      rates_clone->picDisassemblyCount[p] = rates_orig->picDisassemblyCount[p];
      state_clone->mRNATranscrCount[i][p] = state_orig->mRNATranscrCount[i][p];
    }
    state_clone->proteinConc[i] = state_orig->proteinConc[i];;
    state_clone->mRNACytoCount[i] = state_orig->mRNACytoCount[i];
    state_clone->mRNANuclearCount[i] = state_orig->mRNANuclearCount[i];
    state_clone->mRNATranslCytoCount[i] = state_orig->mRNATranslCytoCount[i];
  }

  /* clone queue and empty original */
  clone_queue(&(state_orig->mRNATranscrTimeEnd), &(state_orig->mRNATranscrTimeEndLast),
              &(state_clone->mRNATranscrTimeEnd), &(state_clone->mRNATranscrTimeEndLast));
  clone_queue(&(state_orig->mRNATranslTimeEnd), &(state_orig->mRNATranslTimeEndLast),
              &(state_clone->mRNATranslTimeEnd), &(state_clone->mRNATranslTimeEndLast));
}

void initialize_new_cell_genotype(Genotype *genes, Genotype genes_clone)
{
  int i, p;

  /* initialize hindrance for all TFGENES */
  for (p=0; p < TFGENES; p++) {
    genes->hindrancePositions[p] = genes_clone.hindrancePositions[p];
  }

  /* initialize the non-cisregulatory parts of the genotype */
  for (i=0; i < NGENES; i++) {
    for (p=0; p < MAX_COPIES; p++) {
      genes->activating[i][p]=  genes_clone.activating[i][p];
      genes->PICdisassembly[i][p]=  genes_clone.PICdisassembly[i][p];
    }
    genes->replication_time[i] =  genes_clone.replication_time[i];
    genes->mRNAdecay[i] =  genes_clone.mRNAdecay[i];
    genes->proteindecay[i] =  genes_clone.proteindecay[i];
    genes->translation[i] =  genes_clone.translation[i];
  }
}

void initialize_new_cell_state(CellState *state, CellState state_clone, 
                               GillespieRates *rates, double fraction)
{
  int i, j;
  // reset pointers
  state->mRNATranscrTimeEnd = NULL;
  state->mRNATranscrTimeEndLast = NULL;
  state->mRNATranslTimeEnd = NULL;
  state->mRNATranslTimeEndLast = NULL;
  state->replicationTimeEnd = NULL;
  state->replicationTimeEndLast = NULL;

  fflush(stdout);

  if (state->mRNATranscrTimeEnd != NULL) {
    printf("after delete_queues  mRNATranscrTimeEnd is not NULL\n");
    printf("time=%g, geneID=%d\n", state->mRNATranscrTimeEnd->time, state->mRNATranscrTimeEnd->geneID);
  }

  if (state->replicationTimeEnd != NULL) {
    printf("after delete_queues  replicationTimeEnd is not NULL\n");
    printf("time=%g, geneID=%d\n", state->replicationTimeEnd->time, state->replicationTimeEnd->geneID);
  }

  /* reset the length of TF data structures to zero, we have to reconstruct them */
  state->tfBoundCount = 0;
  state->tfHinderedCount = 0;

  state->in_s_phase = 0;   /* reset S phase state to 0 */
  // TODO: check to see whether we want burn-in at the birth of each new cell ?
  state->burn_in = 0;      /* only do burn-in at beginning of runs */
  state->division_time = TIME_INFINITY;  /* reset division time for this cell */
  state->divisions = state_clone.divisions; /* copy division counter */

  state->cellSize = state_clone.cellSize*fraction;   /* reset cell size */

  /* TODO: check! growth rate, take instantaneous growth rate just before
     division, this will be updated after first new time step */
  // TODO, check that this is always inherited always from mother cell
  state->growthRate = state_clone.growthRate;

  /* keep thermodynamic state the same */ 
  state->RTlnKr = state_clone.RTlnKr;
  state->temperature = state_clone.temperature;

  /* initialize the rate counts */
  for (i=0; i < MAX_COPIES; i++) {
    for (j=0; j < NGENES; j++)
      state->statechangeIDs[ACETYLATION][i][j] = -1;
    rates->acetylationCount[i]=0;
    rates->acetylationCount_operations=0;

    rates->deacetylationCount[i]=0;
    rates->deacetylationCount_operations=0;

    rates->picAssemblyCount[i]=0;
    rates->picAssemblyCount_operations=0;

    rates->transcriptInitCount[i]=0;
    rates->transcriptInitCount_operations=0;

    rates->picDisassemblyCount[i]=0;
    rates->picDisassemblyCount_operations=0;
  }
}

void initialize_new_cell_statechangeIDs(CellState *state,
                                        CellState state_clone,
                                        int type,
                                        int count[MAX_COPIES],
                                        int clone_count[MAX_COPIES],
                                        int copy1,
                                        int copy2,
                                        int geneID) 
{
  int j;
  for (j=0; j < (clone_count[copy1]); j++) {
    int gene = state_clone.statechangeIDs[type][copy1][j];
    if (gene == geneID) {
      state->statechangeIDs[type][0][count[0]] = gene;
      LOG_VERBOSE(" statechangeIDs[%d][0][%2d]=%2d (copy %2d)\n", type, j, state->statechangeIDs[type][0][count[0]], copy1);
      (count[0])++;
    }
  }

  for (j=0; j < clone_count[copy2]; j++) {
    if (state_clone.statechangeIDs[type][copy2][j] == geneID) {
      state->statechangeIDs[type][1][count[1]] = state_clone.statechangeIDs[type][copy2][j];
      LOG_VERBOSE(" statechangeIDs[%d][1][%2d]=%2d (copy %2d)\n", type, j, state->statechangeIDs[type][1][count[1]], copy2);
      (count[1])++;
    }
  }
  count[2] = 0;
  count[3] = 0;
}

void initialize_new_cell_gene(Genotype *genes, Genotype genes_clone, 
                              CellState *state, CellState *state_clone,
                              GillespieRates *rates, GillespieRates *rates_clone,
                              int geneID, int copy1, int copy2, int *lastpos)
{
  int j, k;

  for (k=0; k < CISREG_LEN; k++) {
    genes->cisRegSeq[geneID][0][k] = genes_clone.cisRegSeq[geneID][copy1][k];
    genes->cisRegSeq[geneID][1][k] = genes_clone.cisRegSeq[geneID][copy2][k];
    genes->cisRegSeq[geneID][2][k] = genes_clone.cisRegSeq[geneID][copy1][k];
    genes->cisRegSeq[geneID][3][k] = genes_clone.cisRegSeq[geneID][copy2][k];
  }

  /* don't update TF if this geneID is not controlling a TF */
  if (geneID < TFGENES) {
    for (k=0; k < TF_ELEMENT_LEN; k++) {
      genes->transcriptionFactorSeq[geneID][0][k] = genes_clone.transcriptionFactorSeq[geneID][copy1][k];
      genes->transcriptionFactorSeq[geneID][1][k] = genes_clone.transcriptionFactorSeq[geneID][copy2][k];
      genes->transcriptionFactorSeq[geneID][2][k] = genes_clone.transcriptionFactorSeq[geneID][copy1][k];
      genes->transcriptionFactorSeq[geneID][3][k] = genes_clone.transcriptionFactorSeq[geneID][copy2][k];
    }
  }


  if (genes_clone.copies[geneID] - 1 >= copy1) {
    *lastpos = move_gene_copy(copy1, 0, geneID, &genes_clone, genes, state_clone, state, *lastpos);
  }

  if (genes_clone.copies[geneID] - 1 >= copy2) {
    *lastpos = move_gene_copy(copy2, 1, geneID, &genes_clone, genes, state_clone, state, *lastpos);
  }

  fflush(stdout);

  /* copy activation state */
  state->active[geneID][0] = state_clone->active[geneID][copy1];
  state->active[geneID][1] = state_clone->active[geneID][copy2];

  state->active[geneID][2] = ON_WITH_NUCLEOSOME;
  state->active[geneID][3] = ON_WITH_NUCLEOSOME;

  /* copy statechangeID state */
  // TODO: ultimately maybe acteylationCounts etc. should be moved to "CellState"
  // since they are not cached values

  /* ACETYLATION */
  initialize_new_cell_statechangeIDs(state, *state_clone, ACETYLATION,
                                     rates->acetylationCount, rates_clone->acetylationCount,
                                     copy1, copy2, geneID);

  /* DEACETYLATION */
  initialize_new_cell_statechangeIDs(state, *state_clone, DEACETYLATION,
                                     rates->deacetylationCount, rates_clone->deacetylationCount,
                                     copy1, copy2, geneID);
  /* PICASSEMBLY */
  initialize_new_cell_statechangeIDs(state, *state_clone, PICASSEMBLY,
                                     rates->picAssemblyCount, rates_clone->picAssemblyCount,
                                     copy1, copy2, geneID);
  /* TRANSCRIPTINIT */
  initialize_new_cell_statechangeIDs(state, *state_clone, TRANSCRIPTINIT,
                                     rates->transcriptInitCount, rates_clone->transcriptInitCount,
                                     copy1, copy2, geneID);
  /* PICDISASSEMBLY */
  initialize_new_cell_statechangeIDs(state, *state_clone, PICDISASSEMBLY,
                                     rates->picDisassemblyCount, rates_clone->picDisassemblyCount,
                                     copy1, copy2, geneID);

  /* initialize the per-gene (and/or copy) values of the mRNA counts
     these will be updated later */
  for (j=0; j < MAX_COPIES; j++)
    state->mRNATranscrCount[geneID][j] = 0;
  state->mRNACytoCount[geneID] = 0;
  state->mRNANuclearCount[geneID] = 0;
  state->mRNATranslCytoCount[geneID] = 0;

  /* for transcribing mRNAs: split up the genetic material, as mRNA is
     actually physically attached to the DNA during the replication and
     doesn't go randomly to one or other of the cells split up
     mRNATranscrCount, along with FixedTime events */

  FixedEvent *timeEnd = state_clone->mRNATranscrTimeEnd;

  while (timeEnd != NULL) {
    for (j=0; j < MAX_COPIES; j++)
      LOG_VERBOSE("mRNATranscrCount[%2d][%d]=%d\n", 
                  geneID, j, state_clone->mRNATranscrCount[geneID][j]);
    
    int geneIDQueue = timeEnd->geneID;
    int copy = timeEnd->copy;
    float time = timeEnd->time;
    
    if (geneIDQueue == geneID && copy == copy1) { 
      state->mRNATranscrCount[geneID][0]++;
      add_fixed_event(geneID, 0, time, &(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEnd));
      printf("move mRNATranscrTime event at time=%g on gene %d copy=%d to copy=0 of total mRNATranscrCount=%d\n", 
             time, geneID, copy1, state_clone->mRNATranscrCount[geneID][copy1]);
      /* remove from original queue */
      printf("removing mRNATranscrTime event at time=%g on gene %2d (copy %d) from clone of queue\n", time, geneID, copy1);
      delete_fixed_event(geneID, copy1, 0, &(state_clone->mRNATranscrTimeEnd), &(state_clone->mRNATranscrTimeEndLast));
      (state_clone->mRNATranscrCount[geneID][copy1])--;
    } else if (geneIDQueue == geneID && copy == copy2) {  
      state->mRNATranscrCount[geneID][1]++;
      add_fixed_event(geneID, 1, time, &(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEnd));
      printf("move mRNATranscrTime event at time=%g on gene %d copy=%d to copy=1 of total mRNATranscrCount=%d\n", 
             time, geneID, copy2, state_clone->mRNATranscrCount[geneID][copy2]);
      /* remove from original queue */
      printf("removing mRNATranscrTime event at time=%g on gene %2d (copy %d) from clone of queue\n", time, geneID, copy2);
      delete_fixed_event(geneID, copy2, 0, &(state_clone->mRNATranscrTimeEnd), &(state_clone->mRNATranscrTimeEndLast));
      (state_clone->mRNATranscrCount[geneID][copy2])--;
    }
    timeEnd = timeEnd->next;
    fflush(stdout);
  }
}

void realloc_cell_memory(CellState *state, float **koffvalues) 
{
  //state->tfBoundIndexes = realloc(state->tfBoundIndexes, 2*(state->tfBoundCount+1)*sizeof(int));
  state->tfBoundIndexes = realloc(state->tfBoundIndexes, 10*maxbound*sizeof(int));

  //state->tfHinderedIndexes = realloc(state->tfHinderedIndexes, 10*state->tfHinderedCount*sizeof(int));
  state->tfHinderedIndexes = realloc(state->tfHinderedIndexes, 100*maxbound*sizeof(int));
  
  // reallocate koffvalues memory
  //*koffvalues = realloc(*koffvalues, 2*(state->tfBoundCount+1)* sizeof(float)); 
  *koffvalues = realloc(*koffvalues, 100*maxbound*sizeof(float)); 

  if (!state->tfBoundIndexes || !state->tfHinderedIndexes) {
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
  int original_bind_count = mother->bindSiteCount;
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
  initialize_new_cell_genotype(daughter, genes_clone);
  // print_genotype(daughter, daughterID);

  /* initialize the state of the new cell */
  /* also reset size to appropriate scale of cell size */
  initialize_new_cell_state(daughter_state, state_clone, daughter_rates, fraction);
  // TODO: fix to use scale according to the current size as per below
  //initialize_new_cell_state(daughter_state, state_clone, daughter_rates, fraction*state_clone.cellSize);

  /* set the founderID of the daughter cell from the mother cell */
  daughter_state->founderID = state_clone.founderID;
  daughter_state->divisions = 0;    /* daughter cell resets divisions */ 

  printf("daughter cell %03d is founded by cell %03d with %2d divisions\n", 
         daughterID, daughter_state->founderID, daughter_state->divisions);
  LOG_NOCELLID("daughter cell %03d is founded by cell %03d with %2d divisions\n", 
               daughterID, daughter_state->founderID, daughter_state->divisions);

  if (no_replace_mother) {
    initialize_new_cell_genotype(mother, genes_clone);
    //print_genotype(mother, motherID);

    initialize_new_cell_state(mother_state, state_clone, mother_rates, (1-fraction));
    // TODO (see above)
    //initialize_new_cell_state(mother_state, state_clone, mother_rates, (1-fraction)*state_clone.cellSize);
    mother_state->divisions++;      /* update divisions in mother */
    printf("mother   cell %03d is founded by cell %03d with %2d divisions\n", 
         motherID, mother_state->founderID, mother_state->divisions);
    LOG_NOCELLID("mother   cell %03d is founded by cell %03d with %2d divisions\n", 
                 motherID, mother_state->founderID, mother_state->divisions);
  }

  LOG_VERBOSE_NOCELLID("[cell %03d] (mother) total tfBoundCount=%d, tfHinderedCount=%d\n", motherID, 
                       state_clone.tfBoundCount, state_clone.tfHinderedCount);

  //print_all_binding_sites(mother->copies, mother->allBindingSites, mother->bindSiteCount, 
  //                        mother->transcriptionFactorSeq, mother->cisRegSeq, mother->tfsStart); 

  if (verbose) {
    LOG_NOCELLID("[cell %03d]: acetylation counts:\n", motherID);
    total = 0;
    for (j=0; j < MAX_COPIES; j++)  {
      LOG_NOCELLID(" before mother acetylationCount[%2d]=%2d\n", j, rates_clone.acetylationCount[j]);
      total += rates_clone.acetylationCount[j];
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

    LOG_NOCELLID("initialize daughter=%2d, geneID=%2d\n", daughterID, i);
    initialize_new_cell_gene(daughter, genes_clone, 
                             daughter_state, &state_clone,
                             daughter_rates, &rates_clone,
                             i, daughter_copy1, daughter_copy2, &lastpos_daughter);

    if (no_replace_mother) {
      LOG_NOCELLID("initialize mother  =%2d, geneID=%2d\n", motherID, i);
      initialize_new_cell_gene(mother, genes_clone, 
                               mother_state, &state_clone,
                               mother_rates, &rates_clone,
                               i, mother_copy1, mother_copy2, &lastpos_mother);
    }
  }

  if (verbose) {
    total = 0;
    for (j=0; j < MAX_COPIES; j++)  {
      LOG_NOCELLID(" after daughter acetylationCount[%2d]=%2d\n", j, daughter_rates->acetylationCount[j]);
      LOG_NOCELLID(" after mother acetylationCount[%2d]=%2d\n", j, mother_rates->acetylationCount[j]);
      total += daughter_rates->acetylationCount[j];
      total += mother_rates->acetylationCount[j];
    }
    LOG_NOCELLID(" after total acetylation=%d\n", total);

    // check that we have emptied the list of transcribing mRNAs
    total = 0;
    for (i=0; i < NGENES; i++) 
      for (j=0; j < MAX_COPIES; j++) 
        total += state_clone.mRNATranscrCount[i][j];
    
    LOG_NOCELLID("mRNATranscrCount=%2d left in state_clone after moving to mother+daughter\n", total);
  }

  realloc_cell_memory(daughter_state, daughter_koffvalues);
  if (no_replace_mother) {
    realloc_cell_memory(mother_state, mother_koffvalues);
  }

  /* recompute *all* binding sites in daughter, then relabel sites */
  calc_all_binding_sites(daughter->copies, 
                         daughter->cisRegSeq, 
                         daughter->transcriptionFactorSeq, 
                         &(daughter->bindSiteCount),
                         &(daughter->allBindingSites),
                         //daughter->hindrancePositions,
                         mother->hindrancePositions,
                         daughter->tfsPerGene,
                         daughter->tfsStart); 

  //print_all_binding_sites(daughter->copies, daughter->allBindingSites, daughter->bindSiteCount, 
  //                        daughter->transcriptionFactorSeq, daughter->cisRegSeq, daughter->tfsStart); 

  if (no_replace_mother) {
    /* recompute *all* binding sites in mother, then relabel sites */
    calc_all_binding_sites(mother->copies, 
                           mother->cisRegSeq, 
                           mother->transcriptionFactorSeq, 
                           &(mother->bindSiteCount),
                           &(mother->allBindingSites),
                           mother->hindrancePositions,
                           mother->tfsPerGene,
                           mother->tfsStart); 
  }

  //print_all_binding_sites(mother->copies, mother->allBindingSites, mother->bindSiteCount, 
  //                        mother->transcriptionFactorSeq, mother->cisRegSeq, mother->tfsStart); 

  if (no_replace_mother) {
    LOG_NOCELLID("original number of binding sites=%d should = (mother=%d + daughter=%d) = %d\n", 
                 original_bind_count, mother->bindSiteCount, daughter->bindSiteCount, mother->bindSiteCount + daughter->bindSiteCount);
    if (original_bind_count != mother->bindSiteCount + daughter->bindSiteCount) {
      LOG_ERROR_NOCELLID("original number of binding sites=%d  != (mother=%d + daughter=%d) = %d\n", 
                         original_bind_count, mother->bindSiteCount, daughter->bindSiteCount, mother->bindSiteCount + daughter->bindSiteCount);
      exit(0);
    }
  }
  
  /* split up the volume of the cell */
  for (i=0; i < NPROTEINS; i++) {

    // first protein
    daughter_state->proteinConc[i] = fraction * state_clone.proteinConc[i];
    if (no_replace_mother) {
      mother_state->proteinConc[i] = (1-fraction) * state_clone.proteinConc[i];
      LOG_VERBOSE_NOCELLID("daughter=%g (%g), mother=%g (%g) = total protein=%g\n", 
                   daughter_state->proteinConc[i], fraction, mother_state->proteinConc[i], (1-fraction), state_clone.proteinConc[i]);
    }
  }

  /* now loop over all transcribing genes */
  for (i=0; i < NGENES; i++) {
    // mRNAs in cytoplasm (not translating)
    daughter_state->mRNACytoCount[i] = rint(fraction * state_clone.mRNACytoCount[i]);
    if (no_replace_mother) {
      mother_state->mRNACytoCount[i] =  state_clone.mRNACytoCount[i] - daughter_state->mRNACytoCount[i];
      LOG_VERBOSE_NOCELLID("daughter=%d, mother=%d of total mRNACytoCount=%d\n", 
                   daughter_state->mRNACytoCount[i], mother_state->mRNACytoCount[i], state_clone.mRNACytoCount[i]);
    }

    // mRNAs in nucleus
    daughter_state->mRNANuclearCount[i] = rint(fraction * mother_state->mRNANuclearCount[i]);
    if (no_replace_mother) {
      mother_state->mRNANuclearCount[i] =  state_clone.mRNANuclearCount[i] - daughter_state->mRNANuclearCount[i];
      LOG_VERBOSE_NOCELLID("daughter=%d, mother=%d of total mRNANuclearCount=%d\n", 
                   daughter_state->mRNANuclearCount[i], mother_state->mRNANuclearCount[i], state_clone.mRNANuclearCount[i]);
    }

    // split up mRNATranslCytoCount, along with FixedTime events
    split_mRNA(&(state_clone.mRNATranslTimeEnd), &(state_clone.mRNATranslTimeEndLast), state_clone.mRNATranslCytoCount,
               &(daughter_state->mRNATranslTimeEnd), &(daughter_state->mRNATranslTimeEndLast), daughter_state->mRNATranslCytoCount,
               &(mother_state->mRNATranslTimeEnd), &(mother_state->mRNATranslTimeEndLast), mother_state->mRNATranslCytoCount,
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
    LOG_VERBOSE_NOCELLID("tfBoundCount=%d (motherID=%d), tfBoundCount=%d (daughterID=%d)\n", 
                         mother_state->tfBoundCount, motherID, daughter_state->tfBoundCount, daughterID);
    LOG_VERBOSE_NOCELLID("tfHinderedCount=%d (motherID=%d), tfHinderedCount=%d (daughterID=%d)\n", 
                         mother_state->tfHinderedCount, motherID, daughter_state->tfHinderedCount, daughterID);
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

  LOG("snapshot at time=%g:\n x=%g, koff=%g = %d (tfBoundCount) * %g (koff/tfBoundCount)\n transport=%g\n decay=%g\n",
      t, x, rates->koff, state->tfBoundCount, rates->koff/(float)state->tfBoundCount, 
      rates->transport, rates->mRNAdecay);
  LOG_NOFUNC(" rates->salphc=%g\n rates->maxSalphc=%g rates->minSalphc=%g\n", rates->salphc, rates->maxSalphc, rates->minSalphc);
  LOG_NOFUNC(" konrate=%g\n", konrate);
  LOG_NOFUNC(" PICdisassembly=%g\n kon=%g = %d * %g\n",
             rates->picDisassembly, rates->salphc+(konrate), konStates->nkon, (rates->salphc+(konrate))/(float)konStates->nkon);
  
  for (p=0; p < MAX_COPIES; p++) {
    LOG_NOFUNC(" acetylation=%g (copy %d)\n deacetylation=%g (copy %d)\n PIC assembly=%g (copy %d)\n transcriptinit=%g (copy %d)\n",
               (float)rates->acetylationCount[p]*acetylate, p, (float)rates->deacetylationCount[p]*deacetylate, p, 
               (float)rates->picAssemblyCount[p]*PICassembly, p, (float)rates->transcriptInitCount[p]*transcriptinit, p);
  }
  LOG_NOFUNC(" total rates=%g=%g+%g\n", rates->total + (konrate), rates->total, konrate);
  LOG_NOFUNC(" total free=%d + total bound=%d + total hindered=%d = total sites=%d\n", 
             konStates->nkon, state->tfBoundCount, state->tfHinderedCount, genes->bindSiteCount);
  for (i = 0; i < TFGENES; i++) {
    nkon += konStates->konList[i]->site_count;
    LOG_NOFUNC(" unoccupied binding sites=%d available for TF=%d \n", konStates->konList[i]->site_count, i);
  }
  LOG_NOFUNC(" nkon recomputed=%d\n", nkon);
  LOG_NOFUNC("\n");
}

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
    printf("recalibrating cell %3d after burn-in!\n", state->cellID);
    LOG("recalibrating cell %3d after burn-in!\n", state->cellID);
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
  if (critical_size > 0.0 && state->cellSize >= critical_size && !state->in_s_phase)  { /* run until checkpoint reached */
    reach_s_phase(state, genes, *t);
    /* current time plus 30 mins for S phase and a further 30 mins of growth after S phase */
    state->division_time = *t + time_s_phase + time_g2_phase;   
    LOG("at t=%g add division time=%g, cellSize=%g\n", *t, state->division_time, state->cellSize);
  }
  
  // TODO: currently disable printing out the TF occupancy
#if 0
  print_tf_occupancy(state, genes->allBindingSites, *t);
  print_rounding(state, rates, *t);
#endif
  
  *x = expdev(&seed);        /* draw random number */

  /* check for rounding error */
  // TODO 2009-03-03: will remove very soon after more testing
  if (rates->koff < 0.0){
    konrate2 = 0.0;
    for (i=0; i < state->tfBoundCount; i++) konrate2 += koffvalues[i];
    if ((verbose) || konrate2>0.0) {
      LOG_WARNING("koffvalues add up to %g rates->koff=%g < 0\n", konrate2, rates->koff);
    }
    rates->koff = konrate2;
  }
  
  /* do first Gillespie step to chose next event */
  calc_dt(x, dt, rates, konStates, mRNAdecay, genes->mRNAdecay,
          state->mRNACytoCount, state->mRNATranslCytoCount, state->cellID);

  if (*dt < 0.0) {
    LOG_ERROR("dt=%g is negative after calc_dt, t=%g\n", *dt, *t);
  }

  LOG_VERBOSE("next stochastic event due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);
  
  if (!(state->mRNATranscrTimeEndLast)) {
    for (i=0; i < NGENES; i++)
      for (j=0; j < MAX_COPIES; j++) {
        LOG_VERBOSE("%d transcription events on gene %2d (copy %2d)\n", state->mRNATranscrCount[i][j], i, j);
      }
  }
  
  /* first check to see if a fixed event occurs in current t->dt window,
     or in tdevelopment if running for a fixed development time */
  fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);

  event = does_fixed_event_end(state->mRNATranslTimeEnd,
                               state->mRNATranscrTimeEnd,
                               state->replicationTimeEnd,
                               fixed_time);
  
  /* while there are either transcription or translation events
     occuring in current t->dt window */
  while (event > 0) {
    *konrate = (*x)/(*dt);
    
    switch (event) {
    case 1:   /* if a transcription event ends */
      end_transcription(dt, *t, state, transport, rates);
      
      update_protein_conc_cell_size(state->proteinConc, state, genes, *dt,
                                    rates, konStates, *t,
                                    timecoursestart, timecourselast,
                                    state->proteinConc);
      break;
    case 2:            /* if a translation event ends */
      *dt = state->mRNATranslTimeEnd->time - *t;         /* make dt window smaller */
      total=0;  /* number of translation events */
      
      /* count current number of mRNAs that have recently arrived in cytoplasm */
      for (i=0; i<NGENES; i++) total += state->mRNATranslCytoCount[i];

      LOG_VERBOSE("translation event finishes out of %d possible t=%g dt=%g\n",
                  total, *t, *dt); /* bug: dt can be negative */
      
      /* get identity of gene that has just finished translating */
      i=state->mRNATranslTimeEnd->geneID;   
      
      /* there is one less mRNA that has just finished translation */
      (state->mRNATranslCytoCount[i])--;   
      
      /* delete the event that just happened */
      LOG_VERBOSE("delete translation event that just happened at time=%g", *t);
      delete_fixed_event_start(&(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));
      
      /* there is one more mRNA that is now in cytoplasm */
      (state->mRNACytoCount[i])++;
      
      /* update protein concentration */
      update_protein_conc_cell_size(state->proteinConc, state, genes, *dt,
                                    rates, konStates, *t,
                                    timecoursestart, timecourselast,
                                    state->proteinConc);
      
      /* the number of mRNAs in cytoplasm affects binding */
      change_mRNA_cytoplasm(i, genes, state, rates, konStates);
      break;
    case 3:  /* replicate gene */
      *dt = state->replicationTimeEnd->time - *t;         /* make dt window smaller */
      
      replicate_gene(state, genes, rates, konStates, koffvalues, state->replicationTimeEnd->geneID, *t);
      
      /* delete the event that just happened */
      delete_fixed_event_start(&(state->replicationTimeEnd), &(state->replicationTimeEndLast));
      
      update_protein_conc_cell_size(state->proteinConc, state, genes, *dt,
                                    rates, konStates, *t,
                                    timecoursestart, timecourselast,
                                    state->proteinConc);
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
            genes->mRNAdecay, state->mRNACytoCount, state->mRNATranslCytoCount, state->cellID);
    
    LOG_VERBOSE("next stochastic event (2) due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

    fixed_time = no_fixed_dev_time ? (*t+*dt) : fminf(*t+*dt, tdevelopment);    

    /* check to see there aren't more fixed events to do */
    event = does_fixed_event_end(state->mRNATranslTimeEnd, 
                                 state->mRNATranscrTimeEnd, 
                                 state->replicationTimeEnd,
                                 fixed_time);
  } 

  /* no remaining fixed events to do in dt, now do stochastic events */
  
  /* if we haven't already reached end of development with last delta-t */
  if (*t+*dt < tdevelopment || no_fixed_dev_time) {
    
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
        // TODO: remove
        // recompute_rates(rates, state, genes, konStates, koffvalues) ;
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
          if (*x < rates->picDisassembly) {
            disassemble_PIC_event(rates, state, genes, konStates, 
                                  timecoursestart,  timecourselast, *dt, *t, *x);
          } else {
            *x -= rates->picDisassembly;
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
              if (*x < (float) sum_rate_counts(rates->acetylationCount) * acetylate) {
                
                histone_acteylation_event(rates, state, genes, konStates, 
                                          timecoursestart, timecourselast, *dt, *t);
              } else {
                
                *x -= (float) sum_rate_counts(rates->acetylationCount) * acetylate;
                /* 
                 * STOCHASTIC EVENT: histone deacetylation
                 */
                if (*x < (float) sum_rate_counts(rates->deacetylationCount) * deacetylate) {
                  
                  histone_deacteylation_event(rates, state, genes, konStates, 
                                              timecoursestart, timecourselast, *dt, *t);
                } else {
                  *x -= (float) sum_rate_counts(rates->deacetylationCount) * deacetylate;
                  /* 
                   * STOCHASTIC EVENT: PIC assembly
                   */
                  if (*x < (float) sum_rate_counts(rates->picAssemblyCount) * PICassembly) {
                    assemble_PIC_event(rates, state, genes, konStates, 
                                       timecoursestart, timecourselast, *dt, *t);
                  } else {
                    *x -= (float) sum_rate_counts(rates->picAssemblyCount) * PICassembly;
                    /* 
                     * STOCHASTIC EVENT: transcription initiation
                     */
                    if (*x < (float) sum_rate_counts(rates->transcriptInitCount) * transcriptinit) {
                      transcription_init_event(rates, state, genes, konStates, 
                                               timecoursestart, timecourselast, *dt, *t, *x);
                    } else {
                      /*
                       * FALLBACK: shouldn't get here, previous
                       * events should be exhaustive
                       */
                      
                      LOG_ERROR("[cell %03d] t=%g no event assigned: x=%g, rates->total+konrate=%g, recalibrate cell\n", 
                                state->cellID, *t, *x, rates->total + *konrate);

                      log_snapshot(rates, state, genes, konStates, &koffvalues, mRNAdecay,
                                   transport, *konrate, *x, *t);
                      recalibrate_cell(rates, state, genes, konStates, &koffvalues,
                                       mRNAdecay, transport, *dt); 
                      log_snapshot(rates, state, genes, konStates, &koffvalues, mRNAdecay,
                                   transport, *konrate, *x, *t);
                      if (*t > 1000.0) {
                        LOG_ERROR("should probably stop cell run, has been running for > 1000 mins\n");
                        // TODO
                        // state->cellSize = critical_size;
                      }
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
    update_protein_conc_cell_size(state->proteinConc, state, genes, *dt,
                                  rates, konStates,
                                  *t, timecoursestart, timecourselast,
                                  state->proteinConc);
    /* advance to end of development (this exits the outer while loop) */
    *t = tdevelopment;
  }
  return 0;
}

/*
 * develop: run the population of cells for a given number of
 * divisions or run cell(s) for a specific length of time
 */
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
  for (i=0; i < NGENES; i++) {
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }


  for (j = 0; j < POP_SIZE; j++) {

    output=1;
    /* for all cells in this replicate initialize all parts of the
       genotype, *except* for the cis-reg sequence using the initial
       genes[0]  */
    initialize_genotype(&genes[j], genes[0], kdis, j);
    /* if genotype is held constant, start varying the seed *after*
       initialize_genotype, so we can run the same genotype with
       slightly different amounts of noise  */
    if (hold_genotype_constant)
      for (curr_seed=0; curr_seed<dummyrun; curr_seed++) 
        ran1(&seed);
   
    initialize_cell(&state[j], j, genes[j].copies, genes[j].mRNAdecay, initmRNA, initProteinConc, burn_in);

    /* print binding sites */
    if (output_binding_sites) 
      print_all_binding_sites(genes[j].copies, genes[j].allBindingSites, genes[j].bindSiteCount, 
                              genes[j].transcriptionFactorSeq, genes[j].cisRegSeq, genes[j].tfsStart); 


    
    /* set cell temperature and value of RTlnKr constant */
    state[j].temperature = temperature;
    state[j].RTlnKr = GasConstant * temperature * log(Kr);

    /* initialize time courses */
    for (i=0; i < NPROTEINS; i++){
      timecoursestart[j][i] = NULL;
      timecourselast[j][i] = NULL;
    } 

    add_time_points((float) 0.0, state[j].proteinConc, timecoursestart[j], timecourselast[j]);

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
    //insert_with_priority(&(time_queue), &(time_queue_end), j, t[j]);
    int ops;
    ops = insert_with_priority_heap(queue, j, t[j]);
    LOG_NOCELLID("[cell=%03d] inserted at time=%g\n", j, t[j]); 
  }

  while ((no_fixed_dev_time && divisions < max_divisions) ||    /* no fixed dev time, run until # divisions reached */
         (!no_fixed_dev_time && t_next < tdevelopment)) {       /* or, if fixed dev time, run until tdevelopment reached */

    /* get the next cell with the smallest t to advance next */
    //t_next = get_next(&time_queue, &time_queue_end, &cell);
    int ops, retval;

    LOG_VERBOSE_NOCELLID("before choosing next cell (queue len reaper=%03d, main=%03d)\n", empty_queue->n, queue->n);

    if (empty_queue->n == POP_SIZE || queue->n == 0) {
      /* if all cells in population are dead, we exit */
      LOG_NOCELLID("all %03d cells in population are dead, main queue len=%03d\n", empty_queue->n, queue->n);
      break;
    } 
    /* otherwise get next event from queue */
    t_next = get_next_heap(queue, &cell, &ops);
    LOG_VERBOSE_NOCELLID("[cell=%03d] get minimum time=%g, size=%g (queue len reaper=%03d, main=%03d)\n", 
                         cell, t_next, state[cell].cellSize, empty_queue->n, queue->n);

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

      /* to get new daughter cell slot, first check list of empty cells */
      if (empty_queue->n > 0) {
        get_next_heap(empty_queue, &daughterID, &ops);  /* use one of the empty cells */
        LOG("[cell %03d] in empty_queue (length=%03d) used as daughter cell\n", daughterID, empty_queue->n);
      } else {
        daughterID = rint((POP_SIZE-1)*ran1(&seed));  /* otherwise choose one other cell randomly */
        /* removing pending event in daughter cell from queue, if different from the mother 
           only remove if the daughter cell wasn't already removed from queue */
        if (daughterID != motherID) {
          delete_element_heap(queue, daughterID);
          LOG_NOCELLID("removing pending event from queue (length=%3d) replaced by daughter cell=%d\n", queue->n, daughterID);
        } 
      }

      printf("[cell %03d] dividing into mother=%03d and daughter=%03d at t=%g, division=%g, total divisions=%d (t_next=%g)\n", 
             cell, motherID, daughterID, t[cell], current_division_time, divisions, t_next);
      LOG_NOCELLID("[cell %03d] dividing into mother=%03d and daughter=%03d at t=%g, division=%g, total divisions=%d (t_next=%g)\n", 
                   cell, motherID, daughterID, t[cell], current_division_time, divisions, t_next);


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
      // ops = insert_with_priority(&(time_queue), &(time_queue_end), cell, t[cell]);
      ops = insert_with_priority_heap(queue, cell, t[cell]);
    }
  }

  /* output the founder information */
  for (j = 0; j < POP_SIZE; j++) {
    printf("cell %03d derived from founder %03d had %2d divisions%s\n", 
           j, state[j].founderID, state[j].divisions, t[j] == TIME_INFINITY ? " [dead]": "");
    LOG_NOCELLID("cell %03d derived from founder %03d had %2d divisions%s\n", 
                 j, state[j].founderID, state[j].divisions, t[j] == TIME_INFINITY ? " [dead]": "");
  }
  // TODO: check output contents of reaper queue
  bh_dump(empty_queue);

  /* cleanup data structures */
  bh_free(queue);
  bh_free(empty_queue);

  for (j = 0; j < POP_SIZE; j++) {
    free(koffvalues[j]);
    for (i=0; i < NPROTEINS; i++) {
      free(konStates[j].konList[i]->available_sites);
      free(konStates[j].konList[i]);
    }
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
