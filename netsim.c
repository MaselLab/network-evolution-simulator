/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster
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

/* local includes */
#include "random.h"
#include "lib.h"
#include "netsim.h"
/* #include "list.h" */

/*
  version27: each founding indiv has a different genotype, CopyGenotype is obsolete
  version28: add selection (wrong sort)
  v29: transcription initiation quenching rule, take selection out
  v30reject: add delays to quenching and unquenching
  v30: redo selection based on v29, buggy
  v30alt: redo selection based on v27, free up memory of CellState and TimeCourse
  v31: based on v30alt, replace stability only selection with selection for specific opt.
  Hopefully this helps avoid evolution from always-on to always-off (not sufficient yet).
  Move meanmRNA and protein initial states from initialize_cell to be universal
  v32: koffIDs and nkoff have the same info as B2 and y2: remove them
  v33: implement transcription initiation scheme with bursting, PIC assembly & H2A.Z
  v34: disassembly rates bootstrapped from constitutively on Garcia-Martinez subset, not TF subset distribution
  v35: output only a subset of timecourse measurements to prevent excessive output
  change NumSitesInGenome value
  v36: reglen changed to 150
*/

/* initialize constants as static const */
static const int maxelements=500*PLOIDY; 
/* start by allocating maxelements when initializing a genotype, double as needed, reduce at end */
static const int maxbound=500*PLOIDY;
static const int PopSize=1;
static const int nmin=4;
static const float tdevelopment=120.0;
static const float kon=1e-4; /* lower value is so things run faster */
/* kon=0.2225 is based on 1 molecule taking 240seconds=4 minutes
   and 89% of the proteins being in the nucleus*/
static const float kRNA=618.0;
static const float ttranslation=1.0;
static const float ttranscription=1.0;
static const float pact=0.62;
static const float transcriptinit=8.5; /* replace betaon and betaoff */
static const float deacetylate=0.462;
static const float acetylate=0.1155;
static const float PICassembly=0.0277;

static const float startnucleus=0.1;
static const float Kr=10;    /* don't put this less than 1, weird things happen to koff calculation */
static const float GasConstant=8.31447;
static const float cooperativity=1.0;/* dGibbs, relative to 1 additional specific nt */
static const float NumSitesInGenome = 1.8e+6;
static const float selection = 1.0;

static const float mN = 0.1;
static const int Generations=5;

static int output = 0;
static long seed =  28121; /* something is wrong here: changing seed changes nothing */
static int dummyrun=4; /* used to change seed */

/* file output parameters */
static char *output_directory = "output";
int verbose = 0;
FILE *fperrors;
FILE *fp_cellsize;
FILE *fp_growthrate;

void initialize_sequence(char Seq[], 
                         int len)
{
  float x;
  int i;

  if (PLOIDY == 1) {
    for (i=0; i<len; i++){
      x = ran1(&seed);
      if (x<0.25)
        Seq[i] = 'a';
      else if (x<0.5)
        Seq[i] = 'c';
      else if (x<0.75)
        Seq[i] = 'g';
      else Seq[i] = 't';
    }
  } else {
    int current_element = len/(NGENES*PLOIDY);
    int first, second;
    /* printf("current_element=%d\n", current_element); */
    for (i=0; i<len/PLOIDY; i++){
      first = (i / current_element)*2*current_element + i % current_element;
      second = first + current_element;
      /* first = i;
         second = first + len/PLOIDY;
         printf("first=%d, second=%d\n", first, second); */
      x = ran1(&seed);
      if (x<0.25)
        Seq[first] = 'a';
      else if (x<0.5)
        Seq[first] = 'c';
      else if (x<0.75)
        Seq[first] = 'g';
      else Seq[first] = 't';
      /* clone the randomly chosen sequence */
      Seq[second] = Seq[first];
    }
  }

  /* printf("length: %d, sequence is %s\n", strlen(Seq), Seq); */
}


void print_interaction_matrix(AllTFBindingSites *allBindingSites, 
                            int numElements,
                            char transcriptionFactorSeq[NGENES][PLOIDY][TF_ELEMENT_LEN],
                            char cisRegSeq[NGENES][PLOIDY][CISREG_LEN])
{
  int i;

  //printf("cis-reg     gene %2d (copy 0): %.*s\n", i, CISREG_LEN*NGENES*PLOIDY, cisRegSeq[i][0]);

  for (i=0; i < NGENES; i++) {
    printf("TF sequence gene %2d (copy 0): %.*s\n", i, TF_ELEMENT_LEN, transcriptionFactorSeq[i][0]);
    printf("cis-reg     gene %2d (copy 0): %.*s\n", i, CISREG_LEN, cisRegSeq[i][0]);
    if (PLOIDY == 2) {
      printf("TF sequence gene %2d (copy 1): %.*s\n", i, TF_ELEMENT_LEN, transcriptionFactorSeq[i][1]);
      printf("cis-reg     gene %2d (copy 1): %.*s\n", i, CISREG_LEN, cisRegSeq[i][1]);
    }
  } 

  printf("numElements: %3d\n", numElements);
  
  for (i=0; i < numElements; i++) {
    printf("binding site %3d:\n", i);
    printf("       cis-reg region: %3d", allBindingSites[i].cisregID);
    printf("         cis-reg copy: %3d", allBindingSites[i].cisregCopy);
    printf(" (sequence %.*s)\n", CISREG_LEN, cisRegSeq[allBindingSites[i].cisregID][allBindingSites[i].cisregCopy]);
    printf(" transcription-factor: %3d", allBindingSites[i].tfID);
    printf("              TF copy: %3d", allBindingSites[i].tfCopy);
    printf(" (sequence: %.*s)\n", TF_ELEMENT_LEN, transcriptionFactorSeq[allBindingSites[i].tfID][allBindingSites[i].tfCopy]);
    printf("             position: %3d\n", allBindingSites[i].sitePos);
    printf("               strand: %3d\n", allBindingSites[i].strand);
    printf("         Hamming dist: %3d\n", allBindingSites[i].hammingDist); 
    printf("        Hind Position: %3d\n", allBindingSites[i].hindPos); 
    printf("         Start 15 Pos: %3d\n", allBindingSites[i].leftEdgePos);        
  }
}

void initialize_genotype(Genotype *indiv, 
                         float kdis[])
{
  int i, j, p;
  
  initialize_sequence((char *)indiv->cisRegSeq, CISREG_LEN*PLOIDY*NGENES);
  initialize_sequence((char *)indiv->transcriptionFactorSeq, TF_ELEMENT_LEN*PLOIDY*NGENES);
  for(p=0;p<NGENES;p++){
      if(HIND_LENGTH == TF_ELEMENT_LEN){
      indiv->hindrancePositions[p]=0;
      } else {
           indiv->hindrancePositions[p]=rand()%10;
      }
      printf(" %d\n", indiv->hindrancePositions[p]);
  } 
  calc_interaction_matrix(indiv->cisRegSeq, indiv->transcriptionFactorSeq, &(indiv->bindSiteCount), &(indiv->allBindingSites), indiv->hindrancePositions);
  print_interaction_matrix(indiv->allBindingSites, indiv->bindSiteCount, indiv->transcriptionFactorSeq, indiv->cisRegSeq); 
  
  fprintf(fperrors,"activators vs repressors ");
  
  for (i=0; i<NGENES; i++){
    indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    while (indiv->mRNAdecay[i]<0.0)
      indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    indiv->proteindecay[i]=-1.0;
    while (indiv->proteindecay[i]<0.0){
      if (ran1(&seed)<0.08421)
        indiv->proteindecay[i] = 0.0;
      else indiv->proteindecay[i] = exp(0.7874*gasdev(&seed)-3.7665);
    }
    indiv->proteindecay[i] += 0.00578; /* dilution due to cell growth */
    indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    while (indiv->translation[i] < 0.0)
      indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);

    /* make the activations the same in each copy */
    if (ran1(&seed)<pact) {
      indiv->activating[i][0] = 1;
      if (PLOIDY == 2) 
        indiv->activating[i][1] = 1;
    }
    else {
      indiv->activating[i][0] = 0;
      if (PLOIDY == 2) 
        indiv->activating[i][1] = 0;
    }
    fprintf(fperrors,"%d ", indiv->activating[i][0]);
    if (PLOIDY == 2) 
      fprintf(fperrors,"%d ", indiv->activating[i][1]);

    j = trunc(NUM_K_DISASSEMBLY * ran1(&seed));
    
    indiv->PICdisassembly[i][0] = kdis[j];
    if (PLOIDY == 2) 
      indiv->PICdisassembly[i][1] = kdis[j];

  }

  fprintf(fperrors,"\n");
}

void mutate(Genotype *old,
            Genotype *new,
            float m)
{
  int i,k;
  char x;
  
  for (i=0; i<NGENES; i++) {
    for (k=0; k<CISREG_LEN; k++) {
      new->cisRegSeq[i][0][k] = old->cisRegSeq[i][0][k];      
      if (m > ran1(&seed)) {
        x = old->cisRegSeq[i][0][k]; /* because sometimes old and new are the same */
        while (new->cisRegSeq[i][0][k] == x)
          initialize_sequence(&(new->cisRegSeq[i][0][k]),(int) 1);
      }
    }
    for (k=0; k<TF_ELEMENT_LEN; k++) 
      new->transcriptionFactorSeq[i][0][k] = old->transcriptionFactorSeq[i][0][k];  
    new->hindrancePositions[NGENES] = old->hindrancePositions[NGENES];  
    new->mRNAdecay[i] = old->mRNAdecay[i];
    new->proteindecay[i] = old->proteindecay[i];
    new->translation[i] = old->translation[i];
    new->activating[i][0] = old->activating[i][0];
    if (PLOIDY == 2)
      new->activating[i][1] = old->activating[i][1];
    new->PICdisassembly[i][0] = old->PICdisassembly[i][0];
    if (PLOIDY == 2)
      new->PICdisassembly[i][1] = old->PICdisassembly[i][1];
  }
  calc_interaction_matrix(new->cisRegSeq, new->transcriptionFactorSeq, &(new->bindSiteCount), &(new->allBindingSites), new->hindrancePositions);
}

int calc_interaction_matrix_sister(char cisRegSeq[NGENES][PLOIDY][CISREG_LEN],
                                   char transcriptionFactorSeq[NGENES][PLOIDY][TF_ELEMENT_LEN],
                                   int bindSiteCount,
                                   AllTFBindingSites **allBindingSites,
                                   int *maxAlloc,
                                   int cisRegCopy,
                                   int tfCopy,
                                   int hindPos[NGENES])
{
  int i, j, geneID, tfind, match, maxBindingSiteAlloc;

  maxBindingSiteAlloc = *maxAlloc;

  for (geneID=0; geneID < NGENES; geneID++) { /* which cis-reg region */
    for (i=0; i < CISREG_LEN-TF_ELEMENT_LEN; i++) {      /* scan forwards */
      for (tfind=0; tfind < NGENES; tfind++) {
        match=0;
        for (j=i; j < i+TF_ELEMENT_LEN; j++) {
          if (cisRegSeq[geneID][cisRegCopy][j] == transcriptionFactorSeq[tfind][tfCopy][j-i])
            match++;
        }
        if (match >= nmin){
          if (bindSiteCount + 1 >= maxBindingSiteAlloc) {
            maxBindingSiteAlloc = 2*maxBindingSiteAlloc;
            *allBindingSites = realloc(*allBindingSites, maxBindingSiteAlloc*sizeof(AllTFBindingSites));
            if (!(*allBindingSites)) {
              fprintf(fperrors, "realloc of allBindingSites to bindSiteCount = %d failed.\n", maxBindingSiteAlloc);
              exit(1);
            }
            else if (verbose) fprintf(fperrors, "realloc of allBindingSites to bindSiteCount = %d succeeded\n", maxBindingSiteAlloc);
          }
          if(((i - hindPos[tfind]) >=0) && (((i - hindPos[tfind]) + (HIND_LENGTH -1)) < CISREG_LEN)) {
            (*allBindingSites)[bindSiteCount].cisregID = geneID;
            (*allBindingSites)[bindSiteCount].cisregCopy = cisRegCopy; 
            (*allBindingSites)[bindSiteCount].tfID = tfind;
            (*allBindingSites)[bindSiteCount].tfCopy = tfCopy;
            (*allBindingSites)[bindSiteCount].sitePos = i;
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
      for (tfind=0; tfind<NGENES; tfind++) {
        match=0;
        for (j=i; j>i-TF_ELEMENT_LEN; j--)
          if (
             (cisRegSeq[geneID][cisRegCopy][j]=='a' && transcriptionFactorSeq[tfind][tfCopy][i-j]=='t')
             || (cisRegSeq[geneID][cisRegCopy][j]=='t' && transcriptionFactorSeq[tfind][tfCopy][i-j]=='a')
             || (cisRegSeq[geneID][cisRegCopy][j]=='c' && transcriptionFactorSeq[tfind][tfCopy][i-j]=='g')
             || (cisRegSeq[geneID][cisRegCopy][j]=='g' && transcriptionFactorSeq[tfind][tfCopy][i-j]=='c')            
             ) match++;
        if (match >= nmin){
          if (bindSiteCount + 1 >= maxBindingSiteAlloc){
            maxBindingSiteAlloc = 2*maxBindingSiteAlloc;
            *allBindingSites = realloc(*allBindingSites, maxBindingSiteAlloc*sizeof(AllTFBindingSites));
            if (!(*allBindingSites)){
              fprintf(fperrors, "realloc of allBindingSites to bindSiteCount = %d failed.\n", maxBindingSiteAlloc);
              exit(1);
            }
            else if (verbose) fprintf(fperrors, "realloc of allBindingSites to bindSiteCount = %d succeeded\n", maxBindingSiteAlloc);
          }
          if (((i-TF_ELEMENT_LEN+1 - hindPos[tfind]) >=0) && (((i-TF_ELEMENT_LEN+1 - hindPos[tfind])+(HIND_LENGTH-1))< CISREG_LEN)) {
            (*allBindingSites)[bindSiteCount].cisregID = geneID;
            (*allBindingSites)[bindSiteCount].cisregCopy = cisRegCopy; 
            (*allBindingSites)[bindSiteCount].tfID = tfind;
            (*allBindingSites)[bindSiteCount].tfCopy = tfCopy; 
            (*allBindingSites)[bindSiteCount].sitePos = i-TF_ELEMENT_LEN+1;
            (*allBindingSites)[bindSiteCount].strand = 1;
            (*allBindingSites)[bindSiteCount].hammingDist = TF_ELEMENT_LEN-match;
            (*allBindingSites)[bindSiteCount].hindPos = hindPos[tfind];
            (*allBindingSites)[bindSiteCount].leftEdgePos = i-TF_ELEMENT_LEN+1 - hindPos[tfind];
            bindSiteCount++;
          }
        }
      }
    }
  }
  *maxAlloc = maxBindingSiteAlloc;
  return bindSiteCount;
}

void calc_interaction_matrix(char cisRegSeq[NGENES][PLOIDY][CISREG_LEN],
                             char transcriptionFactorSeq[NGENES][PLOIDY][TF_ELEMENT_LEN],
                             int *newBindSiteCount,
                             AllTFBindingSites **allBindingSites,
                             int hindPos[NGENES])
{
  int maxBindingSiteAlloc, bindSiteCount;
  
  maxBindingSiteAlloc = maxelements;
  *allBindingSites = malloc(maxBindingSiteAlloc*sizeof(AllTFBindingSites));
  if (!(*allBindingSites)) {
    fprintf(fperrors,"initial setting of G failed.\n");
    exit(1);
  }
  bindSiteCount = 0;

  bindSiteCount = calc_interaction_matrix_sister(cisRegSeq, 
                                                 transcriptionFactorSeq, 
                                                 bindSiteCount,
                                                 allBindingSites,
                                                 &maxBindingSiteAlloc,
                                                 0, 0, hindPos);

  if (PLOIDY == 2)  {
    
    /* vary the cisRegCopy (1) but keep tfCopy the same (0) (assume no
       gene divergence) */
    bindSiteCount = calc_interaction_matrix_sister(cisRegSeq, 
                                                   transcriptionFactorSeq, 
                                                   bindSiteCount,
                                                   allBindingSites,
                                                   &maxBindingSiteAlloc,
                                                   1, 0, hindPos);

    /* don't currently consider the variant tfCopy combinations (0,1) or (1,1)  */
    /* bindSiteCount = calc_interaction_matrix_sister(cisRegSeq, 
                                                   transcriptionFactorSeq, 
                                                   bindSiteCount,
                                                   allBindingSites,
                                                   &maxBindingSiteAlloc,
                                                   0, 1); */

    /* bindSiteCount = calc_interaction_matrix_sister(cisRegSeq, 
                                                   transcriptionFactorSeq, 
                                                   bindSiteCount,
                                                   allBindingSites,
                                                   &maxBindingSiteAlloc,
                                                   1, 1); */
  }
  /* printf("bindSiteCount = %d\n", bindSiteCount); */

  *allBindingSites = realloc(*allBindingSites, bindSiteCount*sizeof(AllTFBindingSites));
  if (!(*allBindingSites)) {
    fprintf(fperrors, "realloc of allBindingSites down to bindSiteCount = %d failed.\n", bindSiteCount);
    exit(1);
  }
  //printf("bindSiteCount=%d\n", bindSiteCount);
  *newBindSiteCount = bindSiteCount;
}


void add_fixed_event(int i,
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
  newtime->time = t;
  sls_store(newtime, start, last);
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
  newtime->time = t;
  sls_store_end(newtime, start, last);
}

void delete_fixed_event(int geneID,
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
    if (info->geneID==geneID) {
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
  if (done == 0)
    fprintf(fperrors,
            "error: In %d elements, couldn't find element %d to delete in gene %d\n",
            j+1, i, geneID);
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

void initialize_cell(CellState *indiv,
                     /*  int y[NGENES]: AKL 2008-03-21: removed wasn't being used */
                     float mRNAdecay[NGENES],
                     float meanmRNA[NGENES],
                     float initProteinConc[NGENES])
{
  int i, k, totalmRNA;
  float t;

  /* start cell size normalized at 1.0 */
  indiv->cellSize = 1.0;
  
  indiv->mRNATranscrTimeEnd = indiv->mRNATranscrTimeEndLast = NULL;
  indiv->mRNATranslTimeEnd = indiv->mRNATranslTimeEndLast = NULL;
  indiv->tfBoundCount = 0;  /* initialize with nothing bound */
  indiv->tfHinderedCount = 0;
  indiv->tfBoundIndexes = NULL;
  indiv->tfHinderedIndexes = NULL;
  for (i=0; i<NGENES; i++) {
    indiv->active[i][0] = ON_WITH_NUCLEOSOME;
    if (PLOIDY == 2)
      indiv->active[i][1] = ON_WITH_NUCLEOSOME;
    totalmRNA = (int) poidev(meanmRNA[i],&seed);
    indiv->mRNANuclearCount[i] = (int) bnldev(startnucleus, totalmRNA, &seed);
    indiv->mRNACytoCount[i] = totalmRNA - indiv->mRNANuclearCount[i];
    indiv->mRNATranslCytoCount[i] = 0;
    for (k=0; k<indiv->mRNACytoCount[i]; k++) {
      t = expdev(&seed) / mRNAdecay[i];
      if (t < ttranslation) {
        (indiv->mRNACytoCount[i])--;
        (indiv->mRNATranslCytoCount[i])++;
        add_fixed_event(i, ttranslation-t, &(indiv->mRNATranslTimeEnd), &(indiv->mRNATranslTimeEndLast));
      }
    } 
    indiv->mRNATranscrCount[i] = (int) poidev(meanmRNA[i]*ttranscription*mRNAdecay[i],&seed);
    for (k=0; k < indiv->mRNATranscrCount[i]; k++)
      add_fixed_event(i, ran1(&seed)*ttranscription, &(indiv->mRNATranscrTimeEnd), &(indiv->mRNATranscrTimeEndLast));
    indiv->proteinConc[i] = initProteinConc[i];
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
  int i, j;
  
  r = numer = 0.0;

  /* loop over all genes */
  for (i=0; i < NGENES; i++) {
    /* if currently transcribing add it to the rates */
    /* TODO-DONE: check my interpretation of nkonsum */
    if (konStates->nkonsum[i]>0) {
      ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] * t;
      if (fabs(ct)<10^-6) ect=ct;
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

  if (verbose) fprintf(fperrors,"t=%g f=%g df=%g\n", t, *f, *df);
}

void calc_kon_rate(float t,
                   KonStates *konStates,
                   float *konrate)
{
  float r,ct,ect;
  int i;
  
  r = 0.0;
  for (i=0; i<NGENES; i++) {
    if (konStates->nkonsum[i] > 0) {
      ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]*t;
      if (fabs(ct)<10^-6) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) konStates->nkonsum[i])*konStates->konvalues[i][KON_DIFF_INDEX]*ect;
    }
  }
  *konrate = kon*r/t;
  if (verbose) fprintf(fperrors,"r=%g t=%g konrate=%g\n",r,t,*konrate);
}

/* must have already updated proteinConc first */
void change_mRNA_cytoplasm(int i,
                           Genotype *genes,
                           CellState *state,
                           GillespieRates *rates,
                           KonStates *konStates)
{
  float salphc; 
  
  /* TODO: check this */

  /* number of mRNAs in cytoplasm affects ?? */
  salphc = (float) (state->mRNACytoCount[i]) * genes->translation[i] / genes->proteindecay[i];
  rates->salphc += konStates->nkonsum[i]*kon*(salphc - konStates->konvalues[i][KON_SALPHC_INDEX]);
  rates->maxSalphc += konStates->nkonsum[i]*kon*(fmaxf(state->proteinConc[i], salphc) - fmaxf(state->proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]));
  rates->minSalphc += konStates->nkonsum[i]*kon*(fminf(state->proteinConc[i], salphc) - fminf(state->proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]));    
  konStates->konvalues[i][KON_DIFF_INDEX] = (state->proteinConc[i] - salphc) / genes->proteindecay[i];
  konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
}

void calc_koff(int k,
               AllTFBindingSites *allBindingSites,
               CellState *state,
               float *koff)
{
  float Gibbs;  /*free energy in kJ/mol*/
  int posdiff, front, back, i, j;
  
  front = back = 0;
  Gibbs = (((float) allBindingSites[k].hammingDist)/3.0 - 1.0) * state->RTlnKr; /* subject to revision of TF_ELEMENT_LEN */
  for (j=0; j < state->tfBoundCount; j++) {
    if (allBindingSites[k].cisregID==allBindingSites[state->tfBoundIndexes[j]].cisregID &&
        allBindingSites[k].cisregCopy==allBindingSites[state->tfBoundIndexes[j]].cisregCopy &&
        !(k==state->tfBoundIndexes[j])) {
      posdiff = allBindingSites[k].leftEdgePos - allBindingSites[state->tfBoundIndexes[j]].leftEdgePos;
      //printf("diff=%d\n", posdiff);
      if (abs(posdiff) < HIND_LENGTH) {/*Phey*/
        fprintf(fperrors,
                "error: steric hindrance has been breached with site %d (on copy %d of gene %d), %d away from site %d (on copy %d of gene %d)\n",
                k, allBindingSites[k].cisregCopy, allBindingSites[k].cisregID, posdiff, 
                state->tfBoundIndexes[j], allBindingSites[state->tfBoundIndexes[j]].cisregCopy, 
                allBindingSites[state->tfBoundIndexes[j]].cisregID);
      }
      if (abs(posdiff) < 20) {
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
                       float *koffvalues)
{
  int posdiff,j,i;
  float diff;
  
  /* for all the currently bound sites */
  for (j = 0; j < state->tfBoundCount; j++) {
    
    /* we are on the same cisreg sequence and we aren't the same TF binding  */
    if (allBindingSites[indexChanged].cisregID == allBindingSites[state->tfBoundIndexes[j]].cisregID && 
        allBindingSites[indexChanged].cisregCopy == allBindingSites[state->tfBoundIndexes[j]].cisregCopy && 
        !(indexChanged==state->tfBoundIndexes[j])) {

      /* how close are we on the sequence */
      posdiff = allBindingSites[indexChanged].leftEdgePos - allBindingSites[state->tfBoundIndexes[j]].leftEdgePos;
      //printf("diff3=%d\n", posdiff);
      if (abs(posdiff) < HIND_LENGTH) { /* within 6: bad: shouldn't happen Phey*/
        fprintf(fperrors,
                "error: steric hindrance 2 has been breached with site %d %d away from site %d\n",
                indexChanged, posdiff, state->tfBoundIndexes[j]);
      }
      if (abs(posdiff) < 20) {  /* within 20, adjust koff */

        /* save old value */
        diff = -koffvalues[j];

        /* recompute koffvalues */
        calc_koff(state->tfBoundIndexes[j], allBindingSites, state, &(koffvalues[j]));

        /* calculating how koff changes  */
        diff += koffvalues[j];

        /* adjust rates by difference */
        rates->koff += diff;
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
  int i, k;
  
  k = 0;
  /* TODO: check ! */
  /* find the index of the site that binds in konIDs  */
  /* while (!(konStates->konIDs[k][SITEID_INDEX] == siteID) && k < konStates->nkon) {
    k++;
    } */
  
  while (!(konStates->konList[TFID]->available_sites[k] == siteID) && k < konStates->konList[TFID]->site_count) {
    k++;
  }

  if (verbose)
    fprintf(fperrors, ">>> remove site %d konList (k=%d of %d total sites for TF %d, grandtotal=%d)\n", 
            siteID, k, konStates->konList[TFID]->site_count, TFID, konStates->nkon);

  /* make sure that we have enough unoccupied sites left */
  /* if (k < konStates->nkon) {    */
  if (k < konStates->konList[TFID]->site_count && k < konStates->nkon) { 
    /* adjust rates */
    rates->salphc -= kon*salphc;
    rates->maxSalphc -= kon*fmaxf(proteinConcTFID, salphc);
    rates->minSalphc -= kon*fminf(proteinConcTFID, salphc);

    /* one less site available for binding of total */
    (konStates->nkon)--;

    /* also per gene */
    (konStates->nkonsum[TFID])--;

    (konStates->konList[TFID]->site_count)--;

    /* TODO: check */
    /* move the last element end of array into space vacated by site k */
    /* konStates->konIDs[k][SITEID_INDEX] = konStates->konIDs[konStates->nkon][SITEID_INDEX];
       konStates->konIDs[k][TFID_INDEX] = konStates->konIDs[konStates->nkon][TFID_INDEX]; */

    konStates->konList[TFID]->available_sites[k] = konStates->konList[TFID]->available_sites[konStates->konList[TFID]->site_count];
    /*konStates->konIDs[k][TFID_INDEX] = konStates->konIDs[konStates->nkon][TFID_INDEX]; */

  } else {
    if (verbose)
      fprintf(fperrors, "||| couldn't remove site %d from TF %d in konList (k=%d)\n", siteID, TFID, k);
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
  rates->maxSalphc += fmaxf(proteinConcTFID, salphc);
  rates->minSalphc += fminf(proteinConcTFID, salphc);

  /* add back siteID to pool of available sites */
  /* konStates->konIDs[konStates->nkon][SITEID_INDEX] = siteID;
     konStates->konIDs[konStates->nkon][TFID_INDEX] = TFID; */
  
  konStates->konList[TFID]->available_sites[konStates->konList[TFID]->site_count] = siteID;
  /* konStates->konIDs[konStates->nkon][TFID_INDEX] = TFID; */

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
                        int activating[NGENES][PLOIDY],
                        int *on)
{
  int i,off;
  
  *on=off=0;
  for (i=0; i < tfBoundCount; i++) {
    if (geneID==allBindingSites[tfBoundIndexes[i]].cisregID &&
        geneCopy==allBindingSites[tfBoundIndexes[i]].cisregCopy)
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
                     int activating[NGENES][PLOIDY])
{
  int i;
  
  for (i=0; i < tfBoundCount; i++)
    if (geneID==allBindingSites[tfBoundIndexes[i]].cisregID && 
        geneCopy==allBindingSites[tfBoundIndexes[i]].cisregCopy &&
        (activating[allBindingSites[tfBoundIndexes[i]].tfID][geneCopy])) 
      return (1);
  return (0);
}

/* only appropriate if nothing is  bound ie CalcFromInitialState and everything is in state 2
   v31 and earlier have some parts of code needed with prebound stuff,
   rates->mRNAdecay and [7] are done in CalcDt*/
void calc_from_state(Genotype *genes,
                     CellState *state,
                     GillespieRates *rates,
                     KonStates *konStates,
                     float transport[],
                     float mRNAdecay[])
/* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit */
{
  int i, k;
  float salphc, proteinConcTFID; 

  for (i=0; i<NGENES; i++) {
    salphc = (float) (state->mRNACytoCount[i]) * genes->translation[i] / genes->proteindecay[i];
    konStates->konvalues[i][KON_DIFF_INDEX] = (state->proteinConc[i] - salphc) / genes->proteindecay[i];
    konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = genes->proteindecay[i];
    konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
    konStates->nkonsum[i]=0;  
    konStates->konList[i]->site_count = 0;
  }
  state->tfBoundCount=0;

  /* for (i=0; i<7; i++) rates[i]=0.0;    */
  rates->koff=0.0;
  rates->transport=0.0;
  rates->mRNAdecay=0.0;
  rates->picDisassembly=0.0;
  rates->salphc=0.0;
  rates->maxSalphc=0.0;
  rates->minSalphc=0.0;

  for (k=0; k < genes->bindSiteCount; k++) {
    i = genes->allBindingSites[k].tfID;
    proteinConcTFID = state->proteinConc[i];
    salphc = konStates->konvalues[i][KON_SALPHC_INDEX];
    rates->salphc += salphc;
    rates->maxSalphc += fmaxf(proteinConcTFID, salphc);
    rates->minSalphc += fminf(proteinConcTFID, salphc);
    /* konStates->konIDs[k][SITEID_INDEX] = k;
       konStates->konIDs[k][TFID_INDEX] = i; */
    
    /* update the list of sites that bind for a particular TF, i */
    konStates->konList[i]->available_sites[konStates->konList[i]->site_count] = k;
    (konStates->konList[i]->site_count)++;
    (konStates->nkonsum[i])++;
  }

  /* initialize konStates->nkon as the total number of binding sites */
  konStates->nkon = genes->bindSiteCount;

  for (i=0; i<NGENES; i++) {
    transport[i] = kRNA * (float) (state->mRNANuclearCount[i]);
    rates->transport += transport[i];
    
    /* start all genes in acteylated state */
    state->statechangeIDs[ACTEYLATION][0][i] = i;
    if (PLOIDY == 2)
      state->statechangeIDs[ACTEYLATION][1][i] = i;
  }
  /* for (i=4; i<7; i++) rates[i] *= kon; */
  rates->salphc *= kon;
  rates->maxSalphc *= kon;
  rates->minSalphc *= kon;

  rates->acetylationCount[0]=NGENES;
  rates->deacetylationCount[0]=0;
  rates->picAssemblyCount[0]=0;
  rates->transcriptInitCount[0]=0;
  rates->picDisassemblyCount[0]=0;

  if (PLOIDY == 2) {
    rates->acetylationCount[1]=NGENES;
    rates->deacetylationCount[1]=0;
    rates->picAssemblyCount[1]=0;
    rates->transcriptInitCount[1]=0;
    rates->picDisassemblyCount[1]=0;
  }

}

/* returns:
 *  0 if there is no fixed event occuring before time t
 *  1 if a transcription event happens before time t
 *  2 if a translation event happens before time t
 */
int does_fixed_event_end(FixedEvent *mRNATranslTimeEnd,
                         FixedEvent *mRNATranscrTimeEnd,
                         float t)
{
  if (mRNATranslTimeEnd==NULL) {
    if (mRNATranscrTimeEnd==NULL) return(0);
    else{
      if (mRNATranscrTimeEnd->time<t) return(1);
      else return(0);
    }
  } else {
    if (mRNATranscrTimeEnd==NULL){
      if (mRNATranslTimeEnd->time<t) return(2);
      else return(0);
    } else {
      if (mRNATranscrTimeEnd->time < mRNATranslTimeEnd->time){
        if (mRNATranscrTimeEnd->time < t) return(1);
        else return(0);
      } else {
        if (mRNATranslTimeEnd->time < t) return(2);
        else return(0);
      }
    }
  }
}

void calc_dt(float *x,
             float *dt,
             GillespieRates *rates,
             KonStates *konStates,
             float mRNAdecay[],
             float mRNAdecayrates[],
             int mRNACytoCount[],
             int mRNATranslCytoCount[])
{
  float tbound1,tbound2;
  int i;

  /* reset the total rate for current step */
  rates->total=0.0;
  
  /* reset mRNA decay rate */
  rates->mRNAdecay=0.0;

  /* update mRNAdecay rate based on the total number of mRNAs in both
     cytoplasm (mRNACytoCount) and ones that have only just recently arrived
     (mRNATranslCytoCount) */
  for (i=0; i<NGENES; i++){
    mRNAdecay[i] = mRNAdecayrates[i] * ((float) mRNACytoCount[i] + (float) mRNATranslCytoCount[i]);
    rates->mRNAdecay += mRNAdecay[i];
  }
  /* AKL: should this go from 1-6 ? rather than 1-5?  Answer: no
   * because 5, 6 are cached values, not part of rates used by
   * Gillespie algorithm.
   */

  /* recompute and cache the total rate in data structure */
  rates->total += rates->koff;
  rates->total += rates->transport;
  rates->total += rates->mRNAdecay;
  rates->total += rates->picDisassembly;
  rates->total += rates->salphc;
  /* OLD, was: for (i=0; i<5; i++) rates->total += rates[i]; */

  /* 
   * convert the counts back into rates using the constants 
   * TODO-DONE: check
   */
  rates->total += (float) rates->acetylationCount[0] * acetylate;
  rates->total += (float) rates->deacetylationCount[0] * deacetylate;
  rates->total += (float) rates->picAssemblyCount[0] * PICassembly;
  rates->total += (float) rates->transcriptInitCount[0] * transcriptinit;    

  if (PLOIDY == 2) {
    rates->total += (float) rates->acetylationCount[1] * acetylate;
    rates->total += (float) rates->deacetylationCount[1] * deacetylate;
    rates->total += (float) rates->picAssemblyCount[1] * PICassembly;
    rates->total += (float) rates->transcriptInitCount[1] * transcriptinit;    
  }

  /* TODO: check logic here of tbound1,2 */  
  tbound1 = *x/(rates->total + rates->maxSalphc);
  tbound2 = *x/(rates->total + rates->minSalphc);
  if (verbose) 
    fprintf(fperrors,"bounds %g %g\n",tbound1,tbound2);

  /* if bounds are the same, simply choose tbound1 */
  if (tbound1==tbound2){
    if (konStates->nkon!=0)
      fprintf(fperrors,
              "error: nkon=%d when it should be zero x=%f rates->maxSalphc=%g rates->minSalphc=%g rates->total=%g\n",
              konStates->nkon, *x, rates->maxSalphc, rates->minSalphc, rates->total);
    *dt = tbound1;
  } else {
    /* otherwise get delta t by solving the equation using Newton-Raphson method */
    *dt = rtsafe(&calc_time, *x, rates, konStates, tbound1, tbound2, (float) 1e-6); 
  }
}

void end_transcription(float *dt,
                       float t,
                       CellState *state,
                       float transport[NGENES],
                       GillespieRates *rates)
{
  int i, total;
  
  /* recompute the delta-t based on difference between now and the
     time of transcription end */
  *dt = state->mRNATranscrTimeEnd->time - t;
  total = 0;
  for (i=0; i < NGENES; i++) 
    total += state->mRNATranscrCount[i];
  
  if (verbose) 
    fprintf(fperrors,"\ntranscription event finishes out of %d possible t=%g dt=%g\n",
            total, t, *dt);

  /* get the gene which is ending transcription */
  i = state->mRNATranscrTimeEnd->geneID;

  /* increase number of mRNAs in nucleus */
  (state->mRNANuclearCount[i])++;

  /* decrease the number of mRNAs undergoing transcription */
  (state->mRNATranscrCount[i])--;

  /* delete the fixed even which has just occurred */
  delete_fixed_event_start(&(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));

  /* add rate kRNA to transport and Gillespie rates
   * TODO: check */
  transport[i] += kRNA;
  rates->transport += kRNA;
}

void transport_event(float x,
                     float transport[NGENES],
                     CellState *state,
                     float endtime,
                     GillespieRates *rates)
{
  int i;
  float konrate2;
  
  i = -1;
  konrate2 = 0.0;  /* TODO-DONE: check, why is this set to zero?, should be set to zero fixed */

  /* choose gene product (mRNA) that gets transported to cytoplasm
     based on weighting in transport[] array */
  while (i < NGENES && x > konrate2){
    i++;
    x -= transport[i];
    /* TODO: isn't konrate2 or something supposed to computed here? */
  }
  if (verbose) {
    fprintf(fperrors, "transport event gene %d from %d copies\n", i, state->mRNANuclearCount[i]);
  }

  /* one less mRNAs in nucleus */
  (state->mRNANuclearCount[i])--;

  /* it has just arrived in cytoplasm, ready to be translated */
  (state->mRNATranslCytoCount[i])++;
  
  /* add the endtime for translation */
  add_fixed_event_end(i, endtime, &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));

  /* decrease transport frequency */
  transport[i] -= kRNA;

  /* if below a threshold, make zero:
   * TODO: check */
  if (transport[i] < 0.1*kRNA) 
    transport[i]=0.0;

  /* adjust rates */
  rates->transport -= kRNA;

  /* do similar threshold check */
  if (rates->transport < 0.1*kRNA) 
    rates->transport=0.0;
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
    if (force) 
      /* don't always print because with a 4->3 transition PIC assembly is not there to be removed */
      fprintf(fperrors, "error removing %d from array of length %d, type=%d\n", toberemoved, *len, type);
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
  
  /* disassemble PIC in OFF state */
  if (state->active[geneID][geneCopy] == OFF_PIC) {
    (state->active[geneID][geneCopy]) = OFF_NO_PIC;
    state->statechangeIDs[DEACTEYLATION][geneCopy][rates->deacetylationCount[geneCopy]] = geneID;
    (rates->deacetylationCount[geneCopy])++;
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
    state->statechangeIDs[ACTEYLATION][geneCopy][rates->acetylationCount[geneCopy]] = geneID;
    (rates->acetylationCount[geneCopy])++;      
  }
  /* OFF_NO_PIC -> ON_NO_PIC */
  if ((transcriptrule) && oldstate==OFF_NO_PIC) {
    state->active[geneID][geneCopy] = ON_NO_PIC;
    remove_from_array(geneID, DEACTEYLATION,  state->statechangeIDs[DEACTEYLATION][geneCopy], &(rates->deacetylationCount[geneCopy]), (int) 1);
    if (numactive){
      state->statechangeIDs[PICASSEMBLY][geneCopy][rates->picAssemblyCount[geneCopy]] = geneID;
      (rates->picAssemblyCount[geneCopy])++;
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
    remove_from_array(geneID, ACTEYLATION, state->statechangeIDs[ACTEYLATION][geneCopy], &(rates->acetylationCount[geneCopy]), (int) 1);
  }
  
  /* ON_NO_PIC -> OFF_NO_PIC */
  if (!(transcriptrule) && oldstate==ON_NO_PIC){          
    state->active[geneID][geneCopy] = OFF_NO_PIC;
    remove_from_array(geneID, PICASSEMBLY, state->statechangeIDs[PICASSEMBLY][geneCopy], &(rates->picAssemblyCount[geneCopy]), (int) 0);
    state->statechangeIDs[DEACTEYLATION][geneCopy][rates->deacetylationCount[geneCopy]] = geneID;
    (rates->deacetylationCount[geneCopy])++;
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
    fprintf(fperrors, "state change from %d to %d in gene %d, copy %d\n", oldstate, state->active[geneID][geneCopy], geneID, geneCopy);
  }
}

void remove_tf_binding(Genotype *genes,
                       CellState *state,
                       GillespieRates *rates,
                       KonStates *konStates,
                       int site,
                       float koffvalues[])
{
  int i, j, k, bound, siteID, geneID, geneCopy, transcriptrule, oldstate, numactive;

  i = 0;

  /* given site 'site', look for the index in the list of bound sites */
  while ((state->tfBoundIndexes[i] != site) && (i < state->tfBoundCount)) 
    i++;
  if (i == state->tfBoundCount) {  /* couldn't find the site */
    fprintf(fperrors, "error: remove_tf_binding could not find site %d with %d possibilities\n Bound sites are\n",
            site, state->tfBoundCount);
    for (j = 0; j < state->tfBoundCount; j++)  {
      fprintf(fperrors, "%d\n", state->tfBoundIndexes[j]);
    }
  }
  else {
    j = 0;
    /* loop through the sterically hindered sites */
    while (j < state->tfHinderedCount) {

      /* TODO: check logic here */
      /* check all sites hindered by binding to location 'site' */
      if (state->tfHinderedIndexes[j][1] == site) {
        k = bound = 0;

        /* is anything else hindering the same site? */
        while (bound == 0 && k < state->tfHinderedCount) {
          if (state->tfHinderedIndexes[j][0] == state->tfHinderedIndexes[k][0] && j != k) 
            bound=1;
          k++;
        }

        /* if nothing else is hindering this site then allow site to be bound */
        if (bound==0) {
          siteID = state->tfHinderedIndexes[j][0];
          if (verbose) 
            fprintf(fperrors,"Site %d leftEdgePos %d on gene %d freed from steric hindrance\n",
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
      } else { /* only increment if we haven't shortened array (TODO: cleanup)  */
        j++;
      }
    }    

    /* reduce the koff rate by the mount */
    rates->koff -= koffvalues[i];

    /* one less bound site */
    (state->tfBoundCount)--;

    /* shift end of array to hole opened up */
    state->tfBoundIndexes[i] = state->tfBoundIndexes[state->tfBoundCount];

    /* likewise with koffvalues */
    koffvalues[i] = koffvalues[state->tfBoundCount];

    /* find the gene and copy whose cisreg region has an unbinding event */
    geneID = genes->allBindingSites[site].cisregID;
    geneCopy = genes->allBindingSites[site].cisregCopy;
    if (verbose) 
      fprintf(fperrors,"Add site %d at leftEdgePos %d on gene %d copy %d freed by unbinding\n",
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
    scan_nearby_sites(site, genes->allBindingSites, state, rates, koffvalues);        
  }
}

void attempt_tf_binding(Genotype *genes,
                        CellState *state,
                        GillespieRates *rates,
                        float **koffvalues,
                        KonStates *konStates,
                        int *maxbound2,
                        int *maxbound3,
                        int site)
{
  int geneID, geneCopy, k, posdiff, site2;

  if (verbose) fprintf(fperrors,
                       "kon1 event at site %d out of %d possible, %d TFs previously bound bindSiteCount=%d\n",
                       site, konStates->nkon, state->tfBoundCount, genes->bindSiteCount);

  fflush(fperrors);

  /* if we have run out of space, double memory  */
  if (state->tfBoundCount >= *maxbound2){
    (*maxbound2) *= 2;

    state->tfBoundIndexes = realloc(state->tfBoundIndexes, (*maxbound2)*sizeof(int));

    /* do the copy */
    *koffvalues = realloc(*koffvalues, (*maxbound2)*sizeof(float));

    /* check return value */
    if (!state->tfBoundIndexes || !(*koffvalues)) {
      fprintf(fperrors, "memory allocation error resetting maxbound2=%d\n", *maxbound2);
      exit(1);
    }
  }

  /* append the site to end of indexes */
  state->tfBoundIndexes[state->tfBoundCount] = site;
  if (verbose) fprintf(fperrors, "remove site %d\n", site);

  fflush(fperrors);

  /* remove the site from the kon pool */
  remove_kon(site,
             genes->allBindingSites[site].tfID,
             rates, 
             konStates->konvalues[genes->allBindingSites[site].tfID][KON_SALPHC_INDEX],
             konStates,
             state->proteinConc[genes->allBindingSites[site].tfID]);

  /* recompute the koffvalues */
  calc_koff(site, genes->allBindingSites, state, &((*koffvalues)[state->tfBoundCount]));
  if (verbose) 
    fprintf(fperrors,"new koff = %g is number %d\n",
            (*koffvalues)[state->tfBoundCount], (state->tfBoundCount+1));

  /* adjust rates by adding the new koffvalue to rates->koff */
  rates->koff += (*koffvalues)[state->tfBoundCount];

  /* TODO: check, this seems the same as above ???  maybe remove?? */
  state->tfBoundIndexes[state->tfBoundCount] = site;

  /* increment number of bound TFs */
  (state->tfBoundCount)++;

  /* adjust co-operative binding in context of new TF */
  scan_nearby_sites(site, genes->allBindingSites, state, rates, *koffvalues);

  /* get the gene that the TF is binding to */
  geneID = genes->allBindingSites[site].cisregID;

  /* get the copy that the TF is binding to */
  geneCopy = genes->allBindingSites[site].cisregCopy;
  
  /* update steric hindrance data structures */
  /* JM: this cycles over all sites, not just bound ones, in order to
     record redundancy in steric hindrance*/
  for (k = 0; k < genes->bindSiteCount; k++) {

    /* if we are on the same gene and not the same binding site */
    if (geneID == genes->allBindingSites[k].cisregID &&
        geneCopy == genes->allBindingSites[k].cisregCopy &&
        !(k==site)) {

      /* check distance from current binding site (k) to the original (site) */
                      //fprintf(fperrors, "site=%d k=%d\n", genes->allBindingSites[site].leftEdgePos, genes->allBindingSites[k].leftEdgePos);
      posdiff = genes->allBindingSites[site].leftEdgePos - genes->allBindingSites[k].leftEdgePos;

      /* if within 6, we prevent binding by adding to steric hindrance */
      if (abs(posdiff) < HIND_LENGTH) {/*Phey*/
                    //fprintf(fperrors, "diff < 15\n");
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
        if (verbose) {
          fprintf(fperrors, "%d steric hindrance sites after %d blocks site %d\n", state->tfHinderedCount, site, k);
          fflush(fperrors);
        }

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
  if (verbose) fprintf(fperrors, "tfBoundCount=%d tfHinderedCount=%d maxbound2=%d maxbound3=%d\n",
                       state->tfBoundCount, state->tfHinderedCount, *maxbound2, *maxbound3);

  /* gene activity may change as a result of binding */
  revise_activity_state(geneID, geneCopy, genes, state, rates);
}

/*
 * time course of [TF]s represented as array of TimeCourse lists.
 */
void add_time_points(float time,
                     float proteinConc[NGENES],
                     TimeCourse **timecoursestart,
                     TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i<NGENES; i++)
    add_time_point(time, proteinConc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

void add_integer_time_points(float time,
                             int proteinConc[NGENES],
                             TimeCourse **timecoursestart,
                             TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i<NGENES; i++)
    add_time_point(time, (float) proteinConc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

/* 
 * need some sort of control in case it declines to essentially zero.
 * Add in discrete, stochastic and/or zero values, but this may create
 * false attractor if time steps are small and rising tide keeps
 * getting rounded down
 */
void update_protein_conc(float proteinConc[],
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
  
  rates->maxSalphc = rates->minSalphc = 0.0;
  for (i=0; i<NGENES; i++) {
    ct = konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]*dt;
    ect = exp(-ct);
    if (fabs(ct)<10^-6) ect1=ct;
    else ect1 = 1-ect;   
    proteinConc[i] = konStates->konvalues[i][KON_SALPHC_INDEX]*ect1 + ect*proteinConc[i];
    konStates->konvalues[i][KON_DIFF_INDEX] = (proteinConc[i] - konStates->konvalues[i][KON_SALPHC_INDEX]) / konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX];
    rates->maxSalphc += ((float) konStates->nkonsum[i]) * fmaxf(proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]);
    rates->minSalphc += ((float) konStates->nkonsum[i]) * fminf(proteinConc[i], konStates->konvalues[i][KON_SALPHC_INDEX]);
  }
  rates->maxSalphc *= kon;
  rates->minSalphc *= kon;
  if ((output) && (*timecourselast)->time < t+dt-0.1)
    add_time_points(t+dt,otherdata/*proteinConc*/,timecoursestart,timecourselast);
}

void calc_num_bound(float proteinConc[],
                    int tfBoundCount)
{
  float sum;
  int i;
  
  sum = 0.0;
  for (i=0; i < NGENES; i++) 
    sum += proteinConc[i];
  if (verbose) 
    /* fprintf(fperrors, "%d bound %g expected\n", tfBoundCount, 0.0003*sum);*/
    /* TODO: check! */
    /* if this is wrong for random sequences adjust Kr accordingly */
    fprintf(fperrors, "%d bound %g expected\n", tfBoundCount, (CISREG_LEN*NGENES*sum)/NumSitesInGenome);
}

int sum_rate_counts(int rate_array[PLOIDY])
{
  if (PLOIDY == 2) {
    return (rate_array[0] + rate_array[1]);
  } else {
    return (rate_array[0]);
  }
}

void get_gene(int rate_array[PLOIDY], int pos, int *geneLoc, int *geneCopy)
{
  
  if (PLOIDY == 2) {
    //int sum = rate_array[0] + rate_array[1];
    if (pos < rate_array[0]) {
      *geneCopy = 0;
      *geneLoc = pos;
    } else {
      *geneCopy = 1;
      *geneLoc = pos - rate_array[0];
    }
  } else {
    *geneCopy = 0;
    *geneLoc = pos;
  }
  /* printf("rate_array[0]=%d, rate_array[1]=%d, pos=%d, geneCopy=%d, geneLoc=%d\n", 
     rate_array[0], rate_array[1], pos, *geneCopy, *geneLoc); */
}


/* -----------------------------------------------------
 * START
 * Functions that handle each possible Gillespie event 
 * ----------------------------------------------------- */
void tf_binding_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                      KonStates *konStates, float *koffvalues, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float konrate, float dt, float t, int maxbound2, int maxbound3)
{
  float x = ran1(&seed) * (rates->salphc + konrate)/kon;
  int i, j = -1;
  int k, l = -1;  /* new */
  float total_konrate2, konrate2_for_TF = 0.0;     
  int siteID;

  /* loop through all *available* binding sites to
   * choose one, the interval is weighted by the
   * frequency of rates */
  /* while (j < konStates->nkon-1 && x > konrate2) { */

  /* loop through all TFs, then choose a particular binding site */
  for (k=0; k < NGENES; k++) {

    /* compute the total rate for all available binding sites for this
     * particular TF: see if we are in the right range
     */

    /* if no sites available for this TF, skip to next gene */
    if (konStates->konList[k]->site_count == 0) {
      if (verbose)
        fprintf(fperrors, "looking at TF: %d, has %d sites available, skipping\n", k, konStates->konList[k]->site_count);
      continue;
    }
    
    /* first, cache the konrate2 for this particular gene */
    konrate2_for_TF = konStates->konvalues[k][KON_SALPHC_INDEX] + 
      konStates->konvalues[k][KON_DIFF_INDEX] * (1-exp(-konStates->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt))/dt;

    /* compute the *total* kon rate for all unbound sites for this TF  */
    total_konrate2 = ((konStates->konList[k]->site_count)) * konrate2_for_TF;  /* remove +1 for the moment */

    if (verbose)
      fprintf(fperrors, "looking at TF: %d, has %d sites available [konrate2: %g, total_konrate2: %g, x: %g]\n", 
              k, konStates->konList[k]->site_count, konrate2_for_TF, total_konrate2, x); 

    /* if we are already in the appropriate TF, now choose a binding site */
    if (!(x > total_konrate2) || (k == NGENES-1)) {
      float konrate2 = 0.0;
      
      if (verbose) 
        fprintf(fperrors, "selecting TF: %d, konrate2: %g, total_konrate2: %g, x: %g\n", k, konrate2_for_TF, total_konrate2, x); 
      
      while (l < konStates->konList[k]->site_count && x > konrate2) {
        /* this will record the last binding site before we
           reach the appropriate binding site  */
        /* j++; */
        l++;

        /* get ID of site */
        siteID = konStates->konList[k]->available_sites[l];
        /* get ID of TF: now already held implictly */
        /* i = konStates->konIDs[j][TFID_INDEX]; */
        /* i = k; */
        
        if (verbose)
          fprintf(fperrors, "l: %d, site: %d, binds to TF: %d, x = %g\n", l, siteID, k, x); 

        /* compute the kon rate  */
        /* konrate2 = konStates->konvalues[k][KON_SALPHC_INDEX] + 
           konStates->konvalues[k][KON_DIFF_INDEX]*(1-exp(-konStates->konvalues[k][KON_PROTEIN_DECAY_INDEX]*dt))/dt; */
        
        /* adjust random number */
        /* x -= konrate2; */
        konrate2 = konrate2_for_TF;
        x -= konrate2;
      }
      /* found it, so break out of the outer for loop */
      break;
    } else {
      x -= total_konrate2; 
      /* printf("progressing to the next TF: %d\n", k+1); */
    }

  }

  /* printf("found it!!!! l: %d, site: %d, binds to TF: %d, konrate2: %g, x: %g\n", l, siteID, k, konrate2_for_TF, x);   */
  /* update protein concentration before doing the binding */
  update_protein_conc(state->proteinConc, dt, 
                      rates, konStates, 
                      t, timecoursestart, timecourselast,
                      state->proteinConc);
  if (verbose) fflush(fperrors);
  
  /* bind the j-th konID in the list which happens to be the binding site in SITEID_INDEX */
  /* attempt_tf_binding(genes, state, rates, &koffvalues, konStates, &maxbound2, &maxbound3, 
     konStates->konIDs[j][SITEID_INDEX]); */
  
  /* bind siteID */
  attempt_tf_binding(genes, state, rates, &koffvalues, konStates, &maxbound2, &maxbound3, siteID);

  /* calculate the number of TFs bound
   * TODO: check this seems to be purely diagnostic
   */
  calc_num_bound(state->proteinConc, state->tfBoundCount);
}

void tf_unbinding_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                      KonStates *konStates, float *koffvalues, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                        float konrate, float dt, float t, float x)
{
  int i, j = -1;
  int site;

  while (j < state->tfBoundCount && x>0) {
    j++;
    x -= koffvalues[j];
  }
  if (j==state->tfBoundCount) {
    float konrate2 = 0.0;
    for (i = 0; i < state->tfBoundCount; i++) 
      konrate2 += koffvalues[i];
    fprintf(fperrors, "warning: koffvalues add up to %g instead of rates->koff=%g\n",
            konrate2, rates->koff);
    rates->koff = konrate2;
    j--; /* a bit of a fudge for rounding error, really should move on to rates->transport, but too complicated for something so minor */
  } 
  site = state->tfBoundIndexes[j];
  if (verbose) fprintf(fperrors, "koff event %d of %d at site %d\n",
                       j,state->tfBoundCount,site);
  if (j < 0) fprintf(fperrors, "error: koff event %d of %d at site %d\n", j, state->tfBoundCount,site);
  
  /* update protein concentration before removing TF */
  update_protein_conc(state->proteinConc, dt, 
                    rates, konStates, 
                    t, timecoursestart, timecourselast, 
                    state->proteinConc);
  
  /* remove TF binding from 'site' */
  remove_tf_binding(genes, state, rates, konStates, site, koffvalues);
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
    fprintf(fperrors, "warning: x=%g > konrate2=%g out of rates->mRNAdecay=%g\n",
            x, konrate2, rates->mRNAdecay);
  }

  /* assume mRNA cytoplasm transport events equally likely */
  x = ran1(&seed)*((float) (state->mRNACytoCount[i] + state->mRNATranslCytoCount[i]));
  update_protein_conc(state->proteinConc, dt,
                    rates, konStates, 
                    t, timecoursestart, timecourselast,
                    state->proteinConc);
  
  /* 
   * decay mRNA in cytoplasm 
   */
  if (x < (float)state->mRNACytoCount[i]) {
    if (verbose) {
      fprintf(fperrors,"mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
              i, state->mRNACytoCount[i], state->mRNATranslCytoCount[i]);
    }
    /* remove the mRNA from the cytoplasm count */
    (state->mRNACytoCount[i])--;  
    change_mRNA_cytoplasm(i, genes, state, rates, konStates); 
    
  } else {
    /* 
     * decay mRNA in process of translating (TODO: check!) 
     */
    x = ran1(&seed)*((float) state->mRNATranslCytoCount[i]);
    if (verbose){
      fprintf(fperrors,
              "mRNA decay event gene %d not from %d copies in cytoplasm but %f from %d copies translating\n",
              i, state->mRNACytoCount[i], trunc(x), state->mRNATranslCytoCount[i]);
    }
    
    /* delete this fixed event: this mRNA will never be translated */
    delete_fixed_event(i, (int) trunc(x), &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));
    
    /* remove the mRNA from the count */
    (state->mRNATranslCytoCount[i])--;
    if (verbose) for (j=0; j < NGENES; j++)
                   fprintf(fperrors, "%d copies of gene %d translating\n", 
                           state->mRNATranslCytoCount[j], j);                  
  }
}

void histone_acteylation_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                               KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                               float dt, float t)
{
  int geneLoc, geneCopy;
  float x = ran1(&seed)*((float) sum_rate_counts(rates->acetylationCount));
  
  /* geneLoc = (int)trunc(x) / PLOIDY;
     geneCopy = (int)trunc(x) % PLOIDY; */
  get_gene(rates->acetylationCount, (int)trunc(x), &geneLoc, &geneCopy);

  /* choose a particular gene to change state */
  int geneID = state->statechangeIDs[ACTEYLATION][geneCopy][geneLoc];

  /* choose a particular gene to change state */
  /* int geneCopy = choose_gene_copy(); */
  
  if (verbose) fprintf(fperrors,"acetylation event gene %d\nstate change from %d to 4\n",
                       geneID, state->active[geneID][geneCopy]);
  if (state->active[geneID][geneCopy] != ON_WITH_NUCLEOSOME)
    fprintf(fperrors, "error: acetylation event attempted from state %d\n", state->active[geneID][geneCopy]);
  update_protein_conc(state->proteinConc, dt,
                      rates, konStates, 
                      t, timecoursestart, timecourselast,
                      state->proteinConc);
  
  /* set state: eject nucleosome, but there is no PIC yet */
  state->active[geneID][geneCopy] = ON_NO_PIC;
  remove_from_array(geneID, ACTEYLATION, state->statechangeIDs[ACTEYLATION][geneCopy], &(rates->acetylationCount[geneCopy]), (int) 1);
  if (is_one_activator(geneID, geneCopy, state->tfBoundIndexes, state->tfBoundCount, 
                       genes->allBindingSites, genes->activating)) {
    state->statechangeIDs[PICASSEMBLY][geneCopy][rates->picAssemblyCount[geneCopy]] = geneID; 
    (rates->picAssemblyCount[geneCopy])++;
  }
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                                 KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                                 float dt, float t)
{
  float x = ran1(&seed)*((float) sum_rate_counts(rates->deacetylationCount));

  int geneCopy; /* = (int)trunc(x) % PLOIDY; */
  int geneLoc; /* = (int)trunc(x) / PLOIDY; */

  get_gene(rates->deacetylationCount, (int)trunc(x), &geneLoc, &geneCopy);

  /* choose a particular gene and copy to change state */
  int geneID = state->statechangeIDs[DEACTEYLATION][geneCopy][geneLoc];
  /* int geneCopy = choose_gene_copy(); */

  if (verbose) fprintf(fperrors,"deacetylation event gene %d, copy %d\nstate change from %d to 1\n",
                       geneID, geneCopy, state->active[geneID]);
  if (state->active[geneID][geneCopy] != OFF_NO_PIC)
    fprintf(fperrors, "error: deacetylation event attempted from state %d\n", state->active[geneID][geneCopy]);
  update_protein_conc(state->proteinConc, dt, 
                      rates, konStates, 
                      t, timecoursestart, timecourselast,
                      state->proteinConc);
  /* set state: nucleosome returns */
  state->active[geneID][geneCopy] = OFF_FULL;
  remove_from_array(geneID, DEACTEYLATION, state->statechangeIDs[DEACTEYLATION][geneCopy], &(rates->deacetylationCount[geneCopy]), (int) 1);
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                      KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                      float dt, float t)
{

  float x = ran1(&seed)*((float) sum_rate_counts(rates->picAssemblyCount));

  int geneCopy; /* = (int)trunc(x) % PLOIDY; */
  int geneLoc; /*= (int)trunc(x) / PLOIDY; */

  get_gene(rates->picAssemblyCount, (int)trunc(x), &geneLoc, &geneCopy);

  /* choose a particular gene and copy to change state */
  int geneID = state->statechangeIDs[PICASSEMBLY][geneCopy][geneLoc];
  /* int geneCopy = choose_gene_copy(); */

  if (verbose) fprintf(fperrors, "PIC assembly event gene %d copy %d\nstate change from %d to 6\n",
                       geneID, geneCopy, state->active[geneID][geneCopy]);

  if (state->active[geneID][geneCopy] != ON_NO_PIC)
    fprintf(fperrors, "error: PIC assembly event attempted from state %d\n", state->active[geneID][geneCopy]);
  update_protein_conc(state->proteinConc, dt,
                      rates, konStates, 
                      t, timecoursestart, timecourselast,
                      state->proteinConc);
  
  /* turn gene fully on: ready for transcription */
  state->active[geneID][geneCopy] = ON_FULL;
  remove_from_array(geneID, PICASSEMBLY, state->statechangeIDs[PICASSEMBLY][geneCopy], &(rates->picAssemblyCount[geneCopy]), (int) 1);
  state->statechangeIDs[TRANSCRIPTINIT][geneCopy][rates->transcriptInitCount[geneCopy]] = geneID;
  (rates->transcriptInitCount[geneCopy])++;
  state->statechangeIDs[PICDISASSEMBLY][geneCopy][rates->picDisassemblyCount[geneCopy]] = geneID;
  (rates->picDisassemblyCount[geneCopy])++;
  rates->picDisassembly += genes->PICdisassembly[geneID][geneCopy];
}

void disassemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genes, 
                           KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                           float dt, float t, float x)
{
  int geneCopy, geneLoc, geneID;
  int j=-1;
  while (j < NGENES*PLOIDY && x>0) {
    j++;

    get_gene(rates->picDisassemblyCount, j, &geneLoc, &geneCopy);

    /* TODO: why do we keep original rand number x, rather than choose a new one like other functions?  */
    x -= genes->PICdisassembly[state->statechangeIDs[PICDISASSEMBLY][geneCopy][geneLoc]][geneCopy];
  }
  if (j==NGENES*PLOIDY) fprintf(fperrors, "error in PIC disassembly\n");
  geneID = state->statechangeIDs[PICDISASSEMBLY][geneCopy][geneLoc];
  if (verbose) fprintf(fperrors, "PIC disassembly event in copy %d of gene %d\n", geneCopy, geneID);
  disassemble_PIC(state, genes, geneID, geneCopy, rates);
}

void transcription_init_event(GillespieRates *rates, CellState *state, Genotype *genes,
                              KonStates *konStates, TimeCourse **timecoursestart, TimeCourse **timecourselast,
                              float dt, float t, float x)
{
  int geneID;

  x /= transcriptinit;

  int geneCopy; /* = (int)trunc(x) % PLOIDY; */
  int geneLoc; /* = (int)trunc(x) / PLOIDY; */

  get_gene(rates->transcriptInitCount, (int)trunc(x), &geneLoc, &geneCopy);

  /* choose the gene and copy that gets transcribed */
  geneID = state->statechangeIDs[TRANSCRIPTINIT][geneCopy][geneLoc];
  /* geneCopy = choose_gene_copy(); */

  if (verbose) fprintf(fperrors,"transcription event gene %d, copy %d\n", geneID, geneCopy);

  if (state->active[geneID][geneCopy] != ON_FULL && state->active[geneID][geneCopy] != OFF_PIC)
    fprintf(fperrors, "error: transcription event attempted from state %d\n", state->active[geneID][geneCopy]);

  update_protein_conc(state->proteinConc, dt,
                      rates, konStates, 
                      t, timecoursestart, timecourselast,
                      state->proteinConc);

  /* now that transcription of gene has been initiated, 
   * we add the time it will end transcription, 
   * which is dt+time of transcription from now */
  add_fixed_event_end(geneID, t+dt+ttranscription, 
                      &(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));

  /* increase the number mRNAs being transcribed */
  (state->mRNATranscrCount[geneID])++;                      
}
/* -----------------------------------------------------
 * END
 * Functions that handle each possible Gillespie event 
 * ----------------------------------------------------- */

float compute_growth_rate(CellState *cell_state, Genotype *genes, float t, float dt) {
  float growth_rate_mRNA, growth_rate_prot;
  float benefit_term, cost_term_prot, cost_term_mRNA;
  int geneID = 9;      /* select a particular TF */
  float Lp = 12064.28; /* mean gene expression */
  float Lm = 1589836;  /* max gene expression */
  float gpeak = 2*9.627*1e-5;

  //float cost = 0.00001;   /* Wagner (2005) */

  /* alternative 1, assuming gpeak = 1 */

  //float cost = 1/(Lm-2*Lp);
  //float gmax = ((Lm-Lp)*(Lm-Lp))/((Lm-2*Lp)*(Lm-2*Lp));
  //float Kmax = Lp*Lp/(Lm-2*Lp);

  /* alternative 2, assuming c = -log(1-s)/L*dt */
  //float cost = -log(1-1e-5);
  //float gmax = cost*(Lm+(Lp*Lp/(Lm-2*Lp)));
  //float Kmax = Lp*Lp/(Lm-2*Lp);

  /* alternative 3 */
  float hc = 2.583*1e-9;
  float cost = hc;
  float gmax = gpeak + 2*cost*Lp + (pow(cost,2)*pow(Lp,2))/gpeak;
  float Kmax = (cost*pow(Lp,2))/gpeak;
  //float gmax = 0.000138787;
  //float Kmax = 2408.23;

  benefit_term = (gmax*cell_state->proteinConc[geneID])/(cell_state->proteinConc[geneID] + Kmax);
  cost_term_prot = - cost*cell_state->proteinConc[geneID];
  cost_term_mRNA = - cost*(genes->translation[geneID]/genes->proteindecay[geneID])*cell_state->mRNACytoCount[geneID];

  /* alternative 4: converted mRNA count into protein amount in above term for cost */
  growth_rate_mRNA = benefit_term + cost_term_mRNA;
  growth_rate_prot = benefit_term + cost_term_prot;

  /* make sure growth rate can't be negative */
  if (growth_rate_mRNA < 0.0)
    growth_rate_mRNA = 0.0;

  if (growth_rate_prot < 0.0)
    growth_rate_prot = 0.0;

  /* printf("protein=%g gmax=%g, Kmax=%g, cost=%g\ncost_prot=%g, cost_mRNA=%g (mRNA=%d), growth rate=%g\n", 
         cell_state->proteinConc[0], gmax, Kmax, cost, cost_term_prot, cost_term_mRNA, 
         cell_state->mRNACytoCount[geneID], growth_rate); */
  
  fprintf(fp_growthrate, "%g %g %g\n", t, growth_rate_mRNA, growth_rate_prot);

  return (growth_rate_mRNA);
}

float compute_growth_rate_dimer(CellState *cell_state, Genotype *genes, float t, float dt) {

  float growth_rate_mRNA, growth_rate_prot;
  float benefit_term, cost_term_prot, cost_term_mRNA;
  int geneID = 9;      /* select a particular TF */
  float L = cell_state->proteinConc[geneID];  /* cache value */

  float Lp = 12000; /* mean gene expression is 12064.28 */
  float gpeak = 2*9.627*1e-5;
  float mean_decay_rate = 0.023;  /* in min^-1 from mean of distribution from Belle et al (2006) */
  float hc = 2.583*1e-9;
  float cost = hc;
  float M0 = (pow(gpeak,2) - hc*pow(gpeak,2) + 2*hc*gpeak*Lp - 2*pow(hc,2)*gpeak*Lp +  pow(hc*Lp,2) - pow(hc,3)*pow(Lp,2))/gpeak;
  float Km = (hc*pow(gpeak,2) - 2*hc*gpeak*Lp + 2*pow(hc,2)*gpeak*Lp + hc*pow(Lp,2) - 2*pow(hc*Lp,2) + pow(hc,3)*pow(Lp,2))/gpeak;

  benefit_term = 0.5 * (Km + L + M0 - sqrt(pow(Km + L + M0, 2) - 4*L*M0));
  cost_term_prot = - cost*L;
  cost_term_mRNA = - cost*(genes->translation[geneID]/genes->proteindecay[geneID])*cell_state->mRNACytoCount[geneID];

  growth_rate_mRNA = benefit_term + cost_term_mRNA;
  growth_rate_prot = benefit_term + cost_term_prot;

  /* printf("protein=%g M0=%g, Km=%g, cost=%g\nbenefit=%g, cost_prot=%g, cost_mRNA=%g (mRNA=%d), growth rate mRNA=%g, growth rate protein=%g\n", 
         cell_state->proteinConc[0], M0, Km, cost, benefit_term, cost_term_prot, cost_term_mRNA, 
         cell_state->mRNACytoCount[geneID], growth_rate_mRNA, growth_rate_prot);  */

  /* make sure growth rate can't be negative */
  if (growth_rate_mRNA < 0.0)
    growth_rate_mRNA = 0.0;

  if (growth_rate_prot < 0.0)
    growth_rate_prot = 0.0;
  
  fprintf(fp_growthrate, "%g %g %g\n", t, growth_rate_mRNA, growth_rate_prot);

  return (growth_rate_mRNA);
}

void update_cell_size(CellState *cell_state, Genotype *genes, float t, float dt) {
  
  float growth_rate = compute_growth_rate_dimer(cell_state, genes, t, dt);
  cell_state->cellSize = (cell_state->cellSize)*exp(growth_rate*dt);
  fprintf(fp_cellsize, "%g %g\n", t, cell_state->cellSize);

}

/*
 * develop: run the cell for a given length of time
 */

void develop(Genotype *genes,
             CellState *state,
             float temperature, /* in Kelvin */
             TimeCourse **timecoursestart,
             TimeCourse **timecourselast)
{
  float t;
  int i, j, k;

  /* UNUSED here remove: int posdiff; */

  int maxbound2, maxbound3;  /* TODO: check these */
  int site;      /* site where TF gets removed (TODO: check)  */
  int geneID;    /* ID of gene which is undergoing transitions in transcription state */
  int event;     /* boolean to keep track of whether FixedEvent has ended */

  float *koffvalues;   /* rates of unbinding */
  int total;           /* total possible translation events, TODO: rename */
  float transport[NGENES];  /* transport rates of each mRNA */
  float mRNAdecay[NGENES];  /* mRNA decay rates */
  float x;                  /* random number */
  float dt;                 /* delta-t */

  /* cached information about available binding sites for efficiency */
 
  KonStates *konStates =  malloc(sizeof(KonStates));

  GillespieRates *rates = malloc(sizeof(GillespieRates));

  /* stores corresponding geneIDs for [de]acteylation, PIC[dis]assembly, transcriptinit */
  /* int statechangeIDs[5][NGENES]; */

  float f, df, konrate, konrate2, diff, sum, ct, ect;

  /* initialize time courses */
  for (i=0; i<NGENES; i++){
    timecoursestart[i] = NULL;
    timecourselast[i] = NULL;
  }
  
  add_time_points((float) 0.0, state->proteinConc, timecoursestart, timecourselast);

  /* set cell temperature and value of RTlnKr constant */
  state->temperature = temperature;
  state->RTlnKr = GasConstant * temperature * log(Kr);
   
  /* number of possible binding sites */
  /* konStates->konIDs = malloc(2*genes->bindSiteCount*sizeof(int)); */
  /* konStates->konList = malloc(sizeof(KonList)); */

  for (i=0; i<NGENES; i++){
    konStates->konList[i] = malloc(sizeof(KonList));
    konStates->konList[i]->available_sites = malloc(genes->bindSiteCount*sizeof(int));
  }
  
  /* TODO: ? */
  maxbound2 = maxbound;
  maxbound3 = 10*maxbound;

  state->tfBoundIndexes = realloc(state->tfBoundIndexes, maxbound2*sizeof(int));
  koffvalues = malloc(maxbound2*sizeof(float)); 
  state->tfHinderedIndexes = realloc(state->tfHinderedIndexes, 2*maxbound3*sizeof(int));
  
  /*if (!konStates->konvalues || !state->tfBoundIndexes || !koffvalues ||
    !state->tfHinderedIndexes || !konStates->konIDs) { */
  if (!konStates->konvalues || !state->tfBoundIndexes || !koffvalues ||
      !state->tfHinderedIndexes || !konStates->konList) {
    fprintf(fperrors,"memory allocation error at start of develop\n");
    exit(1);
  }

  /* initialize transcriptional state of genes */
  calc_from_state(genes, state, rates, konStates, transport, mRNAdecay);

  t = 0.0;  /* time starts at zero */
  

  while (t < tdevelopment) {  /* run until development stops */
  
    x=expdev(&seed);        /* draw random number */

    /* TODO: what is this doing? */
    if (rates->koff < 0.0){
      konrate2 = 0.0;
      for (i=0; i < state->tfBoundCount; i++) konrate2 += koffvalues[i];
      if ((verbose) || konrate2>0.0)
        fprintf(fperrors,"warning: koffvalues add up to %g rates->koff=%g < 0\n",
                konrate2, rates->koff);
      rates->koff = konrate2;
    }
    
    /* do first Gillespie step to chose next event */

    calc_dt(&x, &dt, rates, konStates, mRNAdecay, genes->mRNAdecay,
           state->mRNACytoCount, state->mRNATranslCytoCount);

    if (verbose) 
      fprintf(fperrors,"next stochastic event due at t=%g dt=%g x=%g\n", t+dt, dt, x);

    if (!(state->mRNATranscrTimeEndLast)) {
      for (i=0;i<NGENES;i++)
        if (verbose) fprintf(fperrors, "%d transcription events\n", state->mRNATranscrCount[i]);
    }
    
    /* first check to see if a fixed event occurs in current t->dt window */
    event = does_fixed_event_end(state->mRNATranslTimeEnd,
                                 state->mRNATranscrTimeEnd,
                                 fminf(t+dt, tdevelopment));

    
    /* while there are either transcription or translation events
       occuring in current t->dt window */
    while (event > 0) {
      konrate = x/dt;

      if (event==1) {  /* if a transcription event ends */
        end_transcription(&dt, t, state, transport, rates);

        update_protein_conc(state->proteinConc, dt,
                            rates, konStates, t,
                            timecoursestart, timecourselast,
                            state->proteinConc);
      } else {           /* if a translation event ends */
        /* TODO: could have else if (event==2) */
        dt = state->mRNATranslTimeEnd->time - t;         /* make dt window smaller */
        total=0;  /* number of translation events */

        /* count current number of mRNAs that have recently arrived in cytoplasm */
        for (i=0; i<NGENES; i++) total += state->mRNATranslCytoCount[i];
        if (verbose) 
          fprintf(fperrors,"\ntranslation event finishes out of %d possible t=%g dt=%g\n",
                  total, t, dt); /* bug: dt can be negative */

        /* get identity of gene that has just finished translating */
        i=state->mRNATranslTimeEnd->geneID;   

        /* there is one less mRNA that has just finished translation */
        (state->mRNATranslCytoCount[i])--;   

        /* delete the event that just happened */
        delete_fixed_event_start(&(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));

        /* there is one more mRNA that is now in cytoplasm */
        (state->mRNACytoCount[i])++;

        /* update protein concentration */
        update_protein_conc(state->proteinConc, dt,
                          rates, konStates, t,
                          timecoursestart, timecourselast,
                          state->proteinConc);
        
        /* the number of mRNAs in cytoplasm affects binding, TODO: check*/
        change_mRNA_cytoplasm(i, genes, state, rates, konStates);
      }

      /* advance time by the dt */
      t += dt;
      x -= dt*konrate;
      if (verbose) 
        fprintf(fperrors, "dt=%g t=%g fixed event old x=%g new x=%g\n", dt, t, x+dt*konrate, x);

      /* re-compute a new dt */
      calc_dt(&x, &dt, rates, konStates, mRNAdecay, 
             genes->mRNAdecay, state->mRNACytoCount, state->mRNATranslCytoCount);

      if (verbose) 
        fprintf(fperrors,"next stochastic event (2) due at t=%g dt=%g x=%g\n", t+dt, dt, x);

      /* check to see there aren't more fixed events to do */
      event = does_fixed_event_end(state->mRNATranslTimeEnd, 
                                   state->mRNATranscrTimeEnd, 
                                   fminf(tdevelopment, t+dt));
    } 

    /* no remaining fixed events to do in dt, now do stochastic events */

    /* if we haven't already reached end of development with last delta-t */
    if (t+dt < tdevelopment) {

      /* compute total konrate (which is constant over the Gillespie step) */
      if (konStates->nkon==0) konrate = (-rates->salphc);
      else calc_kon_rate(dt, konStates, &konrate); 

      /* 
       * choose a new uniform?? (TODO) random number weighted by the
       * probability of all Gillespie events, note that konrate is
       * *not* included in rates->total, so it needs to be added here
       */
      x = ran1(&seed)*(rates->total + konrate);  

      if (verbose){
        fprintf(fperrors,"\nx=%g\tfBoundCount=%g = %d * %g\ntransport=%g\ndecay=%g\n",
                x, rates->koff, state->tfBoundCount, rates->koff/(float)state->tfBoundCount, 
                rates->transport, rates->mRNAdecay);
        fprintf(fperrors,"PICdisassembly=%g\nkon=%g = %d * %g\n",
                rates->picDisassembly, rates->salphc+konrate, konStates->nkon, (rates->salphc+konrate)/(float)konStates->nkon);
        fprintf(fperrors,"acetylation=%g (copy 0)\ndeacetylation=%g (copy 0)\nPIC assembly=%g (copy 0)\ntranscriptinit=%g (copy 0)\n",
                (float)rates->acetylationCount[0]*acetylate, (float)rates->deacetylationCount[0]*deacetylate, 
                (float)rates->picAssemblyCount[0]*PICassembly, (float)rates->transcriptInitCount[0]*transcriptinit);
        if (PLOIDY == 2) {
          fprintf(fperrors,"acetylation=%g (copy 1)\ndeacetylation=%g (copy 1)\nPIC assembly=%g (copy 1)\ntranscriptinit=%g (copy 1)\n",
                  (float)rates->acetylationCount[1]*acetylate, (float)rates->deacetylationCount[1]*deacetylate, 
                  (float)rates->picAssemblyCount[1]*PICassembly, (float)rates->transcriptInitCount[1]*transcriptinit);
        }
        fprintf(fperrors,"total=%g=%g+%g\n\n", rates->total + konrate, rates->total, konrate);
      }
      /* JM: kon generally could be handled better, with more direct
       * references to konStates->nkonsum, probably a bit vulnerable to rounding
       * error
       */

      /* 
       * STOCHASTIC EVENT: a TF unbinds (koff) 
       */
      if (x < rates->koff) {  
        tf_unbinding_event(rates, state, genes, konStates, koffvalues,
                           timecoursestart, timecourselast, konrate, dt, t, x);
      } else {
        x -= rates->koff;  
        /* 
         * STOCHASTIC EVENT: a transport event
         */
        if (x < rates->transport) {     
          update_protein_conc(state->proteinConc, dt, 
                              rates, konStates, 
                              t, timecoursestart, timecourselast, 
                              state->proteinConc);
          transport_event(x, transport, state, t+dt+ttranslation, rates);
        } else {
            
          x -= rates->transport;
          /* 
           * STOCHASTIC EVENT: an mRNA decay event
           */
          if (x < rates->mRNAdecay) {  
            mRNA_decay_event(rates, state, genes, konStates, mRNAdecay,
                             timecoursestart, timecourselast, dt, t, x);
          } else {
            x -= rates->mRNAdecay;
            /* 
             * STOCHASTIC EVENT: PIC disassembly
             */
            if (x < rates->picDisassembly) {
              disassemble_PIC_event(rates, state, genes, konStates, 
                                    timecoursestart,  timecourselast, dt, t, x);
            } else {
              x -= rates->picDisassembly;
              /* 
               * STOCHASTIC EVENT: TF binding event
               */
              if (x < rates->salphc + konrate) {   /* add variable (salphc) and constant (konrate) */
                tf_binding_event(rates, state, genes, konStates, koffvalues,
                                 timecoursestart, timecourselast, konrate, dt, t, 
                                 maxbound2, maxbound3);
              } else {
                x -= (rates->salphc + konrate);
                /* 
                 * STOCHASTIC EVENT: histone acetylation
                 */
                if (x < (float) sum_rate_counts(rates->acetylationCount) * acetylate) {
                     
                  histone_acteylation_event(rates, state, genes, konStates, 
                                            timecoursestart, timecourselast, dt, t);
                } else {
                   
                  x -= (float) sum_rate_counts(rates->acetylationCount) * acetylate;
                  /* 
                   * STOCHASTIC EVENT: histone deacetylation
                   */
                  if (x < (float) sum_rate_counts(rates->deacetylationCount) * deacetylate) {
              
                    histone_deacteylation_event(rates, state, genes, konStates, 
                                                timecoursestart, timecourselast, dt, t);
                  } else {
                    x -= (float) sum_rate_counts(rates->deacetylationCount) * deacetylate;
                    /* 
                     * STOCHASTIC EVENT: PIC assembly
                     */
                    if (x < (float) sum_rate_counts(rates->picAssemblyCount) * PICassembly) {
                      assemble_PIC_event(rates, state, genes, konStates, 
                                         timecoursestart, timecourselast, dt, t);
                    } else {
                      x -= (float) sum_rate_counts(rates->picAssemblyCount) * PICassembly;
                      /* 
                       * STOCHASTIC EVENT: transcription initiation
                       */
                      if (x < (float) sum_rate_counts(rates->transcriptInitCount) * transcriptinit) {
                        transcription_init_event(rates, state, genes, konStates, 
                                                 timecoursestart, timecourselast, dt, t, x);
                      } else {
                        /*
                         * FALLBACK: shouldn't get here, previous
                         * events should be exhaustive
                         */
                          
                        fprintf(fperrors, "error: no event assigned\n");
                      }
                    }
                  }
                }
              }
            }       
          }
        }
      }
      
      /* update cell size */
      update_cell_size(state, genes, t, dt);
      
      /* Gillespie step: advance time to next event at dt */
      t += dt;
      if (verbose) fprintf(fperrors, "dt=%g t=%g\n", dt, t);
 
      /* output time courses during run to capture output before segfault */
      for (i=0; i < NGENES; i++) {
        float *lopt;
        /*print_time_course(timecoursestart[i], i, lopt);*/
      }
 
    } else {
      /* we will reach the end of development in dt */
      if (verbose) fprintf(fperrors, "finish at t=%g dt=%g\n", t, dt);

      /* do remaining dt */
      dt = tdevelopment - t;

      /* final update of protein concentration */
      update_protein_conc(state->proteinConc, dt,
                        rates, konStates,
                        t, timecoursestart, timecourselast,
                        state->proteinConc);
      /* advance to end of development (this exits the outer while loop) */
      t = tdevelopment;
    }
  }
  free(koffvalues);
  for (i=0; i<NGENES; i++) {
    free(konStates->konList[i]->available_sites);
    free(konStates->konList[i]);
  }
  free(konStates);
  free(rates);
}

void dev_stability_only_lopt(float lopt[],
                          TimeCourse **timecoursestart)
{
  int i;
  TimeCourse *start;
  float dt1,dt2;
  
  for (i=0;i<NGENES;i++){
    lopt[i]=0.0;
    start = timecoursestart[i];
    dt1=0.0;
    while (start->next){
      dt2 = (start->next->time - start->time)/2.0;
      lopt[i] += start->concentration * (dt1+dt2);     
      dt1=dt2;
      start = start->next;
    }
    lopt[i] += start->concentration * dt1;
    lopt[i] /= tdevelopment;
    fprintf(fperrors,"lopt[%d] = %g\n",i,lopt[i]);
  }
}

void calc_fitness(float lopt[],
                 float *w,
                 TimeCourse **timecoursestart, 
                 float s)
{
  float d,dt1,dt2,x;
  int i;
  TimeCourse *start[NGENES];
  
  for (i=0;i<NGENES;i++) start[i]=timecoursestart[i];
  dt1=0.0;
  *w=0.0;
  while (start[0]){
    d = 0.0;
    if (start[0]->next) dt2 = (start[0]->next->time - start[0]->time)/2.0;
    else dt2=0.0;
    for (i=0;i<NGENES;i++){
      x = (start[i]->concentration - lopt[i]) / lopt[i];
      d += x*x;
    }
    d = sqrt(d);
    *w += exp(-s*d) * (dt1+dt2);     
    if (verbose) fprintf(fperrors,"t=%g dt=%g+%g d=%g w=%g\n",start[0]->time,dt1,dt2,d,*w/tdevelopment);
    dt1=dt2;
    for (i=0;i<NGENES;i++) start[i] = start[i]->next;
  }
  *w /= tdevelopment;
  fprintf(fperrors,"s=%g w=%g\n",s,*w);
}

void print_time_course(TimeCourse *start,
                       int i,
                       float lopt[])
{
  FILE *fpout;
  char filename[80];
  
  sprintf(filename, "%s/protein%d.dat", output_directory, i);
  if ((fpout = fopen(filename,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n",filename);
  while (start){
    fprintf(fpout,"%g %g\n", start->time, start->concentration);
    start = start->next;
  }
  /* fprintf(fpout, "%g %g\n", tdevelopment, lopt[i]);  */
  fclose(fpout);  
}

int main(int argc, char *argv[])
{
  FILE *fpout, *fpkdis;
  char fperrors_name[80];
  char fp_cellsize_name[80];
  char fp_growthrate_name[80];
  int i, j, k, gen;
  CellState state;
  Genotype indivs[PopSize];
  TimeCourse *timecoursestart[NGENES]; /* array of pointers to list starts */
  TimeCourse *timecourselast[NGENES];
  TimeCourse *start;
  float fitness[PopSize], sumfit, lopt[NGENES], initmRNA[NGENES], initProteinConc[NGENES], x, kdis[NUM_K_DISASSEMBLY];

  int c, directory_success;
  int hold_genotype_constant = 0;
  int curr_seed;
 
  /* parse command-line options */
  while ((c = getopt (argc, argv, "hvgd:r:")) != -1) {
    switch (c)
      {
      case 'd':
        output_directory = optarg;
        break;
      case 'r':
        dummyrun = atoi(optarg);
        break;
      case 'g':
        hold_genotype_constant = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        fprintf(stderr, "%s [-d DIRECTORY] [-r DUMMYRUN] [-h] [-g]\n", argv[0]);
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

  /* create output files for  growth rate */
  sprintf(fp_growthrate_name, "%s/growthrate.dat", output_directory);
  if ((fp_growthrate = fopen(fp_growthrate_name,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", fp_growthrate_name);

  sumfit = 0.0;

  /* slight hack to initialize seed  */

  if (!hold_genotype_constant)
    for (curr_seed=0; curr_seed<dummyrun; curr_seed++) ran1(&seed);

  /* initialize protein concentrations */
  for (i=0; i<NGENES; i++) {
    lopt[i] = exp(1.25759*gasdev(&seed)+7.25669);
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
   
    initialize_cell(&state, indivs[j].mRNAdecay, initmRNA, initProteinConc);
    /* 
     *  AKL 2008-03-21: removed indivs[j].y: wasn't being used
     *  initialize_cell(&state,indivs[j].y,indivs[j].mRNAdecay,initmRNA,initProteinConc); 
     */
    develop(&indivs[j], &state, (float) 293.0, timecoursestart, timecourselast);
    /* dev_stability_only_lopt(lopt, timecoursestart); */
    fprintf(fperrors,"indiv %d\n",j);
    /* calc_fitness(lopt, &(fitness[j]), timecoursestart, selection);
       sumfit += fitness[j];*/
    for (i=0; i < NGENES; i++) {
      if ((output) && j==PopSize-1) print_time_course(timecoursestart[i], i, lopt);
      if (verbose) fprintf(fperrors, "deleting gene %d\n", i);
      delete_time_course(timecoursestart[i]);
      timecoursestart[i] = timecourselast[i] = NULL;
    }
    free_mem_CellState(&state);
  }

  for (j = 0; j < PopSize; j++) {
    free(indivs[j].allBindingSites);
  }

  /*
    for (gen=1;gen<=Generations*PopSize;gen++){
    fprintf(fperrors,"meanfit = %g\n",sumfit / (float) PopSize);
    x = ran1(&seed) * sumfit;
    j = -1;
    while (x>0){
      x -= fitness[j+1];
      j++;
    }
    k = (int) trunc(ran1(&seed)*PopSize);
    fprintf(fperrors,"gen=%d indiv %d reproduces fitness %g replaces indiv %d\n",
      gen,j,fitness[j] * (float)PopSize / sumfit,k);
    fflush(fperrors);
    mutate(&(indivs[j]),&(indivs[k]),mN/(float)PopSize);
    fprintf(fperrors,"finished mutating\n");fflush(fperrors);
    initialize_cell(&state,indivs[k].y,indivs[k].mRNAdecay,initmRNA,initProteinConc);
    develop(&indivs[k],&state,(float) 293.0,timecoursestart,timecourselast);
    fprintf(fperrors,"finished development\n");fflush(fperrors);
    //    dev_stability_only_lopt(lopt,timecoursestart);
    sumfit -= fitness[k];
    calc_fitness(lopt,&(fitness[k]),timecoursestart,selection);
    fflush(fperrors);
    sumfit += fitness[k];
    for (i=0; i<NGENES; i++){
      if (output && gen==Generations*PopSize) print_time_course(timecoursestart[i],i,lopt);
      if (verbose) fprintf(fperrors,"deleting gene %d\n",i);
      delete_time_course(timecoursestart[i]);
      timecoursestart[i] = timecourselast[i] = NULL;
    }
    free_mem_CellState(&state);   
  }
  */
  fclose(fperrors);
  fclose(fp_cellsize);
  fclose(fp_growthrate);
}
