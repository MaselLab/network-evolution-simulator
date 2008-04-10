/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- 
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
#include "random.h"
//#include "list.h"

/*
  version27: each founding indiv has a different genotype, CopyGenotype is obsolete
  version28: add selection (wrong sort)
  v29: transcription initiation quenching rule, take selection out
  v30reject: add delays to quenching and unquenching
  v30: redo selection based on v29, buggy
  v30alt: redo selection based on v27, free up memory of CellState and TimeCourse
  v31: based on v30alt, replace stability only selection with selection for specific opt.
  Hopefully this helps avoid evolution from always-on to always-off (not sufficient yet).
  Move meanmRNA and protein initial states from InitializeCell to be universal
  v32: koffIDs and nkoff have the same info as B2 and y2: remove them
  v33: implement transcription initiation scheme with bursting, PIC assembly & H2A.Z
  v34: disassembly rates bootstrapped from constitutively on Garcia-Martinez subset, not TF subset distribution
  v35: output only a subset of timecourse measurements to prevent excessive output
  change NumSitesInGenome value
  v36: reglen changed to 150
*/
#define MAXIT 100

int maxelements=500; 
//start by allocating maxelements when initializing a genotype, double as needed, reduce at end
int maxbound=50;
#define NGenes 10
#define reglen 150
#define elementlen 6
#define Nkdis 133
int dummyrun=4; //used to change seed
int PopSize=1;
int nmin=4;
float tdevelopment=120.0;
float kon=1e-4; // lower value is so things run faster
/*kon=0.2225 is based on 1 molecule taking 240seconds=4 minutes
  and 89% of the proteins being in the nucleus*/
float kRNA=618.0;
float ttranslation=1.0;
float ttranscription=1.0;
float pact=0.62;
float transcriptinit=8.5; //replace betaon and betaoff
float deacetylate=0.462;
float acetylate=0.1155;
float PICassembly=0.0277;

float startnucleus=0.1;
float Kr=10;//don't put this less than 1, weird things happen to koff calculation
float GasConstant=8.31447;
float cooperativity=1.0;//dGibbs, relative to 1 additional specific nt
float NumSitesInGenome = 1.8e+6;
float selection = 1.0;
int verbose = 1;
int output = 0;
float mN = 0.1;
int Generations=5;

long seed =  28121; //something is wrong here: changing seed changes nothing
FILE *fperrors;

/* random number generator from Numerical Recipes*/
float ran1(long *seed);
/*Normally distributed random number with mean zero and variance 1*/
float gasdev(long *seed);
float poidev(float xm, long *seed);
/*Returns an exponentially distributed, positive, random deviate of unit mean*/
float expdev(long *seed);
float bnldev(float pp, int n, long *seed);

/* enum for 'konvalues' indices
 * are rates of binding with:
 * element KON_DIFF          is (proteinConc - salphc)/proteindecay
 * element KON_PROTEIN_DECAY is proteindecay
 * element KON_SALPHC        is salphc
 */
enum { KON_DIFF_INDEX = 0, KON_PROTEIN_DECAY_INDEX = 1, KON_SALPHC_INDEX = 2 };

/*
 * enum for 'CellState'->active indices
 */
enum { OFF_FULL = 1,           /* repr>activ, still nuclesome, and PIC */
       ON_WITH_NUCLEOSOME = 2, /* activ>repr, still nucleosome */
       OFF_NO_PIC = 3,         /* repr>activ, no nucleosome, but no PIC */
       ON_NO_PIC = 4,          /* activ>repr, no nucleosome, but no PIC  */
       OFF_PIC = 5,            /* repr>activ, no nucleosome, and PIC */
       ON_FULL = 6 };          /* activ>repr, no nucleosome, and a PIC: ready to go!  */

/*
 * enum for konIDs
 */
enum { SITEID_INDEX = 0, TFID_INDEX = 1 };

/*
 * New data structure for keeping track of rates for Gillespie
 * algorithm 
 *
 * Events with exponentially-distributed waiting times are:
 * - TF association and dissociation, 
 * - transcription initiation
 * - mRNA transport and decay
 * - acetylation and deacetylation, 
 * - and PIC assembly and disassembly.
 *
 * was formerly float rates[8];  and int rates2[5] 
 * rates: total rates for 0-koff,1-transport,2-mRNAdecay,3-PIC disassembly,
 *  4-salphc,5-max(L,salphc),6-min(L,salphc) c>0 and 7-total incl. rates2
 *  3-PIC disassembly used to be 3-transcriptinit
 * rates2: #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 
 *  3-transcriptinit, 4=PIC disassembly
 */
struct GillespieRates
{
  float koff;              /* rates[0] */
  float transport;         /* rates[1] */
  float mRNAdecay;         /* rates[2] */
  float picDisassembly;    /* rates[3] */
  float salphc;            /* rates[4] */

  /* the following are cached values to avoid recomputation not
     actually part of the Gillespie rates */
  float maxSalphc;         /* rates[5] */
  float minSalphc;         /* rates[6] */

  /* number of genes in the following states */
  int acetylationCount;       /* rates2[0] */
  int deacetylationCount;     /* rates2[1] */
  int picAssemblyCount;       /* rates2[2] */
  int transcriptInitCount;    /* rates2[3] */
  int picDisassemblyCount;    /* rates2[4] */

  /* total, including rates2 */
  float total;                /* rates[7] */
};

/*
 * konStates : new composite data structure to cache information about
 * available binding sites to avoid re-computation.  This currently
 * just groups the previous data structures: 
 *
 * konIDs, konvalues, nkonsum
 *
 * but will ultimately be refactored again for greater efficiency
 */
struct KonStates {
  /* number of currently *available* binding sites */
  int nkon;

  /* elem 0 is siteID, elem 1 is TF that binds */
  int (*konIDs)[2];

  /* number of available binding sites for a given TF */
  int nkonsum[NGenes];

  /* konvalues are rates of binding with:
   * element 0 is other funny term (TODO ?)
   * element 1 is c
   * element 2 is salphc
   * Other index is which TF binds 
   * The kon term is left out of all of them for computational efficiency
   */
  float konvalues[NGenes][3];
};


/* Newton-Raphson root-finding method with bisection steps, out of
 * Numerical Recipes function bracketed by x1 and x2. Returns root
 * within accuracy +/-xacc funcd is function of interest, returning
 * both function value and first deriv.*/
float rtsafe(void (*funcd)(float, float, struct GillespieRates *, struct KonStates *, float *, float *), 
             float x, struct GillespieRates *rates, struct KonStates *konStates, float x1, float x2, float xacc)
{
  int j,done;
  float df,dx,dxold,f,fh,fl,xtemp;
  float temp,xh,xl,rts;
  
  (*funcd)(x1, x, rates, konStates, &fl, &df);
  (*funcd)(x2, x, rates, konStates, &fh, &df); /* note df isn't used here */
  if (fabs(fl) < 1e-9) return x1;
  if (fabs(fh) < 1e-9) return x2;
  if ((fl > 0.0 && fh > 0.0) || (fl <0.0 && fh < 0.0)){
    if (verbose) fprintf(fperrors,"warning in rtsafe: root should be bracketed\n");
    if (fabs(fl) < fabs(fh)) return x1; else return x2;
  }
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts, x, rates, konStates, &f, &df);
  done = 0;
  for (j=1;j<=MAXIT;j++){
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) {
      done = 1;// modified: otherwise this bisection can mean 2 identical function calls for j=1
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;      
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    if (rts==0.0){
      if (x1<x2) rts = fminf(2.0*x1,(xl+xh)/2.0);
      else rts = fminf(2.0*x2,(xl+xh)/2.0);
      fprintf(fperrors,"warning: dt=0 reset to %g\n",rts);
    }
    if (j>1 || done==0) (*funcd)(rts, x, rates, konStates, &f, &df);
    if (f < 0.0) xl=rts;
    else xh=rts;   
  }
  fprintf(fperrors,"error in rtsafe: too many iterations\n");
  return 0.0;
}

void InitializeSeq(char Seq[], 
                   int len)
{
  float x;
  int i;
  
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
  /* printf("length: %d, sequence is %s\n", strlen(Seq), Seq); */
}

struct Genotype
{
  char cisRegSeq[NGenes][reglen];
  char transcriptionFactorSeq[NGenes][elementlen];
  int bindSiteCount;
  struct TFInteractionMatrix *interactionMatrix;
/* int (*interactionMatrix)[5];
 5 elements are
    identity of cis-regulatory region
    identity of TF that binds
    position 0 to 386 (first position, always in forwards strand)  
    strand 0 (forward) or 1 (backward)
    Hamming distance 0 to n-n_min 
*/ 
  float mRNAdecay[NGenes];
  float proteindecay[NGenes];
  float translation[NGenes];
  int activating[NGenes]; // 1 is activating, 0 is repressing
  float PICdisassembly[NGenes];
};

struct TFInteractionMatrix
{
  int cisregID;     /* 0 */
  int tfID;         /* 1 */
  int sitePos;      /* 2 */
  int strand;       /* 3 */
  int hammingDist;  /* 4 */
};


void CalcInteractionMatrix(char [NGenes][reglen],
                           char [NGenes][elementlen],
                           int *,
                           struct TFInteractionMatrix **);

void PrintInteractionMatrix(struct TFInteractionMatrix *interactionMatrix, 
                            int numElements,
                            char transcriptionFactorSeq[NGenes][elementlen],
                            char cisRegSeq[NGenes][reglen])
{
  int i;

  for (i=0; i<NGenes; i++) {
    printf("transcriptionFactorSeq number %2d: %.*s\n", i, elementlen, transcriptionFactorSeq[i]);
    printf("cis-reg                number %2d: %.*s\n", i, reglen, cisRegSeq[i]);
  } 
  
  for (i=0; i<numElements; i++) {
    printf("binding site %3d:\n", i);
    printf("       cis-reg region: %3d", interactionMatrix[i].cisregID);
    printf(" (sequence %.*s\n)", reglen, cisRegSeq[interactionMatrix[i].cisregID]);
    printf(" transcription-factor: %3d", interactionMatrix[i].tfID);
    printf(" (sequence: %.*s)\n", elementlen, transcriptionFactorSeq[interactionMatrix[i].tfID]);
    printf("             position: %3d\n", interactionMatrix[i].sitePos);
    printf("               strand: %3d\n", interactionMatrix[i].strand);
    printf("         Hamming dist: %3d\n", interactionMatrix[i].hammingDist);
  }
}

void InitializeGenotype(struct Genotype *indiv, 
                        float kdis[])
{
  int i, j;
  
  InitializeSeq((char *)indiv->cisRegSeq,reglen*NGenes);
  InitializeSeq((char *)indiv->transcriptionFactorSeq,elementlen*NGenes); 
  CalcInteractionMatrix(indiv->cisRegSeq,indiv->transcriptionFactorSeq,&(indiv->bindSiteCount),&(indiv->interactionMatrix));
  /* PrintInteractionMatrix(indiv->interactionMatrix, indiv->bindSiteCount, indiv->transcriptionFactorSeq, indiv->cisRegSeq);  */
  fprintf(fperrors,"activators vs repressors ");
  for (i=0; i<NGenes; i++){
    indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    while (indiv->mRNAdecay[i]<0.0)
      indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    indiv->proteindecay[i]=-1.0;
    while (indiv->proteindecay[i]<0.0){
      if (ran1(&seed)<0.08421)
        indiv->proteindecay[i] = 0.0;
      else indiv->proteindecay[i] = exp(0.7874*gasdev(&seed)-3.7665);
    }
    indiv->proteindecay[i] += 0.00578; //dilution due to cell growth
    indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    while (indiv->translation[i] < 0.0)
      indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    if (ran1(&seed)<pact) indiv->activating[i] = 1;
    else indiv->activating[i] = 0;
    fprintf(fperrors,"%d ",indiv->activating[i]);
    j = trunc(Nkdis * ran1(&seed));
    indiv->PICdisassembly[i] = kdis[j];
  }
  fprintf(fperrors,"\n");
}

void Mutate(struct Genotype *old,
            struct Genotype *new,
            float m)
{
  int i,k;
  char x;
  
  for (i=0; i<NGenes; i++){
    for (k=0; k<reglen; k++){
      new->cisRegSeq[i][k] = old->cisRegSeq[i][k];      
      if (m > ran1(&seed)){
        x = old->cisRegSeq[i][k]; //because sometimes old and new are the same
        while (new->cisRegSeq[i][k] == x)
          InitializeSeq(&(new->cisRegSeq[i][k]),(int) 1);
      }
    }
    for (k=0; k<elementlen; k++) 
      new->transcriptionFactorSeq[i][k] = old->transcriptionFactorSeq[i][k];    
    new->mRNAdecay[i] = old->mRNAdecay[i];
    new->proteindecay[i] = old->proteindecay[i];
    new->translation[i] = old->translation[i];
    new->activating[i] = old->activating[i];
    new->PICdisassembly[i] = old->PICdisassembly[i];
  }
  CalcInteractionMatrix(new->cisRegSeq,new->transcriptionFactorSeq,&(new->bindSiteCount),&(new->interactionMatrix));
}


void CalcInteractionMatrix(char cisRegSeq[NGenes][reglen],
                           char transcriptionFactorSeq[NGenes][elementlen],
                           int *newBindSiteCount,
                           struct TFInteractionMatrix **interactionMatrix)
{
  int i,j,geneind,tfind,match,maxy,bindSiteCount;
  
  maxy = maxelements;
  *interactionMatrix = malloc(maxy*sizeof(struct TFInteractionMatrix));
  if (!(*interactionMatrix)){
    fprintf(fperrors,"initial setting of G failed.\n");
    exit(1);
  }
  bindSiteCount = 0;
  for (geneind=0; geneind<NGenes; geneind++) { /* which cis-reg region */
    for (i=0; i<reglen-elementlen; i++) {      /* scan forwards */
      for (tfind=0; tfind<NGenes; tfind++) {
        match=0;
        for (j=i; j<i+elementlen; j++) {
          if (cisRegSeq[geneind][j] == transcriptionFactorSeq[tfind][j-i])
            match++;
        }
        if (match>=nmin){
          if (bindSiteCount+1 >= maxy) {
            maxy = 2*maxy;
            *interactionMatrix = realloc(*interactionMatrix, maxy*sizeof(struct TFInteractionMatrix));
            if (!(*interactionMatrix)) {
              fprintf(fperrors,"realloc of G to bindSiteCount = %d failed.\n",maxy);
              exit(1);
            }
            else if (verbose) fprintf(fperrors,"realloc of G to bindSiteCount = %d succeeded\n",maxy);
          }
          (*interactionMatrix)[bindSiteCount].cisregID = geneind;
          (*interactionMatrix)[bindSiteCount].tfID = tfind;
          (*interactionMatrix)[bindSiteCount].sitePos = i;
          (*interactionMatrix)[bindSiteCount].strand = 0;
          (*interactionMatrix)[bindSiteCount].hammingDist = elementlen-match;
          bindSiteCount++;
        }
      }
    }
    for (i=reglen-1; i>=elementlen-1; i--) {  /* scan backwards */
      for (tfind=0; tfind<NGenes; tfind++) {
        match=0;
        for (j=i; j>i-elementlen; j--)
          if (
             (cisRegSeq[geneind][j]=='a' && transcriptionFactorSeq[tfind][i-j]=='t')
             || (cisRegSeq[geneind][j]=='t' && transcriptionFactorSeq[tfind][i-j]=='a')
             || (cisRegSeq[geneind][j]=='c' && transcriptionFactorSeq[tfind][i-j]=='g')
             || (cisRegSeq[geneind][j]=='g' && transcriptionFactorSeq[tfind][i-j]=='c')            
             ) match++;
        if (match>=nmin){
          if (bindSiteCount+1>=maxy){
            maxy = 2*maxy;
            *interactionMatrix = realloc(*interactionMatrix, maxy*sizeof(struct TFInteractionMatrix));
            if (!(*interactionMatrix)){
              fprintf(fperrors,"realloc of G to bindSiteCount = %d failed.\n",maxy);
              exit(1);
            }
          }
          (*interactionMatrix)[bindSiteCount].cisregID = geneind;
          (*interactionMatrix)[bindSiteCount].tfID = tfind;
          (*interactionMatrix)[bindSiteCount].sitePos = i-elementlen+1;
          (*interactionMatrix)[bindSiteCount].strand = 1;
          (*interactionMatrix)[bindSiteCount].hammingDist = elementlen-match;
          bindSiteCount++;
        }
      }
    }
  }
  *interactionMatrix = realloc(*interactionMatrix, bindSiteCount*5*sizeof(int));
  if (!(*interactionMatrix)) {
    fprintf(fperrors, "realloc of G down to bindSiteCount = %d failed.\n", bindSiteCount);
    exit(1);
  }
  *newBindSiteCount = bindSiteCount;
}

/* transcription/translation delays are sorted linked lists.
   Deleting the head each time, and tack new stuff on the end
   Linked lists are easy to create pre-sorted.
*/

struct FixedEvent
{
  int geneID;
  float time;
  struct FixedEvent *next;
};

struct CellState
{
  int mRNACytoCount[NGenes];          /* mRNAs in cytoplasm */
  int mRNANuclearCount[NGenes];       /* mRNAs in nucleus */
  int mRNATranslCytoCount[NGenes];    /* mRNAs are in the cytoplasm, but only recently */
  struct FixedEvent *mRNATranslTimeEnd;    /* times when mRNAs move to cytoplasm (mRNACytoCount) */
  struct FixedEvent *mRNATranslTimeEndLast; 
  int mRNATranscrCount[NGenes];             /* mRNAs which haven't finished transcription yet */
  struct FixedEvent *mRNATranscrTimeEnd;    /* times when transcripts ? move to nucleus (mRNANuclearCount) */
  struct FixedEvent *mRNATranscrTimeEndLast;
  float proteinConc[NGenes];
  int tfBoundCount;
  int *tfBoundIndexes;
  /* tfBoundIndexes lists just bound TFs according to binding site index in G */
  int tfHinderedCount;
  int (*tfHinderedIndexes)[2];
  /*1st elem tfHinderedIndexes lists binding site indices that cannot be bound due to steric hindrance
    2nd elem gives corresponding index of inhibiting TF in G (TODO: what is this?)
  */
  int active[NGenes];
  /* gives the state of each of the genes, according to figure
     1 is fully off, 2 meets TF criteria
     3 is off but w/o nucleosome, 4 is on but w/o PIC
     5 is on but w/o TF criteria, 6 is fully on
  */
};

void FreeMemCellState(struct CellState *state)
{
  struct FixedEvent *start,*info;

  start = state->mRNATranslTimeEnd;  
  while (start){
    info = start;
    start = start->next;
    free(info);  
  }
  start = state->mRNATranscrTimeEnd;  
  while (start){
    info = start;
    start = start->next;
    free(info);  
  }
  free(state->tfBoundIndexes);
  free(state->tfHinderedIndexes);
}

void sls_store(struct FixedEvent *i, 
               struct FixedEvent **start, 
               struct FixedEvent **last)
{
  struct FixedEvent *old, *p;
  
  p= *start;
  if (!*last) { /*first element in list*/
    i->next = NULL;
    *last = i;
    *start = i;
    return;
  }
  old=NULL;
  while (p) {
    if (p->time < i->time){
      old=p;
      p = p->next;
    }
    else {
      if (old) { /*goes in the middle*/
        old->next = i;
        i->next =p;
        return;
      }
      i->next = p; /*new first element*/
      *start = i;
      return;
    }
  }
  (*last)->next = i; /*put on end*/
  i->next = NULL;
  *last = i;
}

struct TimeCourse
{
  float concentration;
  float time;
  struct TimeCourse *next;
};     

void DeleteTimeCourse(struct TimeCourse *start2)
{
  struct TimeCourse *info,*start;
  
  start = start2; 
  while (start){
    info = start;
    start = start->next;
    free(info);
  }
}

void display2(struct TimeCourse *start)
{
  struct TimeCourse *info;

  info = start;
  while (info){
    fprintf(fperrors,"time %g conc %g\n",info->time,info->concentration);
    info = info->next;
  }
}
      
void sls_store_end(struct FixedEvent *i, 
                   struct FixedEvent **start, 
                   struct FixedEvent **last)
{
  i->next = NULL;
  if (!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

void sls_store_end2(struct TimeCourse *i, 
                    struct TimeCourse **start, 
                    struct TimeCourse **last)
{
  i->next = NULL;
  if (!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

void display(struct FixedEvent *start)
{
  struct FixedEvent *info;

  info = start;
  while (info){
    fprintf(fperrors,"gene %d time %f\n",info->geneID,info->time);
    info = info->next;
  }
}

void AddFixedEvent(int i,
                   float t,
                   struct FixedEvent **start,
                   struct FixedEvent **last)
{
  struct FixedEvent *newtime;
  
  newtime = (struct FixedEvent *)malloc(sizeof(struct FixedEvent));
  if (!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->geneID = i;
  newtime->time = t;
  sls_store(newtime, start, last);
}

void AddTimePoint(float time,
                  float conc,
                  struct TimeCourse **start,
                  struct TimeCourse **last)
{
  struct TimeCourse *newtime;
  
  newtime = (struct TimeCourse *)malloc(sizeof(struct TimeCourse));
  if (!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->time = time;
  newtime->concentration = conc;
  sls_store_end2(newtime, start, last);
}

void AddFixedEventEnd(int i,
                      float t,
                      struct FixedEvent **start,
                      struct FixedEvent **last)
{
  struct FixedEvent *newtime;
  
  newtime = (struct FixedEvent *)malloc(sizeof(struct FixedEvent));
  if (!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->geneID = i;
  newtime->time = t;
  sls_store_end(newtime, start, last);
}

void DeleteFixedEvent(int geneID,
                      int i,
                      struct FixedEvent **start,
                      struct FixedEvent **last)
{
  struct FixedEvent *info,*lastinfo;
  int j,done;
  
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

void DeleteFixedEventStart(struct FixedEvent **start,
                           struct FixedEvent **last)
{
  struct FixedEvent *info;
  
  info = *start;
  *start = info->next;
  if (*last == info) *last = NULL;
  free(info);
}

void InitializeCell(struct CellState *indiv,
                    /*  int y[NGenes]: AKL 2008-03-21: removed wasn't being used */
                    float mRNAdecay[NGenes],
                    float meanmRNA[NGenes],
                    float initProteinConc[NGenes])
{
  int i, k, totalmRNA;
  float t;
  
  indiv->mRNATranscrTimeEnd = indiv->mRNATranscrTimeEndLast = NULL;
  indiv->mRNATranslTimeEnd = indiv->mRNATranslTimeEndLast = NULL;
  indiv->tfBoundCount = 0;  /* initialize with nothing bound */
  indiv->tfHinderedCount = 0;
  indiv->tfBoundIndexes = NULL;
  indiv->tfHinderedIndexes = NULL;
  for (i=0; i<NGenes; i++) {
    indiv->active[i] = ON_WITH_NUCLEOSOME;
    totalmRNA = (int) poidev(meanmRNA[i],&seed);
    indiv->mRNANuclearCount[i] = (int) bnldev(startnucleus, totalmRNA, &seed);
    indiv->mRNACytoCount[i] = totalmRNA - indiv->mRNANuclearCount[i];
    indiv->mRNATranslCytoCount[i] = 0;
    for (k=0; k<indiv->mRNACytoCount[i]; k++) {
      t = expdev(&seed) / mRNAdecay[i];
      if (t < ttranslation) {
        (indiv->mRNACytoCount[i])--;
        (indiv->mRNATranslCytoCount[i])++;
        AddFixedEvent(i, ttranslation-t, &(indiv->mRNATranslTimeEnd), &(indiv->mRNATranslTimeEndLast));
      }
    } 
    indiv->mRNATranscrCount[i] = (int) poidev(meanmRNA[i]*ttranscription*mRNAdecay[i],&seed);
    for (k=0; k < indiv->mRNATranscrCount[i]; k++)
      AddFixedEvent(i, ran1(&seed)*ttranscription, &(indiv->mRNATranscrTimeEnd), &(indiv->mRNATranscrTimeEndLast));
    indiv->proteinConc[i] = initProteinConc[i];
  }
}

/* could perhaps be a little faster with option to skip *df calculation for first 2 calls */
void CalcT (float t, 
            float x, 
            struct GillespieRates *rates,
            struct KonStates *konStates,
            float *f, 
            float *df)
{
  float r, denom, numer, ct, ect;
  int i, j;
  
  r = numer = 0.0;

  /* loop over all genes */
  for (i=0; i < NGenes; i++) {
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

void CalckonRate (float t,
                  struct KonStates *konStates,
                  float *konrate)
{
  float r,ct,ect;
  int i;
  
  r = 0.0;
  for (i=0; i<NGenes; i++){
    if (konStates->nkonsum[i]>0){
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
void ChangemRNACyto(int i,
                    struct Genotype *genes,
                    struct CellState *state,
                    struct GillespieRates *rates,
                    struct KonStates *konStates)
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

void Calckoff(int k,
              struct TFInteractionMatrix *interactionMatrix,
              struct CellState *state,
              float *koff,
              float RTlnKr,
              float temperature)
{
  float Gibbs;  /*free energy in kJ/mol*/
  int posdiff, front, back, i, j;
  
  front = back = 0;
  Gibbs = (((float) interactionMatrix[k].hammingDist)/3.0 - 1.0) * RTlnKr; /* subject to revision of elementlen */
  for (j=0; j < state->tfBoundCount; j++) {
    if (interactionMatrix[k].cisregID==interactionMatrix[state->tfBoundIndexes[j]].cisregID && !(k==state->tfBoundIndexes[j])) {
      posdiff = interactionMatrix[k].sitePos - interactionMatrix[state->tfBoundIndexes[j]].sitePos;
      if (abs(posdiff) < 6) {
        fprintf(fperrors,
                "error: steric hindrance has been breached with site %d %d away from site %d\n",
                k, posdiff, state->tfBoundIndexes[j]);
      }
      if (abs(posdiff) < 20) {
        if (posdiff>0) front++; else back++;
      }
    }
    if ((front) && (back)) 
      j=state->tfBoundCount;
  }
  if (front>0) Gibbs -= cooperativity*RTlnKr/3;
  if (back>0) Gibbs -= cooperativity*RTlnKr/3;
  *koff = NumSitesInGenome*kon*0.25/exp(-Gibbs/(GasConstant*temperature));
  //  fprintf(fperrors,"RTlnKr=%g front=%d back=%d H=%d Gibbs=%g koff=%g\n",RTlnKr,front,back,interactionMatrix[k][4],Gibbs,*koff);
  /* 25% protein in nucleus is implicit in formula above */
}

/* when TF binding changes, adjust cooperativity at neighbouring sites */
void ScanNearby(int indexChanged,
                struct TFInteractionMatrix *interactionMatrix,
                struct CellState *state,
                struct GillespieRates *rates,
                float koffvalues[],
                float RTlnKr,
                float temperature)
{
  int posdiff,j,i;
  float diff;
  
  /* for all the currently bound sites */
  for (j = 0; j < state->tfBoundCount; j++) {
    
    /* we are on the same cisreg sequence and we aren't the same TF binding  */
    if (interactionMatrix[indexChanged].cisregID == interactionMatrix[state->tfBoundIndexes[j]].cisregID && !(indexChanged==state->tfBoundIndexes[j])) {

      /* how close are we on the sequence */
      posdiff = interactionMatrix[indexChanged].sitePos - interactionMatrix[state->tfBoundIndexes[j]].sitePos;
      if (abs(posdiff) < 6) { /* within 6: bad: shouldn't happen */
        fprintf(fperrors,
                "error: steric hindrance 2 has been breached with site %d %d away from site %d\n",
                indexChanged, posdiff, state->tfBoundIndexes[j]);
      }
      if (abs(posdiff) < 20) {  /* within 20, adjust koff */

        /* save old value */
        diff = -koffvalues[j];

        /* recompute koffvalues */
        Calckoff(state->tfBoundIndexes[j], interactionMatrix, state, &(koffvalues[j]), RTlnKr, temperature);

        /* calculating how koff changes  */
        diff += koffvalues[j];

        /* adjust rates by difference */
        rates->koff += diff;
      }
    }
  }
}

void Removekon(int siteID,
               int TFID,
               struct GillespieRates *rates,
               float salphc,
               struct KonStates *konStates,
               float proteinConcTFID)
{
  int i, k;
  
  k = 0;
  /* TODO: check ! */
  /* find the index of the site that binds in konIDs  */
  while (!(konStates->konIDs[k][SITEID_INDEX] == siteID) && k < konStates->nkon) {
    k++;
  }

  /* make sure that we have enough unoccupied sites left */
  if (k < konStates->nkon) {   
    /* adjust rates */
    rates->salphc -= kon*salphc;
    rates->maxSalphc -= kon*fmaxf(proteinConcTFID, salphc);
    rates->minSalphc -= kon*fminf(proteinConcTFID, salphc);

    /* one less site available for binding of total */
    (konStates->nkon)--;

    /* also per gene */
    (konStates->nkonsum[TFID])--;

    /* TODO: check */
    /* move the last element end of array into space vacated by site k */
    konStates->konIDs[k][SITEID_INDEX] = konStates->konIDs[konStates->nkon][SITEID_INDEX];
    konStates->konIDs[k][TFID_INDEX] = konStates->konIDs[konStates->nkon][TFID_INDEX];
  }
  //else do nothing: there is likely a redundancy in steric hindrance, hence no site to remove
}

void Addkon(float proteinConcTFID,
            float salphc,
            int TFID,
            int siteID,
            struct GillespieRates *rates,
            struct KonStates *konStates)
{

  /* update rates because new site is now available */
  rates->salphc += kon*salphc;
  rates->maxSalphc += fmaxf(proteinConcTFID, salphc);
  rates->minSalphc += fminf(proteinConcTFID, salphc);

  /* add back siteID to pool of available sites */
  konStates->konIDs[konStates->nkon][SITEID_INDEX] = siteID;
  konStates->konIDs[konStates->nkon][TFID_INDEX] = TFID;

  /* one more site available */
  (konStates->nkonsum[TFID])++;
  (konStates->nkon)++;
}

// tests whether criterion for transcription is met
int CalcTranscription(int geneID,
                      int *tfBoundIndexes,
                      int tfBoundCount,
                      struct TFInteractionMatrix *interactionMatrix,
                      int activating[],
                      int *on)
{
  int i,off;
  
  *on=off=0;
  for (i=0; i < tfBoundCount; i++) {
    if (geneID==interactionMatrix[tfBoundIndexes[i]].cisregID) {
      if (activating[interactionMatrix[tfBoundIndexes[i]].tfID]) (*on)++;
      else off++;
    }
  }
  if ((float)off <= 0.33442*(float)(*on) + 0.31303) 
    return(1);
  else 
    return(0);
}

int IsOneActivator(int geneID,
                   int *tfBoundIndexes,
                   int tfBoundCount,
                   struct TFInteractionMatrix *interactionMatrix,
                   int activating[])
{
  int i;
  
  for (i=0;i<tfBoundCount;i++)
    if (geneID==interactionMatrix[tfBoundIndexes[i]].cisregID && (activating[interactionMatrix[tfBoundIndexes[i]].tfID])) return(1);
  return(0);
}

/* only appropriate if nothing is  bound ie CalcFromInitialState and everything is in state 2
   v31 and earlier have some parts of code needed with prebound stuff,
   rates->mRNAdecay and [7] are done in CalcDt*/
void CalcFromState(struct Genotype *genes,
                   struct CellState *state,
                   struct GillespieRates *rates,
                   struct KonStates *konStates,
                   float transport[],
                   float mRNAdecay[],
                   float RTlnKr,
                   float temperature,
                   int statechangeIDs[][NGenes])
/* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit */
{
  int i, k;
  float salphc, proteinConcTFID; 

  for (i=0; i<NGenes; i++) {
    salphc = (float) (state->mRNACytoCount[i]) * genes->translation[i] / genes->proteindecay[i];
    konStates->konvalues[i][KON_DIFF_INDEX] = (state->proteinConc[i] - salphc) / genes->proteindecay[i];
    konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX] = genes->proteindecay[i];
    konStates->konvalues[i][KON_SALPHC_INDEX] = salphc;
    konStates->nkonsum[i]=0;  
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
    i = genes->interactionMatrix[k].tfID;
    proteinConcTFID = state->proteinConc[i];
    salphc = konStates->konvalues[i][KON_SALPHC_INDEX];
    rates->salphc += salphc;
    rates->maxSalphc += fmaxf(proteinConcTFID, salphc);
    rates->minSalphc += fminf(proteinConcTFID, salphc);
    konStates->konIDs[k][SITEID_INDEX] = k;
    konStates->konIDs[k][TFID_INDEX] = i;
    (konStates->nkonsum[i])++;
  }

  /* initialize konStates->nkon as the total number of binding sites */
  konStates->nkon = genes->bindSiteCount;

  for (i=0; i<NGenes; i++) {
    transport[i] = kRNA * (float) (state->mRNANuclearCount[i]);
    rates->transport += transport[i];
    statechangeIDs[0][i] = i;
  }
  /* for (i=4; i<7; i++) rates[i] *= kon; */
  rates->salphc *= kon;
  rates->maxSalphc *= kon;
  rates->minSalphc *= kon;

  rates->acetylationCount=NGenes;

  /* for (i=1; i<5; i++) rates2[i]=0; */
  rates->deacetylationCount=0;
  rates->picAssemblyCount=0;
  rates->transcriptInitCount=0;
  rates->picDisassemblyCount=0;
}

/* returns:
 *  0 if there is no fixed event occuring before time t
 *  1 if a transcription event happens before time t
 *  2 if a translation event happens before time t
 */
int DoesFixedEventEnd(struct FixedEvent *mRNATranslTimeEnd,
                      struct FixedEvent *mRNATranscrTimeEnd,
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

void CalcDt(float *x,
            float *dt,
            struct GillespieRates *rates,
            struct KonStates *konStates,
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
  for (i=0; i<NGenes; i++){
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
  rates->total += (float) rates->acetylationCount * acetylate;
  rates->total += (float) rates->deacetylationCount * deacetylate;
  rates->total += (float) rates->picAssemblyCount * PICassembly;
  rates->total += (float) rates->transcriptInitCount * transcriptinit;    

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
    *dt = rtsafe(&CalcT, *x, rates, konStates, tbound1, tbound2, (float) 1e-6); 
  }
}

void EndTranscription(float *dt,
                      float t,
                      struct CellState *state,
                      float transport[NGenes],
                      struct GillespieRates *rates)
{
  int i, total;
  
  /* recompute the delta-t based on difference between now and the
     time of transcription end */
  *dt = state->mRNATranscrTimeEnd->time - t;
  total = 0;
  for (i=0; i < NGenes; i++) 
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
  DeleteFixedEventStart(&(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));

  /* add rate kRNA to transport and Gillespie rates
   * TODO: check */
  transport[i] += kRNA;
  rates->transport += kRNA;
}

void TransportEvent(float x,
                    float transport[NGenes],
                    struct CellState *state,
                    float endtime,
                    struct GillespieRates *rates)
{
  int i;
  float konrate2;
  
  i = -1;
  konrate2 = 0.0;  /* TODO-DONE: check, why is this set to zero?, should be set to zero fixed */

  /* choose gene product (mRNA) that gets transported to cytoplasm
     based on weighting in transport[] array */
  while (i < NGenes && x > konrate2){
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
  AddFixedEventEnd(i, endtime, &(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));

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

void RemoveFromArray(int toberemoved,
                     int a[],
                     int *len,
                     int force)
{
  int i;
  
  i=0;
  while (!(a[i]==toberemoved) && i < *len) i++;
  if (i < *len){  
    (*len)--;
    a[i]=a[*len];
  }
  else if (force) fprintf(fperrors, "error removing %d from array of length %d\n", toberemoved, *len);
  // don't always print because with a 4->3 transition PIC assembly is not there to be removed
}

void DisassemblePIC(int *activestate,
                    int geneID,
                    struct GillespieRates *rates,
                    int statechangeIDs[][NGenes],
                    float disassembly)
{
  RemoveFromArray(geneID, statechangeIDs[3], &(rates->transcriptInitCount), (int) 1);
  RemoveFromArray(geneID, statechangeIDs[4], &(rates->picDisassemblyCount), (int) 1);
  rates->picDisassembly -= disassembly;
  
  /* disassemble PIC in OFF state */
  if (*activestate == OFF_PIC) {
    (*activestate) = OFF_NO_PIC;
    statechangeIDs[1][rates->deacetylationCount] = geneID;
    (rates->deacetylationCount)++;
  }
  /* disassemble PIC in ON state */
  if (*activestate == ON_FULL) {
    (*activestate) = ON_NO_PIC;
  }
}

void ReviseActivityState(int geneID,
                         struct Genotype *genes,
                         struct CellState *state,
                         struct GillespieRates *rates,
                         int statechangeIDs[][NGenes])
{
  int transcriptrule, oldstate, numactive;

  transcriptrule = CalcTranscription(geneID, state->tfBoundIndexes, 
                                     state->tfBoundCount,
                                     genes->interactionMatrix,
                                     genes->activating,
                                     &numactive);

  /* get last state of transcription initiation */
  oldstate = state->active[geneID];

  /* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit, 4-PIC disassembly*/

  /*
   * first set of rules:
   * ACTIVATING TFs exceed REPRESSING TFs 
   */

  /* OFF_FULL -> ON_WITH_NUCLEOSOME */
  if ((transcriptrule) && oldstate==OFF_FULL){
    state->active[geneID] = ON_WITH_NUCLEOSOME;
    statechangeIDs[0][rates->acetylationCount] = geneID;
    (rates->acetylationCount)++;      
  }
  /* OFF_NO_PIC -> ON_NO_PIC */
  if ((transcriptrule) && oldstate==OFF_NO_PIC) {
    state->active[geneID] = ON_NO_PIC;
    RemoveFromArray(geneID, statechangeIDs[1], &(rates->deacetylationCount), (int) 1);
    if (numactive){
      statechangeIDs[2][rates->picAssemblyCount] = geneID;
      (rates->picAssemblyCount)++;
    }
  }
  /* OFF_PIC -> ON_FULL */
  if ((transcriptrule) && oldstate==OFF_PIC) {
    state->active[geneID] = ON_FULL;
  }

  /*
   * second set of rules:
   * REPRESSING TFs exceed ACTIVATING TFs 
   */

  /* ON_WITH_NUCLEOSOME -> OFF_FULL */
  if (!(transcriptrule) && oldstate==ON_WITH_NUCLEOSOME) {
    state->active[geneID] = OFF_FULL;
    RemoveFromArray(geneID, statechangeIDs[0], &(rates->acetylationCount), (int) 1);
  }
  
  /* ON_NO_PIC -> OFF_NO_PIC */
  if (!(transcriptrule) && oldstate==ON_NO_PIC){          
    state->active[geneID] = OFF_NO_PIC;
    RemoveFromArray(geneID, statechangeIDs[2], &(rates->picAssemblyCount), (int) 0);
    statechangeIDs[1][rates->deacetylationCount] = geneID;
    (rates->deacetylationCount)++;
  }
  /* ON_FULL -> OFF_PIC  */
  if (!(transcriptrule) && oldstate==ON_FULL) {
    state->active[geneID] = OFF_PIC;
  }

  /* do remaining transitions:
   * OFF_PIC -> OFF_NO_PIC
   * ON_FULL -> ON_NO_PIC 
   */
  if ((state->active[geneID]==OFF_PIC || state->active[geneID]==ON_FULL) && numactive==0)
    DisassemblePIC(&(state->active[geneID]), geneID, rates, statechangeIDs,
                   genes->PICdisassembly[geneID]);

  if (verbose && (oldstate!=state->active[geneID])) {
    fprintf(fperrors, "state change from %d to %d in gene %d\n", oldstate, state->active[geneID], geneID);
  }
}

void RemoveBinding(struct Genotype *genes,
                   struct CellState *state,
                   struct GillespieRates *rates,
                   struct KonStates *konStates,
                   int site,
                   float koffvalues[],
                   float RTlnKr,
                   float temperature,
                   int statechangeIDs[][NGenes])
{
  int i, j, k, bound, siteID, geneID, transcriptrule, oldstate, numactive;

  i = 0;

  /* given site 'site', look for the index in the list of bound sites */
  while ((state->tfBoundIndexes[i] != site) && (i < state->tfBoundCount)) 
    i++;
  if (i == state->tfBoundCount) {  /* couldn't find the site */
    fprintf(fperrors,"error: RemoveBinding could not find site %d with %d possibilities\n Bound sites are\n",
            site,state->tfBoundCount);
    for (j = 0; j < state->tfBoundCount; j++)  {
      fprintf(fperrors,"%d\n",state->tfBoundIndexes[j]);
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
            fprintf(fperrors,"Site %d pos %d on gene %d freed from steric hindrance\n",
                    siteID, genes->interactionMatrix[siteID].sitePos, genes->interactionMatrix[siteID].cisregID);

          /* adjust rates by returning kon to pool */
          Addkon(state->proteinConc[genes->interactionMatrix[siteID].tfID],
                 konStates->konvalues[genes->interactionMatrix[siteID].tfID][KON_SALPHC_INDEX],
                 genes->interactionMatrix[siteID].tfID,
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

    /* find the gene whose cisreg region has an unbinding event */
    geneID = genes->interactionMatrix[site].cisregID;
    if (verbose) 
      fprintf(fperrors,"Add site %d at pos %d on gene %d freed by unbinding\n",
              site, genes->interactionMatrix[site].sitePos, geneID);

    /* adjust kon */
    Addkon(state->proteinConc[genes->interactionMatrix[site].tfID],
           konStates->konvalues[genes->interactionMatrix[site].tfID][KON_SALPHC_INDEX],
           genes->interactionMatrix[site].tfID,
           site,
           rates,
           konStates);

    /* adjust the state of the gene */
    ReviseActivityState(geneID, genes, state, rates, statechangeIDs);

    /* when TF unbinds adjust the co-operativity at close sites */
    ScanNearby(site, genes->interactionMatrix, state, rates, koffvalues, RTlnKr, temperature);        
  }
}

void TFbinds(struct Genotype *genes,
             struct CellState *state,
             struct GillespieRates *rates,
             float **koffvalues,
             struct KonStates *konStates,
             int *maxbound2,
             int *maxbound3,
             int site,
             float RTlnKr,
             float temperature,
             int statechangeIDs[][NGenes])
{
  int geneID, k, posdiff, site2;

  if (verbose) fprintf(fperrors,
                      "kon1 event at site %d out of %d possible, %d TFs previously bound bindSiteCount=%d\n",
                      site, konStates->nkon, state->tfBoundCount, genes->bindSiteCount);

  /* if we have run out of space, double memory  */
  if (state->tfBoundCount >= *maxbound2){
    (*maxbound2) *= 2;
    state->tfBoundIndexes = realloc(state->tfBoundIndexes,(*maxbound2)*sizeof(int));

    /* do the copy */
    *koffvalues = realloc(*koffvalues,(*maxbound2)*sizeof(float));

    /* check return value */
    if (!state->tfBoundIndexes || !(*koffvalues)){
      fprintf(fperrors,"memory allocation error resetting maxbound2=%d\n",*maxbound2);
      exit(1);
    }
  }

  /* append the site to end of indexes */
  state->tfBoundIndexes[state->tfBoundCount] = site;
  if (verbose) fprintf(fperrors,"remove site %d\n", site);

  /* remove the site from the kon pool */
  Removekon(site,
            genes->interactionMatrix[site].tfID,
            rates, 
            konStates->konvalues[genes->interactionMatrix[site].tfID][KON_SALPHC_INDEX],
            konStates,
            state->proteinConc[genes->interactionMatrix[site].tfID]);

  /* recompute the koffvalues */
  Calckoff(site, genes->interactionMatrix, state, &((*koffvalues)[state->tfBoundCount]), RTlnKr, temperature);
  if (verbose) fprintf(fperrors,"new koff = %g is number %d\n",
                       (*koffvalues)[state->tfBoundCount],(state->tfBoundCount+1));

  /* adjust rates by adding the new koffvalue to rates->koff */
  rates->koff += (*koffvalues)[state->tfBoundCount];

  /* TODO: check, this seems the same as above ???  maybe remove?? */
  state->tfBoundIndexes[state->tfBoundCount] = site;

  /* increment number of bound TFs */
  (state->tfBoundCount)++;

  /* adjust co-operative binding in context of new TF */
  ScanNearby(site, genes->interactionMatrix, state, rates, *koffvalues, RTlnKr, temperature);

  /* get the gene that the TF is binding to */
  geneID = genes->interactionMatrix[site].cisregID;

  /* update steric hindrance data structures */
  /* JM: this cycles over all sites, not just bound ones, in order to
     record redundancy in steric hindrance*/
  for (k = 0; k < genes->bindSiteCount; k++) {

    /* if we are on the same gene and not the same binding site */
    if (geneID == genes->interactionMatrix[k].cisregID && !(k==site)) {

      /* check distance from current binding site (k) to the original (site) */
      posdiff = genes->interactionMatrix[site].sitePos - genes->interactionMatrix[k].sitePos;

      /* if within 6, we prevent binding by adding to steric hindrance */
      if (abs(posdiff) < 6) {

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
          fprintf(fperrors,"%d steric hindrance sites after %d blocks site %d\n", state->tfHinderedCount, site, k);
        }

        /* remove the kon from pool */
        Removekon(k,
                  genes->interactionMatrix[k].tfID,
                  rates,
                  konStates->konvalues[genes->interactionMatrix[k].tfID][KON_SALPHC_INDEX],
                  konStates,
                  state->proteinConc[genes->interactionMatrix[k].tfID]);
      }
    }
  }
  if (verbose) fprintf(fperrors,"tfBoundCount=%d tfHinderedCount=%d maxbound2=%d maxbound3=%d\n",
                      state->tfBoundCount,state->tfHinderedCount,*maxbound2,*maxbound3);

  /* gene activity may change as a result of binding */
  ReviseActivityState(geneID, genes, state, rates, statechangeIDs);
}

/*time course of [TF]s represented as array of TimeCourse lists.
 */

void AddTimePoints(float time,
                   float proteinConc[NGenes],
                   struct TimeCourse **timecoursestart,
                   struct TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i<NGenes; i++)
    AddTimePoint(time, proteinConc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

void AddIntTimePoints(float time,
                      int proteinConc[NGenes],
                      struct TimeCourse **timecoursestart,
                      struct TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i<NGenes; i++)
    AddTimePoint(time, (float) proteinConc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

/* need some sort of control in case it declines to essentially zero.
   Add in discrete, stochastic and/or zero values, but this may create false attractor
   if time steps are small and rising tide keeps getting rounded down*/
void UpdateProteinConc(float proteinConc[],
                       float dt,
                       struct GillespieRates *rates,
                       struct KonStates *konStates,
                       float t,
                       struct TimeCourse **timecoursestart,
                       struct TimeCourse **timecourselast,
                       float otherdata[])
{
  int i;
  float ct, ect, ect1;
  
  rates->maxSalphc = rates->minSalphc = 0.0;
  for (i=0; i<NGenes; i++) {
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
    AddTimePoints(t+dt,otherdata/*proteinConc*/,timecoursestart,timecourselast);
}

void CalcNumBound(float proteinConc[],
                  int tfBoundCount)
{
  float sum;
  int i;
  
  sum = 0.0;
  for (i=0; i < NGenes; i++) 
    sum += proteinConc[i];
  if (verbose) 
    /* fprintf(fperrors, "%d bound %g expected\n", tfBoundCount, 0.0003*sum);*/
    /* TODO: check! */
    /* if this is wrong for random sequences adjust Kr accordingly */
    fprintf(fperrors, "%d bound %g expected\n", tfBoundCount, (reglen*NGenes*sum)/NumSitesInGenome);
}

void Develop(struct Genotype *genes,
             struct CellState *state,
             float temperature, /* in Kelvin */
             struct TimeCourse **timecoursestart,
             struct TimeCourse **timecourselast)
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
  float transport[NGenes];  /* transport rates of each mRNA */
  float mRNAdecay[NGenes];  /* mRNA decay rates */
  float x;                  /* random number */
  float dt;                 /* delta-t */

  /* cached information about available binding sites for efficiency */
  struct KonStates *konStates =  malloc(sizeof(struct KonStates));

  struct GillespieRates *rates = malloc(sizeof(struct GillespieRates));

  /* stores corresponding geneIDs for [de]acteylation, PIC[dis]assembly, transcriptinit */
  int statechangeIDs[5][NGenes]; 

  float f, df, konrate, konrate2, diff, RTlnKr, sum, ct, ect;

  /* initialize time courses */
  for (i=0; i<NGenes; i++){
    timecoursestart[i] = NULL;
    timecourselast[i] = NULL;
  }
  AddTimePoints((float) 0.0,state->proteinConc,timecoursestart,timecourselast);

  RTlnKr = GasConstant * temperature * log(Kr);     /* compute constant */

  /* number of possible binding sites */
  konStates->konIDs = malloc(2*genes->bindSiteCount*sizeof(int));

  /* TODO: ? */
  maxbound2 = maxbound;
  maxbound3 = 10*maxbound;

  state->tfBoundIndexes = realloc(state->tfBoundIndexes,maxbound2*sizeof(int));
  koffvalues = malloc(maxbound2*sizeof(float));
  state->tfHinderedIndexes = realloc(state->tfHinderedIndexes,2*maxbound3*sizeof(int));
  if (!konStates->konvalues || !state->tfBoundIndexes || !state->tfHinderedIndexes || !konStates->konIDs){
    fprintf(fperrors,"memory allocation error at start of Develop\n");
    exit(1);
  }

  /* initialize transcriptional state of genes */
  CalcFromState(genes, state, /*&nkon, nkonsum, */ rates, konStates, /* konvalues,
                                                                        konIDs, */
                transport, mRNAdecay, RTlnKr, temperature,
                statechangeIDs);

  t=0.0;  /* time starts at zero */

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

    CalcDt(&x, &dt, /*nkon, nkonsum, */ rates, konStates, /*konvalues,*/ mRNAdecay, genes->mRNAdecay,
           state->mRNACytoCount, state->mRNATranslCytoCount);

    if (verbose) 
      fprintf(fperrors,"next stochastic event due at t=%g dt=%g x=%g\n", t+dt, dt, x);

    if (!(state->mRNATranscrTimeEndLast)){
      for (i=0;i<NGenes;i++)
        if (verbose) fprintf(fperrors,"%d transcription events\n", state->mRNATranscrCount[i]);
    }
    
    /* first check to see if a fixed event occurs in current t->dt window */
    event=DoesFixedEventEnd(state->mRNATranslTimeEnd,
                            state->mRNATranscrTimeEnd,
                            fminf(t+dt, tdevelopment));

    /* while there are either transcription or translation events
       occuring in current t->dt window */
    while (event > 0) {
      konrate = x/dt;

      if (event==1) {  /* if a transcription event ends */
        EndTranscription(&dt, t, state, transport, rates);

        UpdateProteinConc(state->proteinConc, dt,
                          rates, konStates, t,
                          timecoursestart, timecourselast,
                          state->proteinConc);
      } 

      else {           /* if a translation event ends */
                       /* TODO: could have else if (event==2) */
        dt = state->mRNATranslTimeEnd->time - t;         /* make dt window smaller */
        total=0;  /* number of translation events */

        /* count current number of mRNAs that have recently arrived in cytoplasm */
        for (i=0;i<NGenes;i++) total += state->mRNATranslCytoCount[i];
        if (verbose) fprintf(fperrors,"\ntranslation event finishes out of %d possible t=%g dt=%g\n",
                             total, t, dt); // bug: dt can be negative

        /* get identity of gene that has just finished translating */
        i=state->mRNATranslTimeEnd->geneID;   

        /* there is one less mRNA that has just finished translation */
        (state->mRNATranslCytoCount[i])--;   

        /* delete the event that just happened */
        DeleteFixedEventStart(&(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));

        /* there is one more mRNA that is now in cytoplasm */
        (state->mRNACytoCount[i])++;

        /* update protein concentration */
        UpdateProteinConc(state->proteinConc, dt,
                          rates, konStates, t,
                          timecoursestart, timecourselast,
                          state->proteinConc);
        
        /* the number of mRNAs in cytoplasm affects binding, TODO: check*/
        ChangemRNACyto(i, genes, state, rates, konStates);
      }

      /* advance time by the dt */
      t += dt;
      x -= dt*konrate;
      if (verbose) 
        fprintf(fperrors, "dt=%g t=%g fixed event old x=%g new x=%g\n", dt, t, x+dt*konrate, x);

      /* compute a new dt */
      CalcDt(&x, &dt, rates, konStates, mRNAdecay, 
             genes->mRNAdecay, state->mRNACytoCount, state->mRNATranslCytoCount);

      if (verbose) 
        fprintf(fperrors,"next stochastic event (2) due at t=%g dt=%g x=%g\n",t+dt,dt,x);

      /* check to see there aren't more fixed events to do */
      event=DoesFixedEventEnd(state->mRNATranslTimeEnd, state->mRNATranscrTimeEnd, fminf(tdevelopment,t+dt));
    } 

    /* no remaining fixed events to do in dt, now do stochastic events */

    /* if we haven't already reached end of development with last delta-t */
    if (t+dt < tdevelopment) {  

      /* compute total konrate (which is constant over the Gillespie step) */
      if (konStates->nkon==0) konrate = (-rates->salphc);
      else CalckonRate(dt, konStates, &konrate); 

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
        fprintf(fperrors,"acetylation=%g\ndeacetylation=%g\nPIC assembly=%g\ntranscriptinit=%g\n",
                (float)rates->acetylationCount*acetylate, (float)rates->deacetylationCount*deacetylate, 
                (float)rates->picAssemblyCount*PICassembly, (float)rates->transcriptInitCount*transcriptinit);
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
        j = -1;
        while (j<state->tfBoundCount && x>0) {
          j++;
          x -= koffvalues[j];
        }
        if (j==state->tfBoundCount) {
          konrate2 = 0.0;
          for (i = 0; i < state->tfBoundCount; i++) 
            konrate2 += koffvalues[i];
          fprintf(fperrors, "warning: koffvalues add up to %g instead of rates->koff=%g\n",
                  konrate2, rates->koff);
          rates->koff = konrate2;
          j--; // a bit of a fudge for rounding error, really should move on to rates->transport, but too complicated for something so minor
        } 
        site=state->tfBoundIndexes[j];
        if (verbose) fprintf(fperrors, "koff event %d of %d at site %d\n",
                             j,state->tfBoundCount,site);
        if (j<0) fprintf(fperrors, "error: koff event %d of %d at site %d\n", j, state->tfBoundCount,site);

        /* update protein concentration before removing TF */
        UpdateProteinConc(state->proteinConc, dt, 
                          rates, konStates, 
                          t, timecoursestart, timecourselast, 
                          state->proteinConc);

        /* remove TF binding from 'site' */
        RemoveBinding(genes, state, rates, konStates, site,
                      koffvalues, RTlnKr, temperature, statechangeIDs);
        CalcNumBound(state->proteinConc, state->tfBoundCount);
      } else {
        x -= rates->koff;  
        
        /* 
         * STOCHASTIC EVENT: a transport event
         */
        if (x < rates->transport) {     
          UpdateProteinConc(state->proteinConc, dt, 
                            rates, konStates, 
                            t, timecoursestart, timecourselast, 
                            state->proteinConc);
          TransportEvent(x, transport, state, t+dt+ttranslation, rates);
        } else {
          x -= rates->transport;
          /* 
           * STOCHASTIC EVENT: an mRNA decay event
           */
          if (x<rates->mRNAdecay) {  
            i = -1;
            konrate2 = 0.0;

            /* loop through mRNA products, to choose the mRNA with the
               proportionally higher decay rate */
            while (i < NGenes-1 && x > konrate2) {
              i++;
              konrate2 += mRNAdecay[i];
            }
            if (x > konrate2) { /* JM: had some rounding errors with rates->mRNAdecay. Calculate in CalcDt, hopefully fixed now */
              fprintf(fperrors, "warning: x=%g > konrate2=%g out of rates->mRNAdecay=%g\n",
                      x, konrate2, rates->mRNAdecay);
            }

            /* assume mRNA cytoplasm transport events equally likely */
            x = ran1(&seed)*((float) (state->mRNACytoCount[i] + state->mRNATranslCytoCount[i]));
            UpdateProteinConc(state->proteinConc, dt,
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
              ChangemRNACyto(i, genes, state, rates, konStates); 
            
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
              DeleteFixedEvent(i, (int) trunc(x) ,&(state->mRNATranslTimeEnd), &(state->mRNATranslTimeEndLast));
              
              /* remove the mRNA from the count */
              (state->mRNATranslCytoCount[i])--;
              if (verbose) for (j=0; j < NGenes; j++)
                            fprintf(fperrors, "%d copies of gene %d translating\n", 
                                    state->mRNATranslCytoCount[j], j);                  
            }
          } else {
            x -= rates->mRNAdecay;
            /* 
             * STOCHASTIC EVENT: PIC disassembly
             */
            if (x < rates->picDisassembly) {             
              j=-1;
              while (j<NGenes && x>0){
                j++;
                x -= genes->PICdisassembly[statechangeIDs[4][j]];                
              }
              if (j==NGenes) fprintf(fperrors, "error in PIC disassembly\n");
              j=statechangeIDs[4][j];
              if (verbose) fprintf(fperrors, "PIC disassembly event in gene %d\n", j);
              DisassemblePIC(&(state->active[j]), j, rates, statechangeIDs,
                             genes->PICdisassembly[j]);
            } else {
              x -= rates->picDisassembly;
              /* 
               * STOCHASTIC EVENT: TF binding event
               */
              if (x <rates->salphc + konrate) {   /* add variable (salphc) and constant (konrate) */

                x = ran1(&seed) * (rates->salphc + konrate)/kon;
                j = -1;
                konrate2 = 0.0;               

                /* loop through all *available* binding sites to
                 * choose one, the interval is weighted by the
                 * frequency of rates */
                while (j < konStates->nkon-1 && x > konrate2) {

                  /* this will record the last binding site before we
                     reach the appropriate binding site  */
                  j++;

                  /* get ID of TF */
                  i = konStates->konIDs[j][TFID_INDEX];

                  /* compute the kon rate  */
                  konrate2 = konStates->konvalues[i][KON_SALPHC_INDEX] + 
                    konStates->konvalues[i][KON_DIFF_INDEX]*(1-exp(-konStates->konvalues[i][KON_PROTEIN_DECAY_INDEX]*dt))/dt;

                  /* adjust random number */
                  x -= konrate2;
                }
                
                /* update protein concentration before doing the binding */
                UpdateProteinConc(state->proteinConc, dt, 
                                  rates, konStates, 
                                  t, timecoursestart, timecourselast,
                                  state->proteinConc);
                if (verbose) fflush(fperrors);

                /* bind the j-th konID in the list which happens to be the binding site in SITEID_INDEX */
                TFbinds(genes, state, rates, &koffvalues, konStates, &maxbound2, &maxbound3, 
                        konStates->konIDs[j][SITEID_INDEX], RTlnKr, temperature, statechangeIDs);

                /* calculate the number of TFs bound
                 * TODO: check this seems to be purely diagnostic
                 */
                CalcNumBound(state->proteinConc, state->tfBoundCount);
              } else {
                x -= (rates->salphc + konrate);
                /* 
                 * STOCHASTIC EVENT: histone acetylation
                 */
                if (x < (float) rates->acetylationCount * acetylate) {
                  x = ran1(&seed)*((float) rates->acetylationCount);

                  /* choose a particular gene to change state */
                  geneID = statechangeIDs[0][(int)trunc(x)];
                  if (verbose) fprintf(fperrors,"acetylation event gene %d\nstate change from %d to 4\n",
                                       geneID, state->active[geneID]);
                  if (state->active[geneID] != ON_WITH_NUCLEOSOME)
                    fprintf(fperrors, "error: acetylation event attempted from state %d\n", state->active[geneID]);
                  UpdateProteinConc(state->proteinConc,dt,
                                    rates, konStates, 
                                    t, timecoursestart, timecourselast,
                                    state->proteinConc);

                  /* set state: eject nucleosome, but there is no PIC yet */
                  state->active[geneID] = ON_NO_PIC;
                  RemoveFromArray(geneID, statechangeIDs[0], &(rates->acetylationCount), (int) 1);
                  if (IsOneActivator(geneID, state->tfBoundIndexes, state->tfBoundCount, 
                                     genes->interactionMatrix, genes->activating)) {
                    statechangeIDs[2][rates->picAssemblyCount] = geneID; 
                    (rates->picAssemblyCount)++;
                  }
                } else {
                  x -= (float) rates->acetylationCount * acetylate;
                  /* 
                   * STOCHASTIC EVENT: histone deacetylation
                   */
                  if (x < (float) rates->deacetylationCount * deacetylate) {
                    x = ran1(&seed)*((float) rates->deacetylationCount);
                    geneID = statechangeIDs[1][(int)trunc(x)];
                    if (verbose) fprintf(fperrors,"deacetylation event gene %d\nstate change from %d to 1\n",
                                         geneID,state->active[geneID]);
                    if (state->active[geneID] != OFF_NO_PIC)
                      fprintf(fperrors, "error: deacetylation event attempted from state %d\n", state->active[geneID]);
                    UpdateProteinConc(state->proteinConc, dt, 
                                      rates, konStates, 
                                      t, timecoursestart, timecourselast,
                                      state->proteinConc);
                    /* set state: nucleosome returns */
                    state->active[geneID] = OFF_FULL;
                    RemoveFromArray(geneID, statechangeIDs[1], &(rates->deacetylationCount), (int) 1);
                  } else {
                    x -= (float) rates->deacetylationCount * deacetylate;
                    /* 
                     * STOCHASTIC EVENT: PIC assembly
                     */
                    if (x < (float) rates->picAssemblyCount * PICassembly) {
                      x = ran1(&seed)*((float) rates->picAssemblyCount);
                      geneID = statechangeIDs[2][(int)trunc(x)];
                      if (verbose) fprintf(fperrors,"PIC assembly event gene %d\nstate change from %d to 6\n",
                                          geneID,state->active[geneID]);
                      if (state->active[geneID] != ON_NO_PIC)
                        fprintf(fperrors, "error: PIC assembly event attempted from state %d\n", state->active[geneID]);
                      UpdateProteinConc(state->proteinConc, dt,
                                        rates, konStates, 
                                        t, timecoursestart, timecourselast,
                                        state->proteinConc);

                      /* turn gene fully on: ready for transcription */
                      state->active[geneID] = ON_FULL;
                      RemoveFromArray(geneID, statechangeIDs[2], &(rates->picAssemblyCount), (int) 1);                      
                      statechangeIDs[3][rates->transcriptInitCount] = geneID;
                      (rates->transcriptInitCount)++;
                      statechangeIDs[4][rates->picDisassemblyCount] = geneID;
                      (rates->picDisassemblyCount)++;
                      rates->picDisassembly += genes->PICdisassembly[geneID];                                            
                    } else {
                      x -= (float) rates->picAssemblyCount * PICassembly;
                      /* 
                       * STOCHASTIC EVENT: transcription initiation
                       */
                      if (x < (float) rates->transcriptInitCount * transcriptinit) {
                        x /= transcriptinit;
                        geneID = statechangeIDs[3][(int)trunc(x)];
                        if (verbose) fprintf(fperrors,"transcription event gene %d\n",geneID);
                        if (state->active[geneID] != ON_FULL && state->active[geneID] != OFF_PIC)
                          fprintf(fperrors,"error: transcription event attempted from state %d\n", state->active[geneID]);
                        UpdateProteinConc(state->proteinConc,dt,
                                          rates, konStates, 
                                          t, timecoursestart, timecourselast,
                                          state->proteinConc);

                        /* now that transcription of gene has been initiated, 
                         * we add the time it will end transcription, 
                         * which is dt+time of transcription from now */
                        AddFixedEventEnd(geneID, t+dt+ttranscription, 
                                         &(state->mRNATranscrTimeEnd), &(state->mRNATranscrTimeEndLast));

                        /* increase the number mRNAs being transcribed */
                        (state->mRNATranscrCount[geneID])++;                      
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
      /* Gillespie step: advance time to next event at dt */
      t += dt;
      if (verbose) fprintf(fperrors, "dt=%g t=%g\n", dt, t);
    } else {
      /* we will reach the end of development in dt */
      if (verbose) fprintf(fperrors, "finish at t=%g dt=%g\n", t, dt);

      /* do remaining dt */
      dt = tdevelopment - t;

      /* final update of protein concentration */
      UpdateProteinConc(state->proteinConc, dt,
                        rates, konStates,
                        t, timecoursestart, timecourselast,
                        state->proteinConc);
      /* advance to end of development (this exits the outer while loop) */
      t = tdevelopment;
    }
  }
  free(koffvalues);
  free(konStates->konIDs);
  free(konStates);
  free(rates);
}

void DevStabilityOnlyLOpt(float lopt[],
                          struct TimeCourse **timecoursestart)
{
  int i;
  struct TimeCourse *start;
  float dt1,dt2;
  
  for (i=0;i<NGenes;i++){
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

void CalcFitness(float lopt[],
                 float *w,
                 struct TimeCourse **timecoursestart,
                 float s)
{
  float d,dt1,dt2,x;
  int i;
  struct TimeCourse *start[NGenes];
  
  for (i=0;i<NGenes;i++) start[i]=timecoursestart[i];
  dt1=0.0;
  *w=0.0;
  while (start[0]){
    d = 0.0;
    if (start[0]->next) dt2 = (start[0]->next->time - start[0]->time)/2.0;
    else dt2=0.0;
    for (i=0;i<NGenes;i++){
      x = (start[i]->concentration - lopt[i]) / lopt[i];
      d += x*x;
    }
    d = sqrt(d);
    *w += exp(-s*d) * (dt1+dt2);     
    if (verbose) fprintf(fperrors,"t=%g dt=%g+%g d=%g w=%g\n",start[0]->time,dt1,dt2,d,*w/tdevelopment);
    dt1=dt2;
    for (i=0;i<NGenes;i++) start[i] = start[i]->next;
  }
  *w /= tdevelopment;
  fprintf(fperrors,"s=%g w=%g\n",s,*w);
}

void PrintTimeCourse(struct TimeCourse *start,
                     int i,
                     float lopt[])
{
  FILE *fpout;
  char filename[80];
  
  sprintf(filename,"output/protein%d.dat",i);
  if ((fpout = fopen(filename,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n",filename);
  while (start){
    fprintf(fpout,"%g %g\n",start->time,start->concentration);
    start = start->next;
  }
  //  fprintf(fpout, "%g %g\n", tdevelopment, lopt[i]);
  fclose(fpout);  
}

int main(int argc, char *argv[])
{
  FILE *fpout, *fpkdis;
  int i, j, k, gen;
  struct CellState state;
  struct Genotype indivs[PopSize];
  struct TimeCourse *timecoursestart[NGenes]; // array of pointers to list starts
  struct TimeCourse *timecourselast[NGenes];
  struct TimeCourse *start;
  float fitness[PopSize], sumfit, lopt[NGenes], initmRNA[NGenes], initProteinConc[NGenes], x, kdis[Nkdis];
  
  fperrors = fopen("netsimerrors.txt", "w");
  sumfit = 0.0;
  for (j=0; j<dummyrun; j++) ran1(&seed);
  for (i=0; i<NGenes; i++) {
    lopt[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initProteinConc[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }
  fpkdis = fopen("kdis.txt","r");
  for (j = 0; j < Nkdis; j++) {
    fscanf(fpkdis,"%f",&kdis[j]);
  }
  for (j = 0; j < PopSize; j++) {
    if (j==PopSize-1) output=1;
    InitializeGenotype(&indivs[j], kdis);
    InitializeCell(&state, indivs[j].mRNAdecay, initmRNA, initProteinConc);
    /* 
     *  AKL 2008-03-21: removed indivs[j].y: wasn't being used
     *  InitializeCell(&state,indivs[j].y,indivs[j].mRNAdecay,initmRNA,initProteinConc); 
     */
    Develop(&indivs[j], &state, (float) 293.0, timecoursestart, timecourselast);
    //    DevStabilityOnlyLOpt(lopt, timecoursestart);
    fprintf(fperrors,"indiv %d\n",j);
    /*    CalcFitness(lopt, &(fitness[j]), timecoursestart, selection);
          sumfit += fitness[j];*/
    for (i=0; i < NGenes; i++) {
      if ((output) && j==PopSize-1) PrintTimeCourse(timecoursestart[i], i, lopt);
      if (verbose) fprintf(fperrors, "deleting gene %d\n", i);
      DeleteTimeCourse(timecoursestart[i]);
      timecoursestart[i] = timecourselast[i] = NULL;
    }
    FreeMemCellState(&state);
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
    Mutate(&(indivs[j]),&(indivs[k]),mN/(float)PopSize);
    fprintf(fperrors,"finished mutating\n");fflush(fperrors);
    InitializeCell(&state,indivs[k].y,indivs[k].mRNAdecay,initmRNA,initProteinConc);
    Develop(&indivs[k],&state,(float) 293.0,timecoursestart,timecourselast);
    fprintf(fperrors,"finished development\n");fflush(fperrors);
    //    DevStabilityOnlyLOpt(lopt,timecoursestart);
    sumfit -= fitness[k];
    CalcFitness(lopt,&(fitness[k]),timecoursestart,selection);
    fflush(fperrors);
    sumfit += fitness[k];
    for (i=0; i<NGenes; i++){
      if (output && gen==Generations*PopSize) PrintTimeCourse(timecoursestart[i],i,lopt);
      if (verbose) fprintf(fperrors,"deleting gene %d\n",i);
      DeleteTimeCourse(timecoursestart[i]);
      timecoursestart[i] = timecourselast[i] = NULL;
    }
    FreeMemCellState(&state);   
  }
  */
  fclose(fperrors);
}
