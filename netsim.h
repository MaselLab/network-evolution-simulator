/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- 
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>

#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */

#ifndef PLOIDY
#define PLOIDY 2           /* 1 = haploid, 2 = diploid */
#endif
#define NGENES 10          /* number of genes */
#define CISREG_LEN 150     /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 6   /* length of binding element on TF */
#define NUM_K_DISASSEMBLY 133 /* number of differents for PIC disassembly from data file  */
#define HIND_LENGTH 6         /*length of hindrance*/
//#define HIND_OLD 6

extern int verbose;
extern FILE *fperrors;

/* 
 * enum for 'konvalues' indices
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
 * enum for stateChangeIDs
 */
enum { ACTEYLATION = 0, 
       DEACTEYLATION = 1, 
       PICASSEMBLY = 2,
       TRANSCRIPTINIT = 3, 
       PICDISASSEMBLY = 4,
};

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
typedef struct GillespieRates GillespieRates;
struct GillespieRates {
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
  int acetylationCount[PLOIDY];       /* rates2[0] */
  int deacetylationCount[PLOIDY];     /* rates2[1] */
  int picAssemblyCount[PLOIDY];       /* rates2[2] */
  int transcriptInitCount[PLOIDY];    /* rates2[3] */
  int picDisassemblyCount[PLOIDY];    /* rates2[4] */

  /* total, including rates2 */
  float total;                /* rates[7] */
};

/*
 * enum for konIDs
 */
enum { SITEID_INDEX = 0, TFID_INDEX = 1 };

/*
 * enum for KonList
 */
/* enum { SITEID_INDEX = 0, TOTAL_AVAILABLE = 1 }; */
typedef struct KonList KonList;
struct KonList {
  int *available_sites;
  int site_count;
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
typedef struct KonStates KonStates;
struct KonStates {
  /* number of currently *available* binding sites */
  int nkon;

  /* elem 0 is siteID, elem 1 is TF that binds */
  int (*konIDs)[2];

  /* list of structs  */
  KonList *konList[NGENES];

  /* number of available binding sites for a given TF */
  int nkonsum[NGENES];

  /* konvalues are rates of binding with:
   * element 0 is other funny term (TODO ?)
   * element 1 is c
   * element 2 is salphc
   * Other index is which TF binds 
   * The kon term is left out of all of them for computational efficiency
   */
  float konvalues[NGENES][3];
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int cisregID;     /* 0 - cis-reg region*/
  int tfID;         /* 1 - trnascription factor*/
  int sitePos;      /* 2 - start position of binding site*/
  int strand;       /* 3 - strand*/
  int hammingDist;  /* 4 - hamming distance*/
  int cisregCopy;    /* which copy of gene */
  int tfCopy;        /* which copy of gene */ 
  int hindPos;       /*position of BS with in 15bp hindrance (offset)*/
  int leftEdgePos;   /*start position of 15bp hindrance*/
};

typedef struct Genotype Genotype;
struct Genotype {
  char cisRegSeq[NGENES][PLOIDY][CISREG_LEN];
  char transcriptionFactorSeq[NGENES][PLOIDY][TF_ELEMENT_LEN];
  int hindrancePositions[NGENES];     /*offset positions of BS based on TF--heritable*/
  int bindSiteCount;
  AllTFBindingSites *allBindingSites;
/* int (*allBindingSites)[5];
 5 elements are
    identity of cis-regulatory region
    identity of TF that binds
    position 0 to 386 (first position, always in forwards strand)  
    strand 0 (forward) or 1 (backward)
    Hamming distance 0 to n-n_min 
*/ 
  float mRNAdecay[NGENES];
  float proteindecay[NGENES];
  float translation[NGENES];
  int activating[NGENES][PLOIDY]; /* 1 is activating, 0 is repressing */
  float PICdisassembly[NGENES][PLOIDY];
};

/* 
 * transcription/translation delays are sorted linked lists.  Deleting
 * the head each time, and tack new stuff on the end.  Linked lists
 * are easy to create pre-sorted.  */
typedef struct FixedEvent FixedEvent;
struct FixedEvent {
  int geneID;
  float time;
  FixedEvent *next;
};

typedef struct CellState CellState;
struct CellState {
  float cellSize;                     /* size of cell */
  int mRNACytoCount[NGENES];          /* mRNAs in cytoplasm */
  int mRNANuclearCount[NGENES];       /* mRNAs in nucleus */
  int mRNATranslCytoCount[NGENES];    /* mRNAs are in the cytoplasm, but only recently */
  FixedEvent *mRNATranslTimeEnd;    /* times when mRNAs move to cytoplasm (mRNACytoCount) */
  FixedEvent *mRNATranslTimeEndLast; 
  int mRNATranscrCount[NGENES];             /* mRNAs which haven't finished transcription yet */
  FixedEvent *mRNATranscrTimeEnd;    /* times when transcripts ? move to nucleus (mRNANuclearCount) */
  FixedEvent *mRNATranscrTimeEndLast;
  float proteinConc[NGENES];
  int tfBoundCount;
  int *tfBoundIndexes;
  /* tfBoundIndexes lists just bound TFs according to binding site index in G */
  int tfHinderedCount;
  int (*tfHinderedIndexes)[2];
  /*1st elem tfHinderedIndexes lists binding site indices that cannot be bound due to steric hindrance
    2nd elem gives corresponding index of inhibiting TF in G (TODO: what is this?)
  */
  int active[NGENES][PLOIDY];
  /* gives the state of each of the genes, according to figure
     1 is fully off, 2 meets TF criteria
     3 is off but w/o nucleosome, 4 is on but w/o PIC
     5 is on but w/o TF criteria, 6 is fully on
  */

  /* stores corresponding geneIDs for [de]acteylation, PIC[dis]assembly, transcriptinit */
  int statechangeIDs[5][PLOIDY][NGENES]; 
  float RTlnKr;
  float temperature;
};

typedef struct TimeCourse TimeCourse;

struct TimeCourse
{
  float concentration;
  float time;
  TimeCourse *next;
};     

extern void calc_interaction_matrix(char [NGENES][PLOIDY][CISREG_LEN],
                                    char [NGENES][PLOIDY][TF_ELEMENT_LEN],
                                    int *,
                                    AllTFBindingSites **,
                                    int [NGENES]);



#endif /* !FILE_NETSIM_SEEN */
