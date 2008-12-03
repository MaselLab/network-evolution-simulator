/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- 
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>

#ifndef POP_SIZE
#define POP_SIZE 1
#endif

#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6

#ifndef MAX_COPIES
#define MAX_COPIES 4       /* each gene can potentially exist as a tetraploid during replication */
#endif

#ifdef  NO_SEPARATE_GENE
#define TFGENES 10          /* number of genes encoding TFs */
#define NGENES TFGENES      /* total number of genes */
#define SELECTION_GENE 9    /* index of selection gene */
#else
#define TFGENES 10          /* number of genes encoding TFs */
#define NGENES (TFGENES+1)  /* total number of genes: add the extra (non-TF) selection gene to the total (default case) */
#define SELECTION_GENE TFGENES   /* index of selection gene */
#endif

#define CISREG_LEN 150     /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 6   /* length of binding element on TF */
#define NUM_K_DISASSEMBLY 133 /* number of differents for PIC disassembly from data file  */

#ifndef HIND_LENGTH
#define HIND_LENGTH 15         /* default length of hindrance (original was 6) */
#endif

#ifndef LOGGING_OFF
#define LOG(...) { fprintf(fperrors, "[%s: cell %03d] ", __func__, state->cellID); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); }
#define LOG_NOCELLID(...) { fprintf(fperrors, "[%s] ", __func__); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); }
#define LOG_ERROR(...) { fprintf(fperrors, "[%s: cell %03d ERROR] ", __func__, state->cellID); \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#define LOG_ERROR_NOCELLID(...) { fprintf(fperrors, "[%s ERROR] ", __func__); fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#define LOG_WARNING(...) { fprintf(fperrors, "[%s: cell %03d WARNING] ", __func__, state->cellID); \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#define LOG_WARNING_NOCELLID(...) { fprintf(fperrors, "[%s WARNING] ", __func__); fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#define LOG_NOFUNC(...) { fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#define LOG_VERBOSE_NOFUNC(...) if (verbose) { fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#define LOG_VERBOSE(...) if (verbose) { \
    if (state!=NULL)                  \
      fprintf(fperrors, "[%s: cell %03d] ", __func__, state->cellID); \
    else \
      fprintf(fperrors, "[%s] ", __func__);  \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#define LOG_VERBOSE_NOCELLID(...) if (verbose) { \
    fprintf(fperrors, "[%s] ", __func__);                 \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
#else
#define LOG
#define LOG_ERROR
#define LOG_NOFUNC
#define LOG_VERBOSE
#endif

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
enum { ACETYLATION = 0, 
       DEACETYLATION = 1, 
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
  int acetylationCount[MAX_COPIES];       
  int deacetylationCount[MAX_COPIES];     
  int picAssemblyCount[MAX_COPIES];       
  int transcriptInitCount[MAX_COPIES];    
  int picDisassemblyCount[MAX_COPIES];    

  /* total, including above */
  float total;
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
 * available binding sites to avoid re-computation.  This groups the
 * previous data structures:
 *
 * konIDs, konvalues, nkonsum
 */
typedef struct KonStates KonStates;
struct KonStates {
  /* number of currently *available* binding sites */
  int nkon;

  /* list of structs  */
  KonList *konList[NGENES];

  // check
  // TODO: remove, this is redundant with konList above
  //
  /* number of available binding sites for a given TF */
  int nkonsum[NGENES];

  /* konvalues are rates of binding with:
   * element 0 is (protein - salphc)/c
   * element 1 is c
   * element 2 is salphc
   * second index is which TF binds 
   * The kon term is left out of all of them for computational efficiency
   */
  float konvalues[NGENES][3];
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int cisregID;     /* cis-reg region */
  int tfID;         /* transcription factor */
  int strand;       /* strand 0 (forward) or 1 (backward)*/
  int hammingDist;  /* hamming distance */
  int geneCopy;     /* which copy of gene, 0 to MAX_COPIES-1 */
  int leftEdgePos;  /* start position of TF on DNA, always with reference to forward strand */
                    /* note that recognition sitePos, can be computed from: (leftEdgePos+hindPos) */
  int hindPos;      /* position of recognition site within the HIND_LENGTH bp hindrance (offset) */
};

typedef struct Genotype Genotype;
struct Genotype {
  char cisRegSeq[NGENES][MAX_COPIES][CISREG_LEN];
  char transcriptionFactorSeq[NGENES][MAX_COPIES][TF_ELEMENT_LEN];
  int hindrancePositions[NGENES];     /* offset positions of each TF's hindrance area relative to recognition site*/
  int bindSiteCount;
  AllTFBindingSites *allBindingSites;
  int tfsPerGene[NGENES];               /* cache number of TFs per gene */
  int tfsStart[NGENES][MAX_COPIES][2];  /* cache start and end positions of binding siteIDs */
  float mRNAdecay[NGENES];
  float proteindecay[NGENES];
  float translation[NGENES];
  int activating[NGENES][MAX_COPIES]; /* 1 is activating, 0 is repressing */
  float PICdisassembly[NGENES][MAX_COPIES];
  int copies[NGENES];                 /* current per-gene ploidy */
  float replication_time[NGENES];     /* per-gene replication time */
};

/* 
 * transcription/translation delays are sorted linked lists.  Deleting
 * the head each time, and tack new stuff on the end.  Linked lists
 * are easy to create pre-sorted.  */
typedef struct FixedEvent FixedEvent;
struct FixedEvent {
  int geneID;
  int copy;
  float time;
  FixedEvent *next;
};

typedef struct CellState CellState;
struct CellState {
  int cellID;                         /* cell ID */
  int founderID;                      /* keep track of the founder cell */
  int burn_in;                        /* keep track of whether burn-in has finished */
  int in_s_phase;                     /* whether cell has entered S (synthesis) phase */
  int divisions;                      /* total number of divisions cell has undergone */
  float division_time;                /* current division time */ 
  float cellSize;                     /* size of cell */
  float growthRate;                   /* total growth rate in the previous deltat */
  int mRNACytoCount[NGENES];          /* mRNAs in cytoplasm */
  int mRNANuclearCount[NGENES];       /* mRNAs in nucleus */
  int mRNATranslCytoCount[NGENES];    /* mRNAs are in the cytoplasm, but only recently */
  FixedEvent *mRNATranslTimeEnd;      /* times when mRNAs become fully loaded with ribosomes and start producing protein */
  FixedEvent *mRNATranslTimeEndLast; 
  int mRNATranscrCount[NGENES][MAX_COPIES];       /* mRNAs which haven't finished transcription yet */
  FixedEvent *mRNATranscrTimeEnd;     /* times when transcription is complete and an mRNA is available to move to cytoplasm*/
  FixedEvent *mRNATranscrTimeEndLast;

  FixedEvent *replicationTimeEnd;     /* times when gene duplicates */
  FixedEvent *replicationTimeEndLast;    

  float proteinConc[NGENES];
  int tfBoundCount;
  int *tfBoundIndexes;                /* tfBoundIndexes lists just bound TFs according to binding site index in all_binding_sites */
  int tfHinderedCount;
  int (*tfHinderedIndexes)[2];
  /* 1st elem tfHinderedIndexes lists binding site indices that cannot be bound due to steric hindrance
   * 2nd elem gives corresponding index of inhibiting TF in all_binding_sites, so that we know when to release hindrance
   * binding sites can be hindered more than once, then multiple constraints must be lifted before TF binding
   */
  int active[NGENES][MAX_COPIES];
  /* gives the state of each of the genes, according to figure
     1 is fully off, 2 meets TF criteria
     3 is off but w/o nucleosome, 4 is on but w/o PIC
     5 is on but w/o TF criteria, 6 is fully on
  */

  /* stores corresponding geneIDs for [de]acteylation, PIC[dis]assembly, transcriptinit */
  int statechangeIDs[5][MAX_COPIES][NGENES]; 
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

extern const int maxelements; 
/* start by allocating maxelements when initializing a genotype, double as needed, reduce at end */
const int maxbound;
const int PopSize;
const int nmin;
const float kRNA;
const float ttranslation;
const float ttranscription;
const float pact;
const float transcriptinit; /* replace betaon and betaoff */
const float deacetylate;
const float acetylate;
const float PICassembly;
const float startnucleus;
const float Kr;    /* don't put this less than 1, weird things happen to koff calculation */
const float GasConstant;
const float cooperativity;/* dGibbs, relative to 1 additional specific nt */
const float cooperative_distance;  /* distance co-operativity operates, changed from 20 */ 
const float NumSitesInGenome ; /* updated from 1.8e+6 */
const float selection ;

const float mN ;
const int Generations;

float kon; /* lower value is so things run faster */
/* kon=0.2225 is based on 1 molecule taking 240seconds=4 minutes
   and 89% of the proteins being in the nucleus*/
float kon_after_burnin; 
int burn_in;

float tdevelopment;  /* default development time: can be changed at runtime */
float timemax;       /* set an upper limit to development time */
int current_ploidy;    /* ploidy can be changed at run-time: 1 = haploid, 2 = diploid */
int output;
long seed ;         /* something is wrong here: changing seed changes nothing */
int dummyrun;       /* used to change seed */
int recompute_koff; /* toggle whether to recompute certain features at each time to avoid
                       compounding rounding error */
float critical_size ; /* critical size at which cell divides, 
                         set to negative to prevent division  */
float growth_rate_scaling; /* set growth rate scaling factor */

float time_s_phase;  /* length of S phase (default: 30 mins) */
float time_g2_phase; /* length of G2/mitosis phase (default: 30 mins) */
int random_replication_time; /* toggle replication times in S phase being random */

/* file output parameters */
char *output_directory ;
int verbose ;
FILE *fperrors;
FILE *fp_cellsize[POP_SIZE];
FILE *fp_koff[POP_SIZE];
FILE *fp_growthrate[POP_SIZE];
FILE *fp_tfsbound[POP_SIZE];

/* protein aging term: used when c=c'+g=0, set to 1e-4 < mean-3*sd of
   Belle et al. (2006) and small with respect to most growth rates  */
float protein_aging;

/* growth rate parameters globally used during simulation */
float Pp;         
float h;
float gmax;


extern void initialize_parameters();

extern void initialize_growth_rate_parameters();

extern void initialize_sequence(char [], int, int);

extern void print_all_binding_sites(int [NGENES],
                                    AllTFBindingSites *, 
                                    int ,
                                    char [NGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                    char [NGENES][MAX_COPIES][CISREG_LEN],
                                    int [NGENES][MAX_COPIES][2]);

extern void print_tf_occupancy(CellState *,
                               AllTFBindingSites *,
                               float);
extern void initialize_genotype(Genotype *, 
                                float []);


extern void mutate(Genotype *,
                   int,
                   int,
                   float);

extern int calc_all_binding_sites_copy(char [NGENES][MAX_COPIES][CISREG_LEN],
                                       char [NGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                       int ,
                                       AllTFBindingSites **,
                                       int *,
                                       int ,
                                       int ,
                                       int [NGENES]);

extern void calc_all_binding_sites(int [NGENES],
                                   char[NGENES][MAX_COPIES][CISREG_LEN],
                                   char[NGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                   int *,
                                   AllTFBindingSites **,
                                   int [NGENES],
                                   int [NGENES],
                                   int [NGENES][MAX_COPIES][2]);

extern int add_fixed_event(int,
                           int,
                           float,
                           FixedEvent **,
                           FixedEvent **);

extern void add_time_point(float,
                           float,
                           TimeCourse **,
                           TimeCourse **);

extern void add_fixed_event_end(int,
                                int,
                                float,
                                FixedEvent **,
                                FixedEvent **);

extern void delete_fixed_event(int,
                               int,
                               int,
                               FixedEvent **,
                               FixedEvent **);

extern void delete_fixed_event_start(FixedEvent **,
                                     FixedEvent **);

extern void initialize_cell(CellState *,
                            int,
                            int [NGENES],
                            float [NGENES],
                            float [NGENES],
                            float [NGENES],
                            int);

extern void calc_time (float, 
                       float, 
                       GillespieRates *,
                       KonStates *,
                       float *, 
                       float *);

extern void calc_kon_rate(float,
                          KonStates *,
                          float *);

extern void change_mRNA_cytoplasm(int,
                                  Genotype *,
                                  CellState *,
                                  GillespieRates *,
                                  KonStates *);

extern void calc_koff(int,
                      AllTFBindingSites *,
                      CellState *,
                      float *,
                      float);

extern void scan_nearby_sites(int,
                              AllTFBindingSites *,
                              CellState *,
                              GillespieRates *,
                              float *,
                              float);

extern void remove_kon(int,
                       int,
                       GillespieRates *,
                       float,
                       KonStates *,
                       float);

extern void add_kon(float,
                    float,
                    int,
                    int,
                    GillespieRates *,
                    KonStates *);

extern int ready_to_transcribe(int,
                               int,
                               int *,
                               int,
                               AllTFBindingSites *,
                               int [NGENES][MAX_COPIES],
                               int *);

extern int is_one_activator(int,
                            int,
                            int *,
                            int,
                            AllTFBindingSites *,
                            int [NGENES][MAX_COPIES]);

extern void calc_from_state(Genotype *,
                            CellState *,
                            GillespieRates *,
                            KonStates *,
                            float [],
                            float []);

extern int does_fixed_event_end(FixedEvent *,
                                FixedEvent *,
                                FixedEvent *,
                                float);

extern void calc_dt(float *,
                    float *,
                    GillespieRates *,
                    KonStates *,
                    float [],
                    float [],
                    int [],
                    int []);

extern void end_transcription(float *,
                              float,
                              CellState *,
                              float [NGENES],
                              GillespieRates *);

extern void remove_from_array(int,
                              int,
                              int [],
                              int *,
                              int );

extern void disassemble_PIC(CellState *,
                            Genotype *,
                            int,
                            int,
                            GillespieRates *);

extern void revise_activity_state(int,
                                  int,
                                  Genotype *,
                                  CellState *,
                                  GillespieRates *);

extern void remove_tf_binding(Genotype *,
                              CellState *,
                              GillespieRates *,
                              KonStates *,
                              int,
                              float [],
                              float);

extern void attempt_tf_binding(Genotype *,
                               CellState *,
                               GillespieRates *,
                               float **,
                               KonStates *,
                               int *,
                               int *,
                               int,
                               float);

extern void add_time_points(float,
                            float [NGENES],
                            TimeCourse **,
                            TimeCourse **);

extern void add_integer_time_points(float,
                                    int [NGENES],
                                    TimeCourse **,
                                    TimeCourse **);

extern void reach_s_phase(CellState *, Genotype *, float);

extern float compute_tprime(float, float, float, float);

extern float compute_integral(float, float, float, float, float, float, float, float, float);

extern float compute_growth_rate_dimer(float *,
                                       float , 
                                       float ,
                                       float [NGENES],
                                       int [NGENES],
                                       float,
                                       float,
                                       float, 
                                       float,
                                       float,
                                       float,
                                       float,
                                       int);

extern void update_protein_conc_cell_size(float[],
                                          CellState *,
                                          Genotype *,
                                          float,
                                          GillespieRates *,
                                          KonStates *,
                                          float,
                                          TimeCourse **,
                                          TimeCourse **,
                                          float []);

extern void calc_num_bound(float[],
                           int );

extern int sum_rate_counts(int[MAX_COPIES]);

extern void get_gene(int [MAX_COPIES], int, int *, int *);

extern void transport_event(GillespieRates *,
                            CellState *,
                            Genotype *,
                            KonStates *,
                            float [NGENES],
                            TimeCourse **, 
                            TimeCourse **, 
                            float,
                            float,
                            float,
                            float);

extern void tf_binding_event(GillespieRates *, CellState *, Genotype *, 
                             KonStates *, float *, TimeCourse **, TimeCourse **,
                             float, float, float, int, int);

extern void tf_unbinding_event(GillespieRates *, CellState *, Genotype *, 
                               KonStates *, float *, TimeCourse **, TimeCourse **,
                               float, float, float, float);

extern void mRNA_decay_event(GillespieRates *, CellState *, Genotype *, 
                             KonStates *, float *, TimeCourse **, TimeCourse **,
                             float, float, float);

extern void histone_acteylation_event(GillespieRates *, CellState *, Genotype *, 
                                      KonStates *, TimeCourse **, TimeCourse **,
                                      float, float);

extern void histone_deacteylation_event(GillespieRates *, CellState *, Genotype *, 
                                        KonStates *, TimeCourse **, TimeCourse **,
                                        float, float);

extern void assemble_PIC_event(GillespieRates *, CellState *, Genotype *, 
                               KonStates *, TimeCourse **, TimeCourse **,
                               float, float);

extern void disassemble_PIC_event(GillespieRates *, CellState *, Genotype *, 
                                  KonStates *, TimeCourse **, TimeCourse **,
                                  float, float, float);

extern void transcription_init_event(GillespieRates *, CellState *, Genotype *,
                                     KonStates *, TimeCourse **, TimeCourse **,
                                     float, float, float);

extern void shift_binding_site_ids(CellState *, 
                                   KonStates *,
                                   int,
                                   int);

extern void replicate_gene(CellState *,
                           Genotype *,
                           GillespieRates *,
                           KonStates *,
                           float *,
                           int,
                           float);

extern float do_single_timestep(Genotype *, 
                                CellState *, 
                                KonStates *, 
                                GillespieRates *, 
                                float *,
                                float *,
                                float [NGENES],
                                float [NGENES],
                                float *,
                                float *,
                                float *,
                                TimeCourse *[NGENES],
                                TimeCourse *[NGENES],
                                int,
                                int,
                                int) ;
  
extern void develop(Genotype [POP_SIZE],
                    CellState [POP_SIZE],
                    TimeCourse *[POP_SIZE][NGENES],
                    TimeCourse *[POP_SIZE][NGENES], 
                    float, /* in Kelvin */
                    float [NGENES],
                    float [NGENES],
                    float [NUM_K_DISASSEMBLY],
                    int,
                    int,
                    int,
                    int); 

extern void print_time_course(TimeCourse *,
                              int, int);

#endif /* !FILE_NETSIM_SEEN */
