/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-  */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>

#ifndef POP_SIZE
#define POP_SIZE 500         
#endif

#ifndef MAX_MUT_STEP         
#define MAX_MUT_STEP 10   // default 
#endif


#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6
#define TIME_INFINITY 9.99e10
#define RATE_OPERATIONS 1e6   /* number of operations on a rate before recomputing that rate */

#ifndef MAX_COPIES
#define MAX_COPIES 1       /* each gene can exists with four copies during replication */
#endif

//because of gene deletion and duplication, the number of genes are no longer constant
//now declared as global integers: line357-361
#ifdef  NO_SEPARATE_GENE    /* set to true if selection gene is also a TF */
  #ifndef TFGENES
  #define TFGENES 10          /* number of genes encoding TFs */
  #endif
  #ifndef NGENES
  #define NGENES TFGENES      /* total number of cisreg genes */
  #endif
  #ifndef NPROTEINS
  #define NPROTEINS TFGENES   /* total number of types of proteins, may be >#GENES if extracellular signals */
  #endif
  #define SELECTION_GENE_A (TFGENES-2)    /* index of selection gene */
  #define SELECTION_GENE_B (TFGENES-1)
#else                       /* otherwise, by default assuming selection gene is not a TF */
  #ifndef TFGENES             /* number of genes encoding TFs */
  #define TFGENES 3
  #endif
  #ifndef NGENES
  #define NGENES (TFGENES+2)  /* total number of genes: add the extra (non-TF) selection gene to the total (default case) */
  #endif
  #ifndef NPROTEINS           /* number of proteins: TODO: must be equal to the number of genes currently */
  #define NPROTEINS (TFGENES+2)
  #endif
  #define SELECTION_GENE_A (NGENES-2)    /* index of selection gene */
  #define SELECTION_GENE_B (NGENES-1)   /* index of selection gene: always the last two genes */
#endif

#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 6      /* length of binding element on TF */
#define NUM_K_DISASSEMBLY 133 /* number of differents for PIC disassembly from data file  */

#ifndef HIND_LENGTH
#define HIND_LENGTH 15         /* default length of hindrance (original was 6) */
#endif

//for parallelize mutation trials
#define N_para_threads 2

//int TFGENES = 3;
//int NGENES = 5;
//int NPROTEINS = 5;
//int SELECTION_GENE_A = 3;
//int SELECTION_GENE_B = 4;

/* 
 * define macros for logging output and warning/errors 
 */
#ifndef LOGGING_OFF
  #define LOG(...) { fprintf(fperrors, "[%s: cell %03d] ", __func__, state->cell_id); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); }
  #define LOG_NOCELLID(...) { fprintf(fperrors, "[%s] ", __func__); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); }
  #define LOG_ERROR(...) { fprintf(fperrors, "[%s: cell %03d ERROR] ", __func__, state->cell_id); \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_ERROR_NOCELLID(...) { fprintf(fperrors, "[%s ERROR] ", __func__); fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_WARNING(...) { fprintf(fperrors, "[%s: cell %03d WARNING] ", __func__, state->cell_id); \
    fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_WARNING_NOCELLID(...) { fprintf(fperrors, "[%s WARNING] ", __func__); fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_NOFUNC(...) { fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_VERBOSE_NOFUNC(...) if (verbose) { fprintf (fperrors, __VA_ARGS__); fflush(fperrors); }
  #define LOG_VERBOSE(...) if (verbose) { \
    if (state!=NULL)                  \
      fprintf(fperrors, "[%s: cell %03d] ", __func__, state->cell_id); \
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

/*
 * primary data structures for model
 */

/* 
 * enumeration for 'konvalues' indices
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
 * enum for state_change_ids
 */
enum { ACETYLATION_STATE = 0, 
       DEACETYLATION_STATE = 1, 
       PICASSEMBLY_STATE = 2,
       TRANSCRIPTINIT_STATE = 3, 
       PICDISASSEMBLY_STATE = 4,
};

/*
 * Rates for Gillespie algorithm
 *
 * Events with exponentially-distributed waiting times are:
 * - TF association and dissociation, 
 * - transcription initiation
 * - mRNA transport and decay
 * - acetylation and deacetylation, 
 * - and PIC assembly and disassembly.
 *
 */
typedef struct GillespieRates GillespieRates;
struct GillespieRates {
  float koff;              /* rates[0] */
  int koff_operations;     /* keeps track of when rounding errors should be recalc*/
  float transport;         /* rates[1] */
  int transport_operations;
  float mRNAdecay;         /* rates[2] */
  int mRNAdecay_operations;
  float pic_disassembly;    /* rates[3] */
  int pic_disassembly_operations;
  float salphc;            /* rates[4] */
  int salphc_operations;

  /* the following are cached values to avoid recomputation not
     actually part of the Gillespie rates */
  float max_salphc;         /* rates[5] */
  int max_salphc_operations;
  float min_salphc;         /* rates[6] */
  int min_salphc_operations;

  /* number of genes subject ot the following transitions, the rates of the
     relevant process are computed from these numbers */
  int acetylation_num[MAX_COPIES];       
  int deacetylation_num[MAX_COPIES];     
  int pic_assembly_num[MAX_COPIES];       
  int transcript_init_num[MAX_COPIES];    
  int pic_disassembly_num[MAX_COPIES];    

  /* subtotal, including above, but not including konrate */
  float subtotal;
};

/*
 * KonList: stores information about available binding sites for a
 * particular TF
 */
typedef struct KonList KonList;
struct KonList {
  int *available_sites;   /* list of available sites for this TF */ 
  int site_count;         /* number of available binding sites for a
                             given TF. Zero if protein is not a TF. */
};

/*
 * KonStates : composite data structure to cache information about
 * available binding sites to avoid re-computation.  This groups 
 * info for all TFs:
 *
 * konIDs, konvalues, nkonsum
 */
typedef struct KonStates KonStates;
struct KonStates {
  /* total number of currently *available* binding sites */
  int nkon;

  /* list of structs: need one for each protein */
  // TODO: currently KonList has cached information for both TF proteins and non-TF
  // proteins, may want to split this out at some point
  KonList *kon_list[NPROTEINS];

  /* konvalues are rates of binding with:
   * element 0 is (protein - salphc)/c
   * element 1 is c
   * element 2 is salphc
   * second index is which TF binds 
   * The kon term is left out of all of them for computational efficiency
   */
  float konvalues[NPROTEINS][3];
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         /* transcription factor */
  float Koff;        /* replacing hamming_dist */
  int BS_pos;        /* start position of BS on DNA, always with reference to forward strand */                     
  int N_hindered;      /* the number of BS hindered by this TF when it binds to the current BS */  
  //  int cisreg_id;     /* cis-reg region */
  //  int hamming_dist;  /* hamming distance */
  //  int strand;        /* strand 0 (forward) or 1 (backward)*/
  //  int gene_copy;     /* which copy of gene, 0 to MAX_COPIES-1 */
};

typedef struct Genotype Genotype;
struct Genotype {
  char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN];
  char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN];
  char tf_seq_rc[TFGENES][MAX_COPIES][TF_ELEMENT_LEN];  /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
  int binding_sites_num[NGENES];                        /* total number of binding sites */
  int max_hindered_sites[NGENES];                       /* maximal number of BS a BS can hinder*/ 
  int N_act_BS[NGENES];                                 /* total number of binding sites of activating TF */
  int N_rep_BS[NGENES];                                 /* total number of binding sites of repressing TF */
  int N_act;                                            /* number of activating TF*/
  int N_rep;                                            /* number of repressing TF*/
  AllTFBindingSites *all_binding_sites[NGENES];
  int activating[TFGENES][MAX_COPIES];                  /* 1 is activating TF, 0 is repressing */ 
  
  float mRNAdecay[NGENES];                              /* kinetic rates*/
  float proteindecay[NGENES];                           /* kinetic rates*/
  float translation[NGENES];                            /* kinetic rates*/   
  float pic_disassembly[NGENES][MAX_COPIES];  
  
//  int copies[NGENES];                 /* current per-gene ploidy */
  /* cached quantities for efficiency, can be recomputed from the above genotype, not part of model */
//  int sites_per_gene[NGENES];               /* cache number of TFBSs per gene */
//  int site_id_pos[NGENES][MAX_COPIES][2];  /* cache start and end positions of binding siteIDs for each cis-reg sequence for use during replication and division*/
};

/* 
 * transcription/translation delays are sorted linked lists.  Deleting
 * the head each time, and tack new stuff on the end.  Linked lists
 * are easy to create pre-sorted.
 */
typedef struct FixedEvent FixedEvent;
struct FixedEvent {
  int gene_id;
  int copy;
  float time;
  FixedEvent *next;
};

typedef struct CellState CellState;
struct CellState {
  int cell_id;                        /* cell ID */
  int burn_in;                        /* keep track of whether burn-in has finished */
  float cell_size;                    /* size of cell */
  float growth_rate;                  /* total growth rate in the previous deltat */
  int mRNA_cyto_num[NGENES];          /* mRNAs in cytoplasm */
  int mRNA_nuclear_num[NGENES];       /* mRNAs in nucleus */
  int mRNA_transl_cyto_num[NGENES];   /* mRNAs are in the cytoplasm, but only recently */
  FixedEvent *mRNA_transl_time_end;   /* times when mRNAs become fully loaded with ribosomes and start producing protein */
  FixedEvent *mRNA_transl_time_end_last;
  int mRNA_transcr_num[NGENES][MAX_COPIES];  /* mRNAs which haven't finished transcription yet */
  FixedEvent *mRNA_transcr_time_end;  /* times when transcription is complete and an mRNA is available to move to cytoplasm */
  FixedEvent *mRNA_transcr_time_end_last;
  FixedEvent *env0_time_end; // times when env=0 ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.
  FixedEvent *env0_time_end_last;
  FixedEvent *env1_time_end;
  FixedEvent *env1_time_end_last;
   

  float protein_conc[NPROTEINS];
  int tf_bound_num;
  int *tf_bound_indexes;                /* indices in all_binding_sites that are currently bound*/
  int tf_hindered_num;
  int (*tf_hindered_indexes)[2];
  /* 1st elem tf_hindered_indexes lists binding site indices that cannot be bound due to steric hindrance
   * 2nd elem gives corresponding index of inhibiting TF in all_binding_sites, so that we know when to release hindrance
   * binding sites can be hindered more than once, then multiple constraints must be lifted before TF binding
   */
  int active[NGENES][MAX_COPIES];
  /* gives the state of each of the genes, according to figure
   *  1 is fully off, 2 meets TF criteria
   *  3 is off but w/o nucleosome, 4 is on but w/o PIC
   *  5 is on but w/o TF criteria, 6 is fully on
   * see the enum definition above
   */

  /* stores corresponding gene_ids ready for [de]acteylation,
     PIC[dis]assembly, transcriptinit, see enum above*/
  int state_change_ids[5][MAX_COPIES][NGENES]; 
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

/*
 * global variables
 */

/* see netsim.c for documentation for these global constant variables */
extern const int MAXELEMENTS; 
const int MAXBOUND;
const int NMIN;
const float KRNA;
const float TTRANSLATION;
const float TTRANSCRIPTION;
const float PROB_ACTIVATING;
const float TRANSCRIPTINIT; 
const float DEACETYLATE;
const float ACETYLATE;
const float PICASSEMBLY;
const float STARTNUCLEUS;
const float KR;    
const float GASCONSTANT;
const float COOPERATIVITY;
const float COOPERATIVE_DISTANCE; 
const float NUMSITESINGENOME ;
const float mN ;     



/* see netsim.c for documentation for these global variables */
float kon; 
float kon_after_burnin; 
int burn_in;
float tdevelopment;  
float timemax;       
int current_ploidy;  
int output;
long seed ;          
int dummyrun;        
int recompute_koff;  
int recompute_kon;
float growth_rate_scaling; 
float Pp_a;
float Pp_b;         
float h;
float gmax_a;
float gmax_b;
float protein_aging;



/* file output parameters */
char *output_directory ;
int verbose ;
FILE *fperrors;
FILE *fp_cellsize[2];
#if 0 
FILE *fp_koff[2];
FILE *fp_growthrate[2];
FILE *fp_tfsbound[2];
FILE *fp_rounding[2];
#endif

/* function prototypes */

extern void initialize_parameters();

extern void initialize_growth_rate_parameters();

extern void initialize_sequence(char [], int, int, int);

extern void print_genotype(Genotype *, int);

extern void print_all_binding_sites(int [NGENES],
                                    AllTFBindingSites *, 
                                    int ,
                                    char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                                    char [NGENES][MAX_COPIES][CISREG_LEN]
//                                    int [NGENES][MAX_COPIES][2]
									);

extern void print_tf_occupancy(CellState *,
                               AllTFBindingSites *,
                               float);

extern void initialize_genotype(Genotype *, 
                                Genotype *,
                                float [],
                                int);

extern void initialize_genotype_fixed(Genotype *,
                                      float *,
                                      int);

extern void calc_all_binding_sites_copy(Genotype *, int);

extern void calc_all_binding_sites(Genotype *);

extern float calc_ratio_act_to_rep(AllTFBindingSites *, 
                                    int,
                                    int,
                                    int,
                                    int, 
                                    int [NGENES][MAX_COPIES], 
                                    float *);

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
                            float [NPROTEINS],
                            int);

extern void initialize_cell_cache(CellState *,
                                  Genotype,
                                  KonStates *,
                                  float **,
                                  int,
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
                                FixedEvent *,
                                float);

extern void calc_dt(float *,
                    float *,
                    GillespieRates *,
                    KonStates *,
                    float [],
                    float [],
                    int [],
                    int [],
                    int);

extern void end_transcription(float *,
                              float,
                              CellState *,
                              float [NGENES],
                              GillespieRates *);


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
                            float [NPROTEINS],
                            TimeCourse **,
                            TimeCourse **);

extern void add_integer_time_points(float,
                                    int [NPROTEINS],
                                    TimeCourse **,
                                    TimeCourse **);

extern float compute_tprime(float, float, float, float);

extern float compute_integral(float, float, float, float, float, float, float, float);

extern float compute_growth_rate_dimer(float *,
                                       float , 
                                       float ,
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
                                       float,
                                       float,
                                       float,
				       int );

extern void update_protein_conc_cell_size(float[],
                                          CellState *,
                                          Genotype *,
                                          float,
                                          GillespieRates *,
                                          KonStates *,
                                          float,
//                                          TimeCourse **,
//                                          TimeCourse **,
                                          float [],
					  int *);

extern void calc_num_bound(float[],
                           int );

extern int sum_rate_counts(int[MAX_COPIES]);

extern void get_gene(int [MAX_COPIES], int, int *, int *);

extern void transport_event(GillespieRates *,
                            CellState *,
                            Genotype *,
                            KonStates *,
                            float [NGENES],
//                            TimeCourse **, 
//                            TimeCourse **, 
                            float,
                            float,
                            float,
                            int *);

extern void tf_binding_event(GillespieRates *, CellState *, Genotype *, 
                             KonStates *, float *, 
//							 TimeCourse **, TimeCourse **,
                             float, float, float, int, int, int,int *);

extern void tf_unbinding_event(GillespieRates *, CellState *, Genotype *, 
                               KonStates *, float *, 
//							   TimeCourse **, TimeCourse **,
                               float, float, float, float, int, int *,int *);

extern void mRNA_decay_event(GillespieRates *, CellState *, Genotype *, 
                             KonStates *, float *,
//							  TimeCourse **, TimeCourse **,
                             float, float, float, int *);

extern void histone_acteylation_event(GillespieRates *, CellState *, Genotype *, 
                                      KonStates *, 
//									  TimeCourse **, TimeCourse **,
                                      float, float, int *);

extern void histone_deacteylation_event(GillespieRates *, CellState *, Genotype *, 
                                        KonStates *,
//										 TimeCourse **, TimeCourse **,
                                        float, float, int *);

extern void assemble_PIC_event(GillespieRates *, CellState *, Genotype *, 
                               KonStates *, 
//							   TimeCourse **, TimeCourse **,
                               float, float, int *);

extern void disassemble_PIC_event(GillespieRates *, CellState *, Genotype *, 
                                  KonStates *, 
//								  TimeCourse **, TimeCourse **,
                                  float, float, float);

extern void transcription_init_event(GillespieRates *, CellState *, Genotype *,
                                     KonStates *,
//									  TimeCourse **, TimeCourse **,
                                     float, float, float, int *);

extern void shift_binding_site_ids(CellState *, 
                                   KonStates *,
                                   int,
                                   int);


extern void recompute_kon_rates(GillespieRates *,
                                CellState *,
                                Genotype *,
                                KonStates *,
                                int);

extern void recompute_koff_rates(GillespieRates *,
                                 CellState *,
                                 Genotype *,
                                 float *,
                                 float);

extern void recalibrate_cell(GillespieRates *,
                             CellState *,
                             Genotype *,
                             KonStates *,
                             float **,
                             float [NGENES],
                             float [NGENES],
                             float);

extern int do_single_timestep(Genotype *, 
                               CellState *, 
                               KonStates *, 
                               GillespieRates *, 
                               float *,
                               float *,
                               float *,
                               float *,
                               float *,
                               float *,
                               float *,
//                               TimeCourse *[NPROTEINS],
//                               TimeCourse *[NPROTEINS],
                               int,
                               int,
                               int,
                               int,
                               int *) ;
							   
extern void free_fixedevent(CellState *);
 
extern float calc_avg_growth_rate(int,
                                  Genotype *, 
                                    CellState *, 
                                    float [NGENES],
                                    float [NGENES],
                                    float *,
                                    float *,
                                    float [NGENES],
                                    float [NGENES],
                                    float *,
                                    float *,
                                    float *,
                                    KonStates *,
                                    GillespieRates *,
                                    float ,
                                    float ,
                                    float ,
                                    int ); 
                                    
extern int try_fixation(float, float);

extern int mutate(Genotype *);
  
extern void init_run_pop(Genotype [N_para_threads],
                         CellState [N_para_threads],
//                         TimeCourse *[2][NGENES],
//                         TimeCourse *[2][NGENES], 
                         float, /* in Kelvin */
                         float [NUM_K_DISASSEMBLY],
                         int,
                         int);


extern void print_time_course(TimeCourse *,
                              int, int);

extern void print_all_protein_time_courses(TimeCourse *[2][NPROTEINS],
                                          TimeCourse *[2][NPROTEINS]);
                                          
extern void clone_cell(Genotype *,                
                	Genotype *,
			int);

extern void log_snapshot(GillespieRates *,
                         CellState *,
                         Genotype *,
                         KonStates *,
                         float **,
                         float [NGENES],
                         float [NGENES],
                         float ,
                         float ,
                         float );

extern int mod(int, int); //thanks to the stupid % in gcc

#endif /* !FILE_NETSIM_SEEN */
