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
#define POP_SIZE 10000        
#endif

#ifndef MAX_MUT_STEP         
#define MAX_MUT_STEP 10   // default 
#endif


#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6
#define TIME_INFINITY 9.99e10
//#define RATE_OPERATIONS 1e6   /* number of operations on a rate before recomputing that rate */

#ifndef MAX_COPIES
#define MAX_COPIES 1       /* each gene can exists with four copies during replication. This is an old setting*/
#endif

/* Because mutation can change the number of genes, the numbers defined here are used to allocate storage space only.
 * Set the numbers to be 8 folds of the initial ngenes and ntfgenes, so that we can have two whole genome duplications*/

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
  #define TFGENES 80          /* the initial value is set in initiate_genotype*/
  #endif
  #ifndef NGENES
  #define NGENES 96  /* total number of genes: add the extra (non-TF) selection gene to the total (default case) */
  #endif
  #ifndef NPROTEINS           /* number of proteins: TODO: must be equal to the number of genes currently */
  #define NPROTEINS 96
  #endif
//  #define SELECTION_GENE_A (NGENES-2)    /* index of selection gene */
//  #define SELECTION_GENE_B (NGENES-1)   /* index of selection gene: always the last two genes */
#endif

#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 6      /* length of binding element on TF */
#define NUM_K_DISASSEMBLY 133 /* number of differents for PIC disassembly from data file  */

#ifndef HIND_LENGTH
#define HIND_LENGTH 15         /* default length of hindrance (original was 6) */
#endif

//for parallelize mutation trials
#define N_para_threads 2

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
enum { KON_PROTEIN_DECAY_INDEX = 0, KON_SALPHC_INDEX = 1 };

/*
 * enum for 'CellState'->active indices
 */

enum { NUC_NO_PIC = 0,
       NO_NUC_NO_PIC =1,
       PIC_NO_NUC = 3,};
       
/*
 * enum for state_change_ids
 */
//enum { ACETYLATION_STATE = 0, 
//       DEACETYLATION_STATE = 1, 
//       PICASSEMBLY_STATE = 2,
//       TRANSCRIPTINIT_STATE = 3, 
//       PICDISASSEMBLY_STATE = 4,
//};

/*
 * Rates for Gillespie algorithm
 *
 * Events with exponentially-distributed waiting times are:
 * - transcription initiation
 * - mRNA transport and decay
 * - acetylation and deacetylation, 
 * - and PIC assembly and disassembly.
 *
 */
typedef struct GillespieRates GillespieRates;
struct GillespieRates {
  float transport;         /* rates[1] */
  float transport_rate[NGENES];
//  int transport_operations;
  float mRNAdecay;         /* rates[2] */
  float mRNAdecay_rate[NGENES];
//  int mRNAdecay_operations;
  float pic_disassembly;    /* rates[3] */
  float pic_disassembly_rate[NGENES];
//  int pic_disassembly_operations;   
  float acetylation;
  float acetylation_rate[NGENES]; 
  float deacetylation;
  float deacetylation_rate[NGENES];  
  float pic_assembly;
  float pic_assembly_rate[NGENES];  
  int transcript_init; 
  int transcript_init_rate[NGENES];    

  /* subtotal, including above, but not including konrate */
  float subtotal;
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         /* transcription factor */
  float Koff;        /* replacing hamming_dist */
  int BS_pos;        /* start position of BS on DNA, always with reference to forward strand */                     
  int N_hindered;      /* the number of BS hindered by this TF when it binds to the current BS */  
};

typedef struct Genotype Genotype;
struct Genotype {
/* directly subjected to mutation*/
    int ngenes;                                           /* the number of actual loci */
    int ntfgenes;                                         /* the number of actual tf loci */
    int nproteins;                                        /* because of gene duplication, the number of proteins and mRNAs can be different 
                                                           * form the number of genes. This nprotein is the number of elements in protein_pool*/
    int *protein_pool[NPROTEINS][2];                      /* element 1 record how many genes/mRNAs producing this protein, 
                                                           * ele 2 is which genes/mRNAs*/
    int which_protein[NGENES];                            /* in case of gene duplication, this array tells which protein given a gene id */
   
    char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN];
    char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN];
    char tf_seq_rc[TFGENES][MAX_COPIES][TF_ELEMENT_LEN];  /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
    int N_act;                                            /* number of activating TF*/
    int N_rep;                                            /* number of repressing TF*/
    int activating[NPROTEINS][MAX_COPIES];                /* 1 is activating TF, 0 is repressing */ 

    float mRNAdecay[NGENES];                              /* kinetic rates*/
    float proteindecay[NGENES];                           /* kinetic rates*/
    float translation[NGENES];                            /* kinetic rates*/   
    float pic_disassembly[NGENES][MAX_COPIES];            /* kinetic rates*/
                           
 /* binding sites related data, not directly subjected to mutation*/
    int re_calc[NGENES][4];                               /* If there is no mutation to a duplicated gene, don't calc the distribution
                                                           * for it. Element 1: where to copy the distribution from, -1 means re_calc. 
                                                           * Ele 2: 1 for the binding distribution can be copied from this gene.
                                                           * Ele 3: 1 for re_calc the bing sites.
                                                           * Ele 4: 1 for update binding sites info for this gene in clone_cell */    
    int binding_sites_num[NGENES];                        /* total number of binding sites */
    int max_hindered_sites[NGENES];                       /* maximal number of BS a BS can hinder*/ 
    int N_act_BS[NGENES];                                 /* total number of binding sites of activating TF */
    int N_rep_BS[NGENES];                                 /* total number of binding sites of repressing TF */
//    int *N_configurations[NGENES];                        /* maximal numbers of activators bound given n rep bound */ 
//    int max_N_rep_bound[NGENES];                          /* maximal number of repressors bound to a promoter */ 
//    int max_N_act_bound[NGENES];
    AllTFBindingSites *all_binding_sites[NGENES];      

    float fitness;
//    float death;
//    float avg_dt;
//    float max_dt;
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
    float cell_size;                    /* size of cell */
    float growth_rate;                  /* total growth rate in the previous deltat */
    int mRNA_cyto_num[NGENES];          /* mRNAs in cytoplasm */
    int mRNA_nuclear_num[NGENES];       /* mRNAs in nucleus */
    int mRNA_transl_cyto_num[NGENES];   /* mRNAs are in the cytoplasm, but only recently */
    int mRNA_transcr_num[NGENES][MAX_COPIES];  /* mRNAs which haven't finished transcription yet */

    FixedEvent *mRNA_transl_time_end;   /* times when mRNAs become fully loaded with ribosomes and start producing protein */
    FixedEvent *mRNA_transl_time_end_last;  
    FixedEvent *mRNA_transcr_time_end;  /* times when transcription is complete and an mRNA is available to move to cytoplasm */
    FixedEvent *mRNA_transcr_time_end_last;
    FixedEvent *env0_time_end;          /* times when env=0 ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.*/
    FixedEvent *env0_time_end_last;
    FixedEvent *env1_time_end;
    FixedEvent *env1_time_end_last;   

    float Pact[NGENES];
    float protein_conc[NPROTEINS];     /* pooled protein concentration from gene_specific_protein_conc */
    float gene_specific_protein_conc[NGENES]; /* stores the "protein" concentration for each gene.
                                               * can be considered temporary data. Make muation easier to
                                               * deal with. */
    float konvalues[NGENES][2];        /* moved from KonState*/  
    int active[NGENES][MAX_COPIES];       /* gives the state of each of the genes, according to figure
                                          *  NUC_NO_PIC = 0,
                                          *  NO_NUC_NO_PIC =1,
                                          *  PIC_NO_NUC = 3,
                                          */ 
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
//const float PROB_ACTIVATING;
const float TRANSCRIPTINIT; 
const float DEACETYLATE;
const float ACETYLATE;
const float PICASSEMBLY;
//const float STARTNUCLEUS;
const float KR;    
//const float GASCONSTANT;
//const float COOPERATIVITY;
//const float COOPERATIVE_DISTANCE; 
const float NUMSITESINGENOME ;
const float mN ;     



/* see netsim.c for documentation for these global variables */
float kon; 
//float kon_after_burnin; 
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
float Koff[TF_ELEMENT_LEN-4+1];



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
                                float []) ;

extern void initialize_genotype_fixed(Genotype *,
                                      float *);

extern void calc_all_binding_sites_copy(Genotype *, int, int *);

extern void calc_all_binding_sites(Genotype *);

extern float calc_ratio_act_to_rep(AllTFBindingSites *,
                                    int ,
                                    int ,
                                    int ,
                                    int ,
                                    int , 
                                    int [NGENES][MAX_COPIES],                                    
                                    int ,
                                    int ,
                                    int *,
                                    float [NGENES]);

extern float calc_ratio_act_to_rep_approximation(AllTFBindingSites *,
                                    int ,
                                    int ,
                                    int ,
                                    int ,
                                    int , 
                                    int [NGENES][MAX_COPIES],                                    
                                    int ,
                                    float [NGENES]);

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
                            int,  
                            int *[NPROTEINS][2],
                            float [NGENES],
                            float [NGENES],
                            float [NPROTEINS],
                            long int *);
//
//extern void initialize_cell_cache(CellState *,
//                                  Genotype,                                 
//                                  float **,
//                                  int,
//                                  int);

//extern void calc_time (float, 
//                       float, 
//                       GillespieRates *,                      
//                       float *, 
//                       float *);

//extern void calc_kon_rate(float,
//                          KonStates *,
//                          float *);

extern void change_mRNA_cytoplasm(int,
                                  Genotype *,
                                  CellState *,
                                  GillespieRates *);

//extern void calc_koff(int,
//                      AllTFBindingSites *,
//                      CellState *,
//                      float *,
//                      float);

//extern void scan_nearby_sites(int,
//                              AllTFBindingSites *,
//                              CellState *,
//                              GillespieRates *,
//                              float *,
//                              float);
//
//extern void remove_kon(int,
//                       int,
//                       GillespieRates *,
//                       float,
//                       KonStates *,
//                       float);
//
//extern void add_kon(float,
//                    float,
//                    int,
//                    int,
//                    GillespieRates *,
//                    KonStates *);

//extern int ready_to_transcribe(int,
//                               int,
//                               int *,
//                               int,
//                               AllTFBindingSites *,
//                               int [NGENES][MAX_COPIES],
//                               int *);
//
//extern int is_one_activator(int,
//                            int,
//                            int *,
//                            int,
//                            AllTFBindingSites *,
//                            int [NGENES][MAX_COPIES]);

//extern void calc_from_state(Genotype *,
//                            CellState *,
//                            GillespieRates *,                            
//                            float [],
//                            float []);

extern int does_fixed_event_end(FixedEvent *,
                                FixedEvent *,
                                FixedEvent *,
                                FixedEvent *,
                                float);

//extern void calc_dt(float *,
//                    float *,
//                    GillespieRates *, 
//                    Genotype *,
//                    float [],                    
//                    int [],
//                    int [],
//                    float [],
//                    float [],
//                    int);

extern void end_transcription(float *,
                              float,
                              CellState *,
//                              float [NGENES],
                              GillespieRates *,
                              int);


extern void disassemble_PIC(Genotype *,
                            CellState *,                                                        
                            GillespieRates *);

extern void revise_activity_state(int,
                                  int,
                                  Genotype *,
                                  CellState *,
                                  GillespieRates *);

//extern void remove_tf_binding(Genotype *,
//                              CellState *,
//                              GillespieRates *,
//                              KonStates *,
//                              int,
//                              float [],
//                              float);
//
//extern void attempt_tf_binding(Genotype *,
//                               CellState *,
//                               GillespieRates *,
//                               float **,
//                               KonStates *,
//                               int *,
//                               int *,
//                               int,
//                               float);

//extern void add_time_points(float,
//                            float [NPROTEINS],
//                            TimeCourse **,
//                            TimeCourse **);
//
//extern void add_integer_time_points(float,
//                                    int [NPROTEINS],
//                                    TimeCourse **,
//                                    TimeCourse **);

extern float compute_tprime(float, float, float, float);

extern float compute_integral(float, float, float, float, float, float, float, float);

extern float compute_growth_rate_dimer(float *,
                                       Genotype *,
                                       CellState *,
                                       float,
                                       float,
                                       float, 
                                       float,
                                       float,
                                       float,
				       int );

extern void update_protein_conc_cell_size(Genotype *,
                                          CellState *,
                                          GillespieRates *,
                                          float,                                        
                                          float,
//                                          TimeCourse **,
//                                          TimeCourse **,                                         
					  int);

//extern void calc_num_bound(float[],
//                           int );

//extern int sum_rate_counts(int[MAX_COPIES]);
//
//extern void get_gene(int [MAX_COPIES], int, int *, int *);

extern void transport_event(Genotype *,
                            CellState *,
                            GillespieRates *,
//                            TimeCourse **, 
//                            TimeCourse **, 
                            float,
                            float,
                            long int *);

//extern void tf_binding_event(GillespieRates *, CellState *, Genotype *, 
//                             KonStates *, float *, 
////							 TimeCourse **, TimeCourse **,
//                             float, float, float, int, int, int,int *);
//
//extern void tf_unbinding_event(GillespieRates *, CellState *, Genotype *, 
//                               KonStates *, float *, 
////							   TimeCourse **, TimeCourse **,
//                               float, float, float, float, int, int *,int *);

extern void mRNA_decay_event(GillespieRates *, CellState *, Genotype *, long int *);
//							  TimeCourse **, TimeCourse **,
//                             float, float, float, int *);

extern void histone_acteylation_event(GillespieRates *, CellState *, Genotype *, long int *); 
                                       
//									  TimeCourse **, TimeCourse **,
//                                      float, float);

extern void histone_deacteylation_event(GillespieRates *, CellState *, Genotype *, long int *); 
                                        
//										 TimeCourse **, TimeCourse **,
//                                        float, float, int *);

extern void assemble_PIC_event(GillespieRates *, CellState *, Genotype *, long int *); 
                                
//							   TimeCourse **, TimeCourse **,
//                               float, float, int *);

extern void disassemble_PIC_event(Genotype *, CellState *, GillespieRates *, long int * 
//								  TimeCourse **, TimeCourse **,
                                  );

extern void transcription_init_event(GillespieRates *, CellState *, Genotype *,
                                     
//									  TimeCourse **, TimeCourse **,
                                     float, float,long int *);

//extern void shift_binding_site_ids(CellState *, 
//                                   
//                                   int,
//                                   int);


//extern void recompute_kon_rates(GillespieRates *,
//                                CellState *,
//                                Genotype *,
//                                KonStates *,
//                                int);
//
//extern void recompute_koff_rates(GillespieRates *,
//                                 CellState *,
//                                 Genotype *,
//                                 float *,
//                                 float);
//
//extern void recalibrate_cell(GillespieRates *,
//                             CellState *,
//                             Genotype *,                            
//                             float **,
//                             float [NGENES],
//                             float [NGENES],
//                             float);

extern int do_single_timestep(Genotype *, 
                               CellState *,
                               GillespieRates *,                              
                               float *,                    
                               int,
                               int,                                                             
                               int *,
                               long int *) ;
							   
extern void free_fixedevent(CellState *);
 
extern float calc_avg_growth_rate(Genotype *, 
                                    CellState *, 
                                    float [NGENES],
                                    float [NGENES],
                                    GillespieRates *,
                                    float ,
                                    float ,
                                    long int *); 
                                    
extern void try_fixation(Genotype *, Genotype *, int *, int *, long int *);

extern int mutate(Genotype *, float [NUM_K_DISASSEMBLY],long int *);
  
extern void init_run_pop(//Genotype [N_para_threads+1],
//                         CellState [N_para_threads+1],
//                         TimeCourse *[2][NGENES],
//                         TimeCourse *[2][NGENES], 
//                         float, /* in Kelvin */
                         float [NUM_K_DISASSEMBLY],
//                         int,
                         FILE *);


extern void print_time_course(TimeCourse *,
                              int, int);

extern void print_all_protein_time_courses(TimeCourse *[2][NPROTEINS],
                                          TimeCourse *[2][NPROTEINS]);
                                          
extern void clone_cell(Genotype *,                
                	Genotype *,
                        int);

extern void log_snapshot(Genotype *,
                         CellState *,
                         GillespieRates *,
                         float ,
                         float );

extern int mod(int, int); //thanks to the stupid % in gcc

extern void calc_all_rates(Genotype *,
                            CellState *,
                            GillespieRates *,                     
                            int);

extern void end_translation(Genotype *, CellState *, GillespieRates *, float *, float );

extern void do_fixed_event(Genotype *, 
                            CellState *, 
                            GillespieRates *, 
                            float *,
                            float ,
                            int , 
                            int *);

extern int do_Gillespie_event(Genotype*, CellState *, GillespieRates *, float, float, long int *);

extern void calc_configurations(Genotype *, int);

extern void susbtitution(Genotype *,long int *);

extern void insertion(Genotype *,long int *);

extern void partial_deletion(Genotype *,long int *);

extern void whole_gene_deletion(Genotype *,long int *);

extern void gene_duplicaton(Genotype *,long int *);

extern void mut_binding_sequence(Genotype *,long int *);

extern void mut_kinetic_constant(Genotype *, float [NUM_K_DISASSEMBLY],long int *);

extern void draw_mutation(int, char *,long int *);

extern void initialize_cache(Genotype *);

#endif /* !FILE_NETSIM_SEEN */
