/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-  */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN
#endif

#include <stdio.h>
#include "RngStream.h"

/*Simulation mode*/
#define JUST_PLOTTING 0
#define CONTINUE 0
#define NEUTRAL 0
#define RUN_FULL_SIMULATION 1
#define SKIP_INITIAL_GENOTYPE 0
#define SET_BS_MANUALLY 0
#define QUICK_BURN_IN 0

/*Runtime control*/
#ifndef MAX_MUT_STEP         
#define MAX_MUT_STEP 10000
#ifndef BURN_IN
#define BURN_IN 0
#endif
#define N_THREADS 1
#define N_REPLICATES 40
#define OUTPUT_INTERVAL 10

/*Miscellaneous settings*/
#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1.0e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6
#define OUTPUT_RNG_SEEDS 0
#define TIME_INFINITY 9.99e10
#define CAUTIOUS 0

/*Biology and evolution settings*/
#define DIRECT_REG 1
#define RANDOM_INIT_KINETIC_CONST 1
#define RANDOM_COOPERATION_LOGIC 0
#define RdcPdup 0
#define N_SIGNAL_TF 2 // the 1st TF enables basal activity in TFN. The 2nd is the actual signal TF. Planning on adding a 3rd tf as another signal TF. 
#define NO_REGULATION_COST 0
#define NO_REGULATION 0 // this locks the state of transcription factors to NUC_NO_PIC
#define ADJUST_FITNESS 0 // allows manually adjust the fitness of a phenotype
#if ADJUST_FITNESS
#define ADJUST_FITNESS_CONDITION 3 // more or less than 3 copies of regular transcription factor genes receives penaly in fitness
#define ADJUSTMENT 1.0e-3 //selection coefficient
#endif
#define UNLIMITED_MUTATION 0
#define IGNORE_BS_OVERLAPPING 0
#define SIMPLE_SUBSTITUTION 1
#define RANDOMIZE_SIGNAL2 0
#define ALPHA 0.2
#define MAX_RECALC_FITNESS 5
#define RULE_OF_REPLACEMENT 4 /* 0 for z-score, 1 for Wilcoxon, 2 for larger-fitness-fixes, 3 for larger-than-epsilon-fixes, 4 for s>minimal_selection_coefficient */
#if RULE_OF_REPLACEMENT==4
#define minimal_selection_coefficient 1.0e-8
#endif
#ifndef MAX_COPIES
#define MAX_COPIES 2       /* each gene can have at most two copies*/
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
#else                       /* otherwise, by default assuming selection gene is not a TF */
    #ifndef TFGENES             /* number of genes encoding TFs */
    #define TFGENES 20         /* the initial value is set in initiate_genotype*/
    #endif
    #ifndef NGENES
    #define NGENES 25  /* total number of genes: add the (non-TF) selection gene to the total (default case) */
    #endif
    #ifndef NPROTEINS           
    #define NPROTEINS 25
    #endif
    #ifndef EFFECTOR_GENES
    #define EFFECTOR_GENES 5  /* this is the upper limit of effector gene copies*/
    #endif
#endif
#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 8      /* length of binding element on TF */
#define NMIN 6
#define NUM_K_DISASSEMBLY 133 /* number of differents for PIC disassembly from data file  */
#ifndef HIND_LENGTH
#define HIND_LENGTH 6         /* default length of hindrance on each side of the binding site (original was 6) */
                              /* the binding of Lac repressor blockes 12 bp. Record MT 1981*/
#endif
#define MAX_MODE 8  /* MAX_MODE is the max number of tf that can bind to a promoter plus 1*/
#define MAX_BS_IN_CLUSTER 100


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
enum { KON_PROTEIN_DECAY_INDEX = 0, KON_SALPHAC_INDEX = 1 };

/*
 * enum for 'CellState'->active indices
 */

enum { NUC_NO_PIC = 0,
       NO_NUC_NO_PIC =1,
       PIC_NO_NUC = 3,};

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
//  float transport;         /* rates[1] */
//  float transport_rate[NGENES];
  float mRNAdecay;         /* rates[2] */
  float mRNAdecay_rate[NGENES];
  float pic_disassembly;    /* rates[3] */
  float pic_disassembly_rate[NGENES]; 
  float acetylation;
  float acetylation_rate[NGENES]; 
  float deacetylation;
  float deacetylation_rate[NGENES];  
  float pic_assembly;
  float pic_assembly_rate[NGENES];  
  int transcript_init; 
  int transcript_init_rate[NGENES];   
  float subtotal;
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         /* transcription factor */
  float Koff;        /* replacing hamming_dist */
  int mis_match;     /* number of mismatched nuc */
  int BS_pos;        /* start position of BS on DNA, always with reference to forward strand */                     
  int N_hindered;    /* the number of BS hindered by this TF when it binds to the current BS */  
//  int compressed;    /* if this site has been compressed, set the value to 1 */ 
};

typedef struct CompressedBindingSites CompressedBindingSites;
//struct CompressedBindingSites
//{
//    int tf_id;
//    double coeff_on_flux[MAX_MODE]; /* coefficient of on flux of mode 1-5*/
//    double off_flux;
//    int N_hindered;
//    int start_pos;
//    int end_pos;
//};

typedef struct Genotype Genotype;
struct Genotype {

    int ngenes;                                             /* the number of loci */
    int ntfgenes;                                           /* the number of tf loci */
    int nproteins;                                          /* because of gene duplication, the number of proteins and mRNAs can be different 
                                                               from the number of loci. nprotein is the number of elements in protein_pool*/
    int ncopies_under_penalty;                              /* this is the number of gene copies that exceed the MAX_COPIES */  
                                                               
    int which_protein[NGENES];                              /* in case of gene duplication, this array tells the protein corresponding to a given gene id */   
    char cisreg_seq[NGENES][CISREG_LEN];    
    
    /*these apply to protein, not loci*/
    int N_act;                                              /* number of activators*/
    int N_rep;                                              /* number of repressors*/    
    int activating[NPROTEINS];                              /* 1 for activator, 0 for repressor, -1 for non-tf */ 
    char tf_seq[NPROTEINS][TF_ELEMENT_LEN];
    char tf_seq_rc[NPROTEINS][TF_ELEMENT_LEN];                /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
//    float P_dup_per_protein[NPROTEINS];                     /* probability of duplication of the copies for each protein */
    int protein_pool[NPROTEINS][2][NGENES];                 /* element 1 record how many genes/mRNAs producing this protein,ele 2 stores which genes/mRNAs*/
    float koff[NPROTEINS];                                     /* kinetic rates*/ 
    
    /*these apply to loci*/
    float mRNAdecay[NGENES];                                /* kinetic rates*/
    float proteindecay[NGENES];                             /* kinetic rates*/
    float translation[NGENES];                              /* kinetic rates*/   
    float pic_disassembly[NGENES];                          /* kinetic rates*/    
    int min_act_to_transc[NGENES];                          /* 1 for OR GATE, at leat 2 FOR AND GATE */ 
//    int N_failed_dup_and_del[NGENES];                       /* If too many mutations has been wasted on the gene, */
                                                            /* reduce the rate of dup and del on this gene by 10 fold*/
 
    /* binding sites related data, applying to loci*/   
    int cisreg_cluster[NGENES][NGENES];                     /* For genes having the same cis-reg, tf distribution can be shared.
                                                             * Genes having the same cis-reg are clustered.
                                                               1st dim stores cluster ids, 2nd dim stores gene_ids in a cluster.
                                                               cisreg_cluster works with which_cluster*/
    int which_cluster[NGENES];                              /* which_cluster stores the cluster id of a gene*/                                                           
    int recalc_TFBS[NGENES];                                /* whether to recalc the TFBS*/
    int clone_info[NGENES];                                 /* whether to copy info back to this gene in clone_cell*/    
    int binding_sites_num[NGENES];                          /* total number of binding sites */
    int max_hindered_sites[NGENES];                         /* maximal number of BSs a BS can hinder*/ 
    int max_hindered_clusters[NGENES];                      /* maximal number of clusters a cluster can hinder*/   
//    int cluster_num[NGENES];                              /* number of clusters of compressed binding sites */
    int N_act_BS[NGENES];                                   /* total number of binding sites of activating TF */
    int N_rep_BS[NGENES];                                   /* total number of binding sites of repressing TF */
    int avg_N_BS_in_cluster[NGENES][NGENES][2];
//    int *N_configurations[NGENES];                        /* maximal numbers of activators bound given n rep bound */ 
//    int max_N_rep_bound[NGENES];                          /* maximal number of repressors bound to a promoter */ 
//    int max_N_act_bound[NGENES];
    AllTFBindingSites *all_binding_sites[NGENES];   
//    CompressedBindingSites *compressed_binding_sites[NGENES];

    float avg_GR1;
    float avg_GR2;
    float sq_SE_GR1;
    float sq_SE_GR2;
    float fitness;
    float sq_SE_fitness;
    float fitness_measurement[MAX_RECALC_FITNESS*N_REPLICATES];
    int N_ffls[4];     
    int N_non_ffl_motifs;
    float proportion_c1ffls;  
    int N_act_genes;   
};

/* 
 * transcription/translation delays are sorted linked lists.  Deleting
 * the head each time, and tack new stuff on the end.  Linked lists
 * are easy to create pre-sorted.
 */
typedef struct FixedEvent FixedEvent;
struct FixedEvent {
  int gene_id; 
  float time;
  FixedEvent *next;
};

typedef struct CellState CellState;
struct CellState {   
    float cell_size;                    /* size of cell */
    float cell_size_after_burn_in;          
    float growth_rate;                  /* total growth rate in the previous deltat */   
    int mRNA_cyto_num[NGENES];          /* mRNAs in cytoplasm */
//    int mRNA_nuclear_num[NGENES];       /* mRNAs in nucleus */
    int mRNA_transl_cyto_num[NGENES];   /* mRNAs are in the cytoplasm, but only recently */
    int mRNA_transcr_num[NGENES];       /* mRNAs which haven't finished transcription yet */

    FixedEvent *mRNA_transl_init_time_end;   /* times when mRNAs become fully loaded with ribosomes and start producing protein */
    FixedEvent *mRNA_transl_init_time_end_last;  
    FixedEvent *mRNA_transcr_time_end;  /* times when transcription is complete and an mRNA is available to move to cytoplasm */
    FixedEvent *mRNA_transcr_time_end_last;
    FixedEvent *signalB_starts_end;          /* times when env=A ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.*/
    FixedEvent *signalB_starts_end_last;   
    FixedEvent *signalA_starts_end;          /* times when env=A ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.*/
    FixedEvent *signalA_starts_end_last; 
    FixedEvent *burn_in_growth_rate;
    FixedEvent *burn_in_growth_rate_last;  
    FixedEvent *sampling_point_end;
    FixedEvent *sampling_point_end_last; 

    float t_to_update_Pact;
    float interval_to_update_Pact; 
    float Pact[NGENES];
    float last_Pact[NGENES];
//    float equilibrium_tf_conc[NPROTEINS]; /* this is the nuclear concentration of free tf when binding/unbinding to non-specific sites reach equilibrium */  
//    float protein_conc[NPROTEINS];        /* this is the concentration of proteins in cell. For tf, this stores the nuclear concentration, for selection proteins, this is cytosol concentration*/
    float protein_number[NPROTEINS];     /* pooled protein number from gene_specific_protein_conc */
    float gene_specific_protein_number[NGENES]; /* stores the "protein" number for each gene.
                                               * can be considered temporary data. Make muation easier to
                                               * deal with. */
    float konvalues[NGENES][2];        /* moved from KonState*/  
    int active[NGENES];       /* gives the state of each of the genes, according to figure
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

typedef struct Mutation Mutation;
struct Mutation
{
    char mut_type;
    int pos_g;
    int pos_n;
    char nuc_diff[3];    /*the first three elements store the nuc after mutation (max_indel=3)*/
    int kinetic_type;    /*0 for pic_disassembly, 1 for mRNA_decay, 2 for translation, 3 for protein_decay*/
    float kinetic_diff;
};

/*
 * global variables
 */

/* see netsim.c for documentation for these global constant variables */
int MAXELEMENTS; 
const int MAXBOUND;
const float KRNA;
const float TTRANSLATION;
const float TTRANSCRIPTION;
const float TRANSCRIPTINIT; 
const float DEACETYLATE;
const float ACETYLATE;
const float PICASSEMBLY;
//const float KR;
const float NUMSITESINGENOME;

/* see netsim.c for documentation for these global variables */
//float lumped_kon;
float tdevelopment; 
float growth_rate_scaling; 
float Pp_s;       
float h;
float h_extra_copy;
float gmax_a;
float gmax_b;
float gpeak;
//float protein_aging;
int current_ploidy;
int init_TF_genes;
float penalty_of_extra_copies;
float reduction_in_P_dup;
float signal_strength;
float basal_signal_strength;
float background_signal_strength;
float env1_t_signalA;    
float env1_t_signalB;     
float env2_t_signalA;
float env2_t_signalB;
int env1_signalA_as_noise;    
int env2_signalA_as_noise;
int env1_signalA_mismatches; 
int env2_signalA_mismatches;
float env1_occurence;
float env2_occurence;
int init_N_act;
int init_N_rep;
int recalc_new_fitness;
char init_env1;
char init_env2;
int min_act_to_transcr_selection_protein;
float cost_term;
float penalty;
int N_replicates;
float duration_of_burn_in_growth_rate;
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

extern void initialize_sequence(char [], int, int, RngStream);

extern void print_genotype(Genotype *, int);

extern void print_all_binding_sites(int [NGENES],
                                    AllTFBindingSites *, 
                                    int ,
                                    char [TFGENES][TF_ELEMENT_LEN],
                                    char [NGENES][CISREG_LEN]
//                                    int [NGENES][2]
									);

extern void print_tf_occupancy(CellState *,
                               AllTFBindingSites *,
                               float);

extern void initialize_genotype(Genotype *, 
                                float [],
                                RngStream ) ;

extern void initialize_genotype_fixed(Genotype *,
                                      float *,
                                      RngStream );

extern void calc_all_binding_sites_copy(Genotype *, int);

extern void calc_all_binding_sites(Genotype *);

extern void set_binding_sites(Genotype *,int *, int (*)[5], int (*)[5], float (*)[5]);

extern void cluster_BS_cluster(Genotype *, int);

extern double calc_flux(AllTFBindingSites *,int,int,double [MAX_MODE]);

extern float calc_ratio_act_to_rep(AllTFBindingSites *,
                                    int ,
                                    int ,
                                    int ,
                                    int ,
                                    int , 
                                    int [NGENES],                                    
                                    int ,
                                    int ,
                                    int *,
                                    float [NGENES]);

extern float calc_TF_dist_from_compressed_BS(   CompressedBindingSites *,
                                                    int ,
                                                    int ,
                                                    int ,             
                                                    int [NPROTEINS], 
                                                    float [NGENES]);
extern float calc_TF_dist_from_all_BS(  AllTFBindingSites *,
                                        int ,
                                        int ,
                                        int ,
                                        int ,
                                        int , 
                                        int [NPROTEINS],                                    
                                        int ,
                                        float [NGENES],
                                        int );

extern float calc_TF_dist_from_all_BS_simple(   AllTFBindingSites *,
                                                int,
                                                int,
                                                int [NPROTEINS],                                    
                                                float [NGENES],                                               
                                                int);

extern int add_fixed_event(int,                           
                           float,
                           FixedEvent **,
                           FixedEvent **);

extern void add_time_point(float,
                           float,
                           TimeCourse **,
                           TimeCourse **);

extern void add_fixed_event_end(int,                                
                                float,
                                FixedEvent **,
                                FixedEvent **);

extern void delete_fixed_event(int,                               
                               int,
                               FixedEvent **,
                               FixedEvent **);

extern void delete_fixed_event_start(FixedEvent **,
                                     FixedEvent **);

extern void initialize_cell(CellState *,
                            int,
                            int,  
                            int [NPROTEINS][2][NGENES],
                            float [NGENES],
                            int [NGENES],
                            float [NPROTEINS],
                            RngStream);

extern void change_mRNA_cytoplasm(int,
                                  Genotype *,
                                  CellState *,
                                  GillespieRates *);

extern int does_fixed_event_end(CellState*,
                                float);

extern int does_fixed_event_end_plotting(   FixedEvent *,
                                            FixedEvent *,
                                            FixedEvent *,
                                            FixedEvent *,
                                            FixedEvent *,
                                            FixedEvent *,
                                            float,
                                            float);

extern void update_protein_number_cell_size(Genotype *,
                                            CellState *,
                                            GillespieRates *,
                                            float,                                        
                                            float,
                                            char,
                                            int *,
                                            char *,
                                            int,
                                            Mutation *);

extern void end_transcription(float *,
                                float,
                                CellState *,
                                GillespieRates *,
                                Genotype *,
                                int *,
                                char *,
                                int,
                                Mutation *,
                                char *);


extern void disassemble_PIC(Genotype *,
                            CellState *,                                                        
                            GillespieRates *);

extern void revise_activity_state(int,
                                  int,
                                  Genotype *,
                                  CellState *,
                                  GillespieRates *);

extern float compute_tprime(Genotype*, CellState*, float*, float);

extern float compute_integral(Genotype *, CellState *, float *, float, float);

extern float compute_growth_rate_dimer( float *,
                                        Genotype *,
                                        CellState *,
                                        float*,
                                        float,
                                        float,                                       
                                        char,
                                        int *,
                                        char *,
                                        int,
                                        Mutation *);

extern void transport_event(Genotype *,
                            CellState *,
                            GillespieRates *,
                            float,
                            float,
                            RngStream);

extern void mRNA_decay_event(GillespieRates *, CellState *, Genotype *, float, float, RngStream);

extern void histone_acteylation_event(GillespieRates *, CellState *, Genotype *, RngStream);

extern void histone_deacteylation_event(GillespieRates *, CellState *, Genotype *, RngStream); 

extern void assemble_PIC_event(GillespieRates *, CellState *, Genotype *, RngStream);

extern void disassemble_PIC_event(Genotype *, CellState *, GillespieRates *, RngStream);

extern void transcription_init_event(GillespieRates *, CellState *, Genotype *, float, float, RngStream);

extern void do_single_timestep( Genotype *, 
                                CellState *,
                                GillespieRates *,                              
                                float *,                                                  
                                char *,
                                float,
                                float,  
                                int,
                                int,
                                RngStream,
				char *,
                                char *,
                                int,
                                Mutation *,
                                int *,int*) ;

extern void do_single_timestep_plotting(    Genotype *, 
                                            CellState *,
                                            GillespieRates *,                              
                                            float *,                                                  
                                            char *,
                                            float,
                                            float,
                                            int,                                           
                                            float (*)[159],
                                            float [159],
                                            RngStream,                              
                                            int *,
                                            int *) ;
							   
extern void free_fixedevent(CellState *);
 
extern void calc_avg_growth_rate(   Genotype *,                                    
                                    int [NGENES],
                                    float [NPROTEINS],                                  
                                    RngStream [N_THREADS],
				    char *,
                                    char *,
                                    int,
                                    float *,
                                    float *,
                                    Mutation *); 
  
extern int init_run_pop(float [NUM_K_DISASSEMBLY], char*, char *, char *, char *, char *, char *, char *, unsigned long int [6]);

extern void print_time_course(TimeCourse *, FILE *);

extern void print_all_protein_time_courses(int, TimeCourse **, TimeCourse **, FILE *);

extern void log_snapshot(Genotype *,
                         CellState *,
                         GillespieRates *,
                         float ,
                         float );

extern int mod(int, int); 

extern void calc_all_rates(Genotype *,
                            CellState *,
                            GillespieRates *, 
                            float,
                            int);

extern void end_translation_init(   Genotype *, 
                                    CellState *,    
                                    GillespieRates *, 
                                    float *, 
                                    float,
                                    int *,
                                    char *,
                                    int,
                                    Mutation *,
                                    char *);

extern int do_fixed_event(  Genotype *, 
                            CellState *, 
                            GillespieRates *, 
                            float *,
                            float ,
                            int , 
                            float ,
                            float,
                            char *,
                            int,                         
                            int *,
                            char *,
                            int,
                            Mutation *);

extern int do_fixed_event_plotting( Genotype *, 
                                    CellState *, 
                                    GillespieRates *, 
                                    float *,
                                    float ,
                                    int , 
                                    float ,
                                    float,        
                                    char *,                                    
                                    int ,                                    
                                    int *);

extern void do_Gillespie_event(Genotype*, CellState *, GillespieRates *, float, float, RngStream, int *, char *, int, Mutation *);

extern void calc_configurations(Genotype *, int);

extern int mutate(Genotype *, float [NUM_K_DISASSEMBLY], RngStream, Mutation *);

extern void mut_susbtitution(Genotype *, Mutation *, RngStream);

extern void mut_insertion(Genotype *,Mutation *, RngStream);

extern void mut_partial_deletion(Genotype *,Mutation *, RngStream);

extern int mut_whole_gene_deletion(Genotype *,Mutation *, RngStream);

extern int mut_duplicaton(Genotype *,Mutation *, RngStream);

extern void mut_binding_sequence(Genotype *,Mutation *, RngStream);

extern void mut_kinetic_constant(Genotype *, Mutation *, float [NUM_K_DISASSEMBLY], RngStream);

extern void mut_identity(Genotype *, Mutation *, RngStream);

extern void mut_koff(Genotype *, Mutation *, RngStream);

extern int reproduce_mutate(Genotype *, Mutation *,RngStream);

extern void reproduce_susbtitution(Genotype *, Mutation *);

extern void reproduce_insertion(Genotype *,Mutation *);

extern void reproduce_partial_deletion(Genotype *,Mutation *);

extern void reproduce_whole_gene_deletion(Genotype *,Mutation *);

extern void reproduce_gene_duplicaton(Genotype *,Mutation *);

extern void reproduce_mut_binding_sequence(Genotype *,Mutation *);

extern void reproduce_mut_kinetic_constant(Genotype *, Mutation *);

extern void reproduce_mut_identity(Genotype *, Mutation *);

extern void reproduce_mut_koff(Genotype *, Mutation *);

extern void draw_mutation(Genotype *, char *, RngStream);

extern void initialize_cache(Genotype *);

extern void update_protein_pool(Genotype *, int, int, char);

extern void update_cisreg_cluster(Genotype *, int, char);

extern void clone_cell_forward(Genotype *, Genotype *, int);

extern void clone_cell_backward(Genotype *, Genotype *, int);

extern float try_fixation(Genotype *, Genotype *, int, int, int *, RngStream);

extern void calc_avg_growth_rate_plotting(  Genotype *,                               
                                            int [NGENES],
                                            float [NPROTEINS], 
                                            RngStream [N_THREADS]); 

extern void summarize_binding_sites(Genotype *,int);

extern void print_binding_sites_distribution(Genotype *,int, int);

extern int check_concurrence(   float , 
                                FixedEvent *, 
                                FixedEvent *, 
                                FixedEvent *, 
                                FixedEvent *,
                                FixedEvent *,
                                float);

extern void set_env(CellState *, char, float, float);

extern void output_genotype(char *, char *, char *, char *, Genotype *, int);

extern void release_memory(Genotype *,Genotype *, RngStream *, RngStream [N_THREADS]);

extern void calc_fx_dfx(float, int, float*, float*, float*, float*, float*);

extern void resolve_overlapping_sites(Genotype *, int, int [NGENES]);

extern void evolve_N_steps(  Genotype *, 
                            Genotype *,
                            int *, 
                            int, 
                            int *,                 
                            char *,
                            char *,
                            char *,
                            char *, 
                            char *, 
                            char *, 
                            float [NPROTEINS],
                            int [NGENES],
                            float [NUM_K_DISASSEMBLY],
                            Mutation *, 
                            RngStream,
                            RngStream [N_THREADS],
                            int );

extern void run_simulation( Genotype *, 
                            Genotype *,
                            char *,
                            char *,
                            char *,
                            char *,
                            char *, 
                            char *, 
                            char *, 
                            float [NPROTEINS],
                            int [NGENES],
                            float [NUM_K_DISASSEMBLY],
                            Mutation *, 
                            RngStream,
                            RngStream [N_THREADS]);
extern void continue_simulation(Genotype *, 
                                Genotype *,
                                FILE *,
                                int,                               
                                char *,
                                char *,
                                char *,
                                char *,
                                char *, 
                                char *, 
                                char *, 
                                float [NPROTEINS],
                                int [NGENES],
                                float [NUM_K_DISASSEMBLY],
                                Mutation *, 
                                RngStream,
                                RngStream [N_THREADS]);
extern void run_plotting(   Genotype *,
                            Genotype *,
                            int [NGENES],
                            float [NGENES],
                            RngStream [N_THREADS],
                            Mutation *,
                            FILE *,
                            int);

extern void calc_fitness_stats( Genotype *,
                                float (*)[N_REPLICATES],
                                float (*)[N_REPLICATES],
                                int);


extern void evolve_neutrally(   Genotype *,
                                Genotype *,
                                char *,
                                Mutation *,
                                int,
                                float [NUM_K_DISASSEMBLY],
                                RngStream);

extern void replay_mutations(   Genotype *,
                                Genotype *,
                                FILE *,
                                Mutation *,
                                int,
                                RngStream);

extern void find_ffl(Genotype *);

#endif /* !FILE_NETSIM_SEEN */
