/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-  */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>
#include "RngStream.h"

/*Simulation mode*/
#define RUN_PLOTTING 0
#define PLOT_ALTERNATIVE_FITNESS 0
#define NEUTRAL 0
#define RUN_FULL_SIMULATION 1
#define EXTERNAL_SIGNAL 0

/*Runtime control*/  
#define MAX_MUT_STEP 400
#define BURN_IN 200
#define MAX_MUTATIONS 800000
#define MAX_TRIALS 2000
#define N_THREADS 10
#define N_REPLICATES 200
#define OUTPUT_INTERVAL 10
#define SAVING_INTERVAL 10
#define OUTPUT_MUTANT_DETAILS 0
#define N_TIMP_POINTS 90


/*Miscellaneous settings*/
#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1.0e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6
#define OUTPUT_RNG_SEEDS 1
#define TIME_INFINITY 9.99e10
#define TIME_OFFSET 0.01
#define CAUTIOUS 0

/*Biology and evolution settings*/
#define DIRECT_REG 0
#define NO_PENALTY 0
#define FORCE_OR_GATE 0
#if FORCE_OR_GATE
#define FORCE_MASTER_CONTROLLED 0
#endif
#define FORCE_DIAMOND 0
#define FORCE_SINGLE_FFL 0
#define RANDOM_COOPERATION_LOGIC 0
#define N_SIGNAL_TF 1 // number of input signals
#define NO_REGULATION_COST 0
#define NO_REGULATION 0 // this locks the state of TF genes to Rep
#define IGNORE_BS_OVERLAPPING 0
#define SIMPLE_SUBSTITUTION 1
#define RANDOMIZE_SIGNAL2 0
#define MAX_RECALC_FITNESS 10
#define minimal_selection_coefficient 1.0e-8

/* Because mutation can change the number of genes, the numbers defined here are used to allocate storage space only.
 * Set the numbers to be 8 folds of the initial ngenes and ntfgenes, so that we can have two whole genome duplications*/
#ifndef MAX_TF_GENES             /* number of genes encoding TFs */
#define MAX_TF_GENES 20         /* the initial value is set in initiate_genotype*/
#endif
#ifndef MAX_GENES
#define MAX_GENES 25  /* total number of genes: add the (non-TF) selection gene to the total (default case) */
#endif
#ifndef MAX_PROTEINS           
#define MAX_PROTEINS 25
#endif
#ifndef MAX_EFFECTOR_GENES
#define MAX_EFFECTOR_GENES 5  /* this is the upper limit of effector gene copies*/
#endif
#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 8      /* length of binding element on TF */
#define NMIN 6                /* minimal number of nucleotide that matches the binding sequence of a TF in its binding site*/
                              /* DO NOT MAKE NMIN<TF_ELEMENT_LEN/2, OTHERWISE calc_all_binding_sites_copy will make mistake*/  

#ifndef HIND_LENGTH
#define HIND_LENGTH 3         /* default length of hindrance on each side of the binding site,i.e. a tf occupies TF_ELEMENT_LEN+2*HIND_LENGTH */
                              /* the binding of Lac repressor blockes 12 bp. Record MT 1981*/
#endif
#define MAX_BINDING 10  /* MAX_MODE is the max number of tf that can bind to a promoter plus 1*/


/* 
 * define macros for logging warning/errors 
 */
#ifndef LOGGING_OFF
  #define LOG(...) { fprintf(fperror, "%s: ", __func__); fprintf (fperror, __VA_ARGS__) ; fflush(fperror); } 
#else
  #define LOG 
#endif

/*
 * enum for 'CellState'->transcriptional_state indices
 */
enum TRANSCRIPTIONAL_STATE {REPRESSED, INTERMEDIATE, ACTIVE};
enum PROTEIN_IDENTITY {ACTIVATOR=1, REPRESSOR=0, NON_TF=-1};
enum BOOLEAN {NA=-1, NO=0, YES=1};

/*
 * Rates for Gillespie algorithm
 */
typedef struct GillespieRates GillespieRates;
struct GillespieRates {
  float total_mRNA_decay_rate;         
  float mRNA_decay_rate[MAX_GENES];
  float total_active_to_intermediate_rate;    
  float active_to_intermediate_rate[MAX_GENES]; 
  float total_repressed_to_intermediate_rate;
  float repressed_to_intermediate_rate[MAX_GENES]; 
  float total_intermediate_to_repressed_rate;
  float intermediate_to_repressed_rate[MAX_GENES];  
  float total_intermediate_to_active_rate;
  float intermediate_to_active_rate[MAX_GENES];  
  int total_N_gene_transcript_initiated; 
  int transcript_initiation_state[MAX_GENES];   
  float total_Gillespie_rate;
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         
  float Kd;
  int mis_match;     /* number of mismatched nuc */
  int BS_pos;        /* start position of BS on DNA, always with reference to forward strand */                     
  int N_hindered;    /* the number of BS hindered by this TF when it binds to the current BS */  
};

/*
 * Selection contains two tests; each test states the signal whether 
 * expressing the effector is beneficial
 */
typedef struct Selection Selection;
typedef struct Test Test;
struct Test
{
    float t_development;
    float signal_on_strength;
    float signal_off_strength;
    float t_signal_on;
    float t_signal_off;
    char initial_effect_of_effector;
    int fixed_effector_effect;  
    float *external_signal;
    float duration_of_burn_in_growth_rate;
};

struct Selection
{
    Test test1;
    Test test2;
    float test1_weight;
    float test2_weight;
    float temporary_miu_ACT_TO_INT_RATE;
    float temporary_miu_protein_syn_rate;
    float temporary_miu_Kd;
    float temporary_DUPLICATION;
    float temporary_SILENCING;
    int temporary_N_effector_genes;
    int temporary_N_tf_genes;
};

/*Genotype*/
typedef struct Genotype Genotype;
struct Genotype {
    /*basic variables*/
    int ngenes;                                             /* the number of loci */
    int ntfgenes;                                           /* the number of tf loci */
    int nproteins;                                          /* because of gene duplication, the number of proteins and mRNAs can be different 
                                                               from the number of loci. nprotein is the number of elements in protein_pool*/  
    int nTF_families;

    int which_protein[MAX_GENES];                              /* in case of gene duplication, this array tells the protein corresponding to a given gene id */ 
    int which_TF_family[MAX_PROTEINS];                         /* We consider mutation to Kd does not create a new TF (matters when counting motifs), 
                                                             * but the calculation of TF binding distribution demands tfs with different Kd
                                                             * to be treated differently. We call tfs that differ only in Kd a tf family.
                                                             * We use this array and TF_family_pool to track which TF belongs which TF family.
                                                             */
    char cisreg_seq[MAX_GENES][CISREG_LEN];    
    
    /*these apply to protein, not loci*/
    int N_act;                                              /* number of activators*/ 
    int N_rep;                                              /* number of repressors*/ 
    int protein_identity[MAX_PROTEINS];                        /* 1 for activator, 0 for repressor, -1 for non-tf */ 
    char tf_seq[MAX_PROTEINS][TF_ELEMENT_LEN];                 /* concensus binding sequences*/
    char tf_seq_rc[MAX_PROTEINS][TF_ELEMENT_LEN];              /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
    int protein_pool[MAX_PROTEINS][2][MAX_GENES];                 /* element 1 record how many genes/mRNAs producing this protein,ele 2 stores which genes/mRNAs*/
    int TF_family_pool[MAX_PROTEINS][2][MAX_PROTEINS];                            
    float Kd[MAX_PROTEINS]; // Kd needs to be changed to locus specific.

    
    /*these apply to loci*/
    int locus_length[MAX_GENES];                               /* in codon, relates to transcriptional and translational delay */
    int total_loci_length;                                    
    float mRNA_decay_rate[MAX_GENES];                          /* kinetic rates*/
    float protein_decay_rate[MAX_GENES];                       /* kinetic rates*/
    float translation_rate[MAX_GENES];                         /* kinetic rates*/   
    float active_to_intermediate_rate[MAX_GENES];              /* kinetic rates*/    
    int min_N_activator_to_transc[MAX_GENES];                  /* 1 for OR GATE, at leat 2 FOR AND GATE */ 
 
    /* binding sites related data, applying to loci*/   
    int cisreg_cluster[MAX_GENES+1][MAX_GENES];                     /* For genes having the same cis-reg, tf distribution can be shared.
                                                             * Genes having the same cis-reg are clustered.
                                                               1st dim stores cluster ids, 2nd dim stores gene_ids in a cluster.
                                                               cisreg_cluster works with which_cluster*/
    int which_cluster[MAX_GENES];                              /* which_cluster stores the cluster id of a gene*/                                                           
    int recalc_TFBS[MAX_GENES];                                /* whether to recalc the TFBS*/    
    int binding_sites_num[MAX_GENES];                          /* total number of binding sites */    
    int max_unhindered_sites[MAX_GENES][3];                    /* maximal number of binding sites that do not hinder each other. element 1 for activator BS, 2 for repressor BS*/  
    int max_hindered_sites[MAX_GENES];                        /* maximal number of BSs a BS can hinder*/ 
    int N_act_BS[MAX_GENES];                                   /* total number of binding sites of activating TF */
    int N_rep_BS[MAX_GENES];                                   /* total number of binding sites of repressing TF */  
    AllTFBindingSites *all_binding_sites[MAX_GENES];   
    int N_allocated_elements;                               /* maximal number of AllTFBindingSites that have been allocated for each gene*/

    
    /*fitness related variable*/
    float avg_GR1;
    float avg_GR2;
    float sq_SE_GR1;
    float sq_SE_GR2;
    float fitness;
    float sq_SE_fitness;
    float fitness_measurement[MAX_RECALC_FITNESS*N_REPLICATES];
    
    /*Motifs related*/
    int N_motifs[39];  
    int TF_in_core_C1ffl[MAX_GENES][MAX_PROTEINS];
    int gene_in_core_C1ffl[MAX_GENES];
    int N_act_genes; 
    int N_act_genes_reg_by_env;
    int N_act_genes_not_reg_by_env;    
};

/* 
 * transcription/translation delays are sorted linked lists.  Deleting
 * the head each time, and tack new stuff on the end.  Linked lists
 * are easy to create pre-sorted.
 */
typedef struct FixedEvent FixedEvent;
struct FixedEvent {
  int event_id; 
  float time;
  FixedEvent *next;
};

typedef struct CellState CellState;
struct CellState { 
    float t;
    float cumulative_fitness;                    
    float cumulative_fitness_after_burn_in;          
    float instantaneous_fitness;                    
    int mRNA_aft_transl_delay_num[MAX_GENES];          /* mRNAs that have finished the translational delay */

    int mRNA_under_transl_delay_num[MAX_GENES];   /* mRNAs that are still under the translational delay (they do not contribute to protein 
                                                * turnover, but contribute to the cost of translation)  */
    int mRNA_under_transc_num[MAX_GENES];       /* mRNAs which haven't finished transcription */

    FixedEvent *mRNA_transl_init_time_end_head;   /* times when mRNAs become fully loaded with ribosomes and start producing protein */
    FixedEvent *mRNA_transl_init_time_end_tail;  
    FixedEvent *mRNA_transcr_time_end_head;  /* times when transcription is complete and an mRNA is available to move to cytoplasm */
    FixedEvent *mRNA_transcr_time_end_tail;
    FixedEvent *signal_off_head;          /* times when env=A ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.*/
    FixedEvent *signal_off_tail;   
    FixedEvent *signal_on_head;          /* times when env=A ends. Note, this event is not gene- or copy-specific. I just use the structure of FixedEvent for convenience.*/
    FixedEvent *signal_on_tail; 
    FixedEvent *burn_in_growth_rate_head;
    FixedEvent *burn_in_growth_rate_tail;  
    FixedEvent *sampling_point_end_head;
    FixedEvent *sampling_point_end_tail;
    FixedEvent *change_signal_strength_head;
    FixedEvent *change_signal_strength_tail;

    char effect_of_effector;
    int cell_activated;
    float t_to_update_probability_of_binding;
    float P_A[MAX_GENES];
    float P_R[MAX_GENES];
    float P_A_no_R[MAX_GENES];
    float P_NotA_no_R[MAX_GENES];
    float last_P_R[MAX_GENES];
    float last_P_A[MAX_GENES];
    float last_P_A_no_R[MAX_GENES];
    float last_P_NotA_no_R[MAX_GENES];
    float last_event_t;

    float protein_number[MAX_PROTEINS];     /* pooled protein number from gene_specific_protein_conc */
    float gene_specific_protein_number[MAX_GENES]; /* stores the "protein" number for each gene.
                                               * can be considered temporary data. Make muation easier to
                                               * deal with. */  
    float protein_synthesis_index[MAX_GENES];  /*this is N_mRNA*translation_rate/degradation rate.*/
    int transcriptional_state[MAX_GENES];       /*can be REPRESSED, INTERMEDIATE, or ACTIVE */
    
    int error;
};

typedef struct Mutation Mutation;
struct Mutation
{
    char mut_type;
    int which_gene;
    int which_nucleotide;
    char nuc_diff[3];    /*the first three elements store the nuc after mutation (max_indel=3)*/
    int kinetic_type;    /*0 for pic_disassembly, 1 for mRNA_decay, 2 for translation, 3 for protein_decay*/
    float kinetic_diff;
    int N_hit_bound;
};

typedef struct Phenotype Phenotype;
struct Phenotype
{    
    int replicate_id;
    int timepoint;
    int total_time_points;
    float *protein_concentration;
    float *gene_specific_concentration;
    float *instantaneous_fitness;    
};


/*
 * global variables
 */

/* see netsim.c for documentation for these global constant variables */
extern int MAXELEMENTS;
extern float signal_profile_matrix[N_THREADS][200][15];
/*output files*/
extern char error_file[32];
extern char mutation_file[32];
extern char RuntimeSumm[32];
extern char output_file[32];

/* see netsim.c for documentation for these global variables */
extern int N_EFFECTOR_GENES;
extern int N_TF_GENES;
extern float MAX_ACT_TO_INT_RATE;
extern float MIN_ACT_TO_INT_RATE;
extern float MAX_MRNA_DECAY;
extern float MIN_MRNA_DECAY;
extern float MAX_PROTEIN_DECAY;
extern float MIN_PROTEIN_DECAY;
extern float MAX_PROTEIN_SYN_RATE;
extern float MIN_PROTEIN_SYN_RATE;
extern float MAX_KD;
extern float MIN_KD;
extern int MAX_GENE_LENGTH;
extern int MIN_GENE_LENGTH;

/* function prototypes */

char set_base_pair(float);

void initialize_sequence(char [], int, int, RngStream);

void initialize_genotype(Genotype *, int, int, int, int, RngStream) ;

void initialize_genotype_fixed(Genotype *, int, int, int, RngStream);

void initialize_cell(Genotype *, CellState *, Test *, int [MAX_GENES], float [MAX_PROTEINS]);

void calc_all_binding_sites_copy(Genotype *, int);

void calc_all_binding_sites(Genotype *);

void calc_all_binding_sites2(Genotype *);

void calc_all_binding_sites_copy2(Genotype *, int);

void calc_TF_dist_from_all_BS(  AllTFBindingSites *,
                                int ,
                                int ,
                                int ,                                         
                                int [MAX_PROTEINS],                                    
                                int [3],
                                float [MAX_GENES],
                                int,
                                float *,
                                float *,
                                float *,
                                float *);

int add_fixed_event(int,                           
                   float,
                   FixedEvent **,
                   FixedEvent **);

void delete_fixed_event(int,                               
                               int,
                               FixedEvent **,
                               FixedEvent **);

int does_fixed_event_end(CellState*, float);

int does_fixed_event_end_plotting(CellState *,float);

void update_protein_number_and_fitness(Genotype *,
                                        CellState *,
                                        GillespieRates *,
                                        float);

void fixed_event_end_transcription( float *,
                                    CellState *,
                                    GillespieRates *,
                                    Genotype *);
float calc_tprime(Genotype*, CellState*, float*, float, float, int);

float calc_integral(Genotype *, CellState *, float *, float, float);

float calc_fitness( float *,
                    Genotype *,
                    CellState *,
                    float*,                                       
                    float);

int Gillespie_event_mRNA_decay(GillespieRates *, CellState *, Genotype *, RngStream);

void Gillespie_event_repressed_to_intermediate(GillespieRates *, CellState *, Genotype *, RngStream);

void Gillespie_event_intermediate_to_repressed(GillespieRates *, CellState *, Genotype *, RngStream); 

void Gillespie_event_intermediate_to_active(GillespieRates *, CellState *, Genotype *, RngStream);

void Gillespie_event_active_to_intermediate(Genotype *, CellState *, GillespieRates *, RngStream);

void Gillespie_event_transcription_init(GillespieRates *, CellState *, Genotype *, float, RngStream);

void do_single_timestep(Genotype *, CellState *, GillespieRates *, Test *, Phenotype *, RngStream) ;

void do_single_timestep_plotting(Genotype *, 
                                CellState *,
                                GillespieRates *,  
                                Test *,
                                
                                float (*)[N_TIMP_POINTS],
                                float [N_TIMP_POINTS],                                      
                                int [4],
                                float [4],
                                RngStream,
                                int *) ;
							   
void free_fixedevent(CellState *);
 
void calc_avg_growth_rate(  Genotype *, 
                            Selection *,
                            int [MAX_GENES],
                            float [MAX_PROTEINS],                                  
                            RngStream [N_THREADS],                                                                     
                            int,
                            float *,
                            float *); 
  
int init_run_pop(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_GENES], RngStream, RngStream [N_THREADS]);

void calc_all_rates(Genotype *, CellState *, GillespieRates *, Test *, int);

int fixed_event_end_translation_init(   Genotype *, 
                                        CellState *,    
                                        GillespieRates *, 
                                        float *);

int do_fixed_event( Genotype *, 
                    CellState *, 
                    GillespieRates *, 
                    Test *,
                    Phenotype *,
                    float *,              
                    int);

int do_Gillespie_event(Genotype*, CellState *, GillespieRates *, float, RngStream);

void initialize_cache(Genotype *);

void clone_genotype(Genotype *, Genotype *);

float try_fixation(Genotype *, Genotype *, int, int, int *, RngStream);

void calc_avg_growth_rate_plotting( Genotype *, Selection *, int [MAX_GENES],float [MAX_PROTEINS], RngStream [N_THREADS]); 

void summarize_binding_sites(Genotype *,int);

void summarize_binding_sites2(Genotype *,int);

int check_concurrence(  float , 
                        FixedEvent *, 
                        FixedEvent *, 
                        FixedEvent *, 
                        FixedEvent *,
                        FixedEvent *,
                        float,
                        FixedEvent *);

void set_signal(CellState *, Test *, RngStream, int);

void output_genotype(Genotype *, int);

void release_memory(Genotype *,Genotype *, RngStream *);

void calc_fx_dfx(float, int, float, float*, float*, float*, float*, float*);

int evolve_N_steps( Genotype *, 
                    Genotype *,
                    Mutation *, 
                    Selection *,
                    int *, 
                    int, 
                    int *,  
                    int [MAX_GENES],
                    float [MAX_PROTEINS],
                    RngStream,
                    RngStream [N_THREADS],
                    int );

void run_simulation( Genotype *, 
                    Genotype *, 
                    Mutation *,
                    Selection *,
                    Selection *,                    
                    int [MAX_GENES],   
                    float [MAX_PROTEINS],
                    int,
                    int,
                    RngStream,
                    RngStream [N_THREADS]);

void continue_simulation(Genotype *, 
                        Genotype *, 
                        Mutation *,
                        Selection *,
                        Selection *,
                        int,
                        int [MAX_GENES],
                        float [MAX_PROTEINS],                         
                        RngStream,
                        RngStream [N_THREADS]);

void run_plotting(  Genotype *,
                    Genotype *,
                    Mutation *,
                    Selection *,
                    int [MAX_GENES],
                    float [MAX_GENES],
                    int);

void calc_fitness_stats(Genotype *,
                        Selection *,
                        float (*)[N_REPLICATES],
                        float (*)[N_REPLICATES],                             
                        int);

void evolve_neutrally(  Genotype *,
                        Genotype *,                             
                        Mutation *,                              
                        RngStream);

void replay_mutations(  Genotype *,
                        Genotype *,
                        Mutation *,
                        FILE *,
                        int);

void plot_alternative_fitness(  Genotype *,
                                Genotype *,
                                Mutation *,
                                Selection *,
                                int [MAX_GENES],
                                float [MAX_GENES]);

void find_ffl(Genotype *);

void tidy_output_files(char*, char*);

void print_core_c1ffls(Genotype *);

void remove_edges_iteratively(Genotype *);

void modify_topology(Genotype *, Genotype *);

void add_binding_site(Genotype *, int);

void remove_binding_sites(Genotype *, int);

void calc_leaping_interval(Genotype*, CellState*, float *, float, int);

void print_mutatable_parameters(Genotype*,int);
#endif /* !FILE_NETSIM_SEEN */
