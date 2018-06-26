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
#define JUST_PLOTTING 0
#define PLOT_ALTERNATIVE_FITNESS 0
#define NEUTRAL 0
#define RUN_FULL_SIMULATION 1
#define SKIP_INITIAL_GENOTYPE 0
#define QUICK_BURN_IN 0
#define EXTERNAL_SIGNAL 0

/*Runtime control*/  
#define MAX_MUT_STEP 50000
#define BURN_IN_I 0
#define BURN_IN_II 0
#define MAX_MUTATIONS 1600000
#define MAX_TRIALS 2000
#define N_THREADS 5
#define N_REPLICATES 100
#define OUTPUT_INTERVAL 10
#define N_TIMEPOINTS 900 // for plotting
#define OUTPUT_MUTANT_DETAILS 0

/*Miscellaneous settings*/
#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define EPSILON 1.0e-6       /* original code used EPSILON 10^-6 */
#define RT_SAFE_EPSILON 1e-6
#define OUTPUT_RNG_SEEDS 1
#define TIME_INFINITY 9.99e10
#define TIME_OFFSET 0.01
#define CAUTIOUS 0

/*Biology and evolution settings*/
#define POOL_EFFECTORS 1
#define SELECT_SENSITIVITY_AND_PRECISION 0
#define SELECT_ON_DURATION 0
#define REALLY_COMPLETECATE 0
#define PEAK_SEARCH 1
#define CHEMOTAXIS 0
#define DIRECT_REG 1
#define NO_PENALTY 0
#define ADD_2_PATHWAYS 0
#define FORCE_OR_GATE 0
#define RANDOM_COOPERATION_LOGIC 0
#define N_SIGNAL_TF 1 // the 1st TF enables basal activity in TFN. The 2nd is the actual signal TF. 
#define NO_REGULATION_COST 0
#define NO_REGULATION 0 // this locks the state of transcription factors to NUC_NO_PIC
#define REGRESSIVE_MUTATION 1
#define IGNORE_BS_OVERLAPPING 0
#define SIMPLE_SUBSTITUTION 1
#define RANDOMIZE_SIGNAL2 0
#define MAX_RECALC_FITNESS 10
#define minimal_selection_coefficient 1.0e-8
/* Because mutation can change the number of genes, the numbers defined here are used to allocate storage space only.
 * Set the numbers to be 8 folds of the initial ngenes and ntfgenes, so that we can have two whole genome duplications*/
#ifndef MAX_NON_OUTPUT_GENES             /* number of genes encoding TFs */
#define MAX_NON_OUTPUT_GENES 21         /* the initial value is set in initiate_genotype*/
#endif
#ifndef MAX_NON_OUTPUT_PROTEINS
#define MAX_NON_OUTPUT_PROTEINS 21
#endif
#ifndef MAX_OUTPUT_GENES
#define MAX_OUTPUT_GENES 10  /* this is the upper limit of effector gene copies*/
#endif
#ifndef MAX_OUTPUT_PROTEINS
#define MAX_OUTPUT_PROTEINS 10
#endif
#ifndef NGENES
#define NGENES MAX_NON_OUTPUT_GENES+MAX_OUTPUT_GENES+1  /* total number of genes: add the (non-TF) selection gene to the total (default case) */
#endif
#ifndef NPROTEINS           
#define NPROTEINS MAX_NON_OUTPUT_PROTEINS+MAX_OUTPUT_PROTEINS+1
#endif
#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 8      /* length of binding element on TF */
#define NMIN 6                /* minimal number of nucleotide that matches the binding sequence of a TF in its binding site*/
                              /* DO NOT MAKE NMIN<TF_ELEMENT_LEN/2, OTHERWISE calc_all_binding_sites_copy will make mistake*/  
//#define NUM_K_DISASSEMBLY 131 /* number of differents for PIC disassembly from data file  */
#ifndef HIND_LENGTH
#define HIND_LENGTH 3         /* default length of hindrance on each side of the binding site,i.e. a tf occupies TF_ELEMENT_LEN+2*HIND_LENGTH */
                              /* the binding of Lac repressor blockes 12 bp. Record MT 1981*/
#endif
#define MAX_BINDING 10  /* MAX_MODE is the max number of tf that can bind to a promoter plus 1*/

/* 
 * define macros for logging warning/errors 
 */
#ifndef LOGGING_OFF
  #define LOG(...) { fprintf(fperrors, "%s: ", __func__); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); } 
#else
  #define LOG 
#endif

/*
 * primary data structures for model
 */
enum GENE_STATE {REPRESSED, INTERMEDIATE, ACTIVE};
enum PROTEIN_IDENTITY {ACTIVATOR=1, REPRESSOR=0, NON_OUTPUT_PROTEIN=-1,NON_TF=-2};
enum BOOLEAN {NA=-1, NO=0, YES=1};


/*
 * Rates for Gillespie algorithm
 */
typedef struct GillespieRates GillespieRates;
struct GillespieRates {
  float total_mRNA_decay_rate;         
  float mRNA_decay_rate[NGENES];
  float total_active_to_intermediate_rate;   
  float active_to_intermediate_rate[NGENES]; 
  float total_repressed_to_intermediate_rate;
  float repressed_to_intermediate_rate[NGENES]; 
  float total_intermediate_to_repressed_rate;
  float intermediate_to_repressed_rate[NGENES];  
  float total_intermediate_to_active_rate;
  float intermediate_to_active_rate[NGENES];  
  int total_N_gene_transcript_initiated; 
  int transcript_initiation_state[NGENES];   
  float total_Gillespie_rate;
};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         /* transcription factor */
  float Kd;
  int mis_match;     /* number of mismatched nuc */
  int BS_pos;        /* start position of BS on DNA, always with reference to forward strand */                     
  int N_hindered;    /* the number of BS hindered by this TF when it binds to the current BS */  
};

typedef struct Genotype Genotype;
struct Genotype {
    int ngenes;                                             /* the number of loci */
    int nproteins;                                          /* because of gene duplication, the number of proteins and mRNAs can be different 
                                                               from the number of loci. nprotein is the number of elements in protein_pool*/
    int n_output_genes;
    int n_output_proteins;
    int nTF_families;
    int which_protein[NGENES];                              /* in case of gene duplication, this array tells the protein corresponding to a given gene id */ 
    
    char cisreg_seq[NGENES][CISREG_LEN];    
    
    /*these apply to protein, not loci*/
    int N_act;                                              /* number of activators*/
    int N_rep;                                              /* number of repressors*/    
    int protein_identity[NPROTEINS][2];                     /* entry 1 marks activator (1) or repressor (0); */
                                                            /* entry 2 marks output protein (id in output_protein_id) or
                                                             * non-output protein (-1) */
    int output_protein_id[MAX_OUTPUT_PROTEINS];
    char tf_binding_seq[NPROTEINS][TF_ELEMENT_LEN];
    char tf_binding_seq_rc[NPROTEINS][TF_ELEMENT_LEN];                /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
    int protein_pool[NPROTEINS][2][NGENES];                 /* element 1 record how many genes/mRNAs producing this protein,ele 2 stores which genes/mRNAs*/
    int which_TF_family[NPROTEINS];                         /* We consider mutation to Kd does not create a new TF (matters when counting motifs), 
                                                             * but the calculation of TF binding distribution demands tfs with different Kd
                                                             * to be treated differently. We call tfs that differ only in Kd a tf family.
                                                             * We use this array and TF_family_pool to track which TF belongs which TF family.
                                                             */
    int TF_family_pool[NPROTEINS][2][NPROTEINS];                            
    float Kd[NPROTEINS];
    
    /*these apply to loci*/
    int locus_length[NGENES];
    int total_loci_length;
    float mRNA_decay_rate[NGENES];                                /* kinetic rates*/
    float protein_decay_rate[NGENES];                             /* kinetic rates*/
    float protein_syn_rate[NGENES];                              /* kinetic rates*/   
    float active_to_intermediate_rate[NGENES];                          /* kinetic rates*/    
    int min_N_activator_to_transc[NGENES];                          /* 1 for OR GATE, at leat 2 FOR AND GATE */ 
    int locus_specific_TF_behavior[NGENES][NPROTEINS];      /* whether a TF behaves as activator or repressor depends on locus*/
    
    /* binding sites related data, applying to loci*/   
    int cisreg_cluster[NGENES+1][NGENES];                     /* For genes having the same cis-reg, tf distribution can be shared.
                                                             * Genes having the same cis-reg are clustered.
                                                               1st dim stores cluster ids, 2nd dim stores gene_ids in a cluster.
                                                               cisreg_cluster works with which_cluster*/
    int which_cluster[NGENES];                              /* which_cluster stores the cluster id of a gene*/                                                           
    int recalc_TFBS[NGENES];                                /* whether to recalc the TFBS*/       
    int binding_sites_num[NGENES];                          /* total number of binding sites */ 
    int N_allocated_elements;
    int max_unhindered_sites[NGENES][3];                    /* maximal number of binding sites that do not hinder each other. 0 for activator BS, 1 for repressor BS*/  
    int max_hindered_sites[NGENES];                        /* maximal number of BSs a BS can hinder*/ 
    int N_act_BS[NGENES];                                   /* total number of binding sites of activating TF */
    int N_rep_BS[NGENES];                                   /* total number of binding sites of repressing TF */
    AllTFBindingSites *all_binding_sites[NGENES];   

    /*statistics of fitness*/
    float avg_GR1;
    float avg_GR2;
    float sq_SE_GR1;
    float sq_SE_GR2;
    float fitness;
    float sq_SE_fitness;
    float fitness_measurement[MAX_RECALC_FITNESS*N_REPLICATES];
    
    /*measurement of network topology*/
    int N_motifs[11]; 
    int TF_in_core_C1ffl[NGENES][NPROTEINS];
    int gene_in_core_C1ffl[NGENES];
    int N_act_genes;  
    int N_act_effector_genes;
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
    float cumulative_fitness[MAX_OUTPUT_PROTEINS];                    
    float cumulative_fitness_after_burn_in[MAX_OUTPUT_PROTEINS];          
    float instantaneous_fitness[MAX_OUTPUT_PROTEINS];  
    float cumulative_cost;
    int mRNA_aft_transl_delay_num[NGENES];          /* mRNAs that have finished the translational delay */
    int mRNA_under_transl_delay_num[NGENES];   /* mRNAs that are still under the translational delay (they do not contribute to protein 
                                                * turnover, but contribute to the cost of translation)  */
    int mRNA_under_transc_num[NGENES];       /* mRNAs which haven't finished transcription */
    float protein_number[NPROTEINS];     /* pooled protein number from gene_specific_protein_conc */
    float gene_specific_protein_number[NGENES]; /* stores the "protein" number for each gene.
                                               * can be considered temporary data. Make muation easier to
                                               * deal with. */
    int gene_state[NGENES];       /* gives the state of each of the genes, according to figure
                                          *  NUC_NO_PIC = 0,
                                          *  NO_NUC_NO_PIC =1,
                                          *  PIC_NO_NUC = 3,
                                          */ 
    float protein_synthesis_index[NGENES];  /*this is N_mRNA*translation_rate/degradation rate.*/
    int cell_activated;
    float t_to_update_probability_of_binding;
    float P_A[NGENES];
    float P_R[NGENES];
    float P_A_no_R[NGENES];
    float P_NotA_no_R[NGENES];
    float last_P_R[NGENES];
    float last_P_A[NGENES];
    float last_P_A_no_R[NGENES];
    float last_P_NotA_no_R[NGENES];
    float last_event_t;  

    /*linked table to the timing of fixed events*/
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


    /*measurement of pulse*/   
    float Duration_pulse;
    float T_pulse_on;
    float T_pulse_off;
    int first_pulse;
    int Pulse_is_on; 
    float *sampled_response;
    int N_samples;
    float cumulative_basal_benefit;
    float cumulative_advanced_benefit;
    float cumulative_damage;   
    /*chemotaxis*/
    float cumulative_benefit;    
    int recording_basal_response;
    int found_gradient;
    int moving;
    int position;                       /*0: in poor media, 1: in rich media, 2: rich to poor*/
    float cumulative_t_in_bias;
    float threshold_response_to_bias[MAX_OUTPUT_PROTEINS];
    float basal_response[MAX_OUTPUT_PROTEINS];
//    float response_amplification;
    /***/
    float sensitivity[3]; /* element 1 is a running record, 2 is the sensitivity for singal change 1, 3 is the sensitivity for signal change 2*/
    float precision[3];  
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
    int which_gene;
    int which_protein;
    int which_nucleotide;
    char nuc_diff[3];    /*the first three elements store the nuc after mutation (max_indel=3)*/
    int kinetic_type;    /*0 for pic_disassembly, 1 for mRNA_decay, 2 for translation, 3 for protein_decay*/
    float kinetic_diff;
    int N_hit_bound;
};

/* global variables: output files*/
char error_file[32];
char mutation_file[32];
char RuntimeSumm[32];
char output_file[32];
FILE *fperrors;

/*mutation*/
float SUBSTITUTION; 
float DUPLICATION;   
float SILENCING;      
float MUT_Kd;
float MUT_ACT_to_INT;
float MUT_mRNA_decay;
float MUT_protein_decay;
float MUT_protein_syn_rate;
float MUT_identity;
float MUT_binding_seq;
float MUT_LOCUS_LENGTH_RATE;
float MUT_cooperation;
float MUT_effector_to_TF;
float MUT_locus_specific_tf_behavior;
float mutational_regression_rate;
float sigma_ACT_TO_INT_RATE; 
float sigma_mRNA_decay; 
float sigma_protein_decay; 
float sigma_protein_syn_rate; 
float sigma_Kd;
float miu_ACT_TO_INT_RATE;
float miu_mRNA_decay;
float miu_protein_decay;
float miu_protein_syn_rate;
float miu_Kd;
const float MAX_ACT_TO_INT_RATE;
const float MIN_ACT_TO_INT_RATE;
const float MAX_MRNA_DECAY;
const float MIN_MRNA_DECAY;
const float MAX_PROTEIN_DECAY;
const float MIN_PROTEIN_DECAY;
const float MAX_PROTEIN_SYN_RATE;
const float MIN_PROTEIN_SYN_RATE;
const float MAX_KD;
const float MIN_KD;
const float NS_Kd;
const int MAX_GENE_LENGTH;
const int MIN_GENE_LENGTH;


/* function prototypes */
char set_base_pair(float);

extern void initialize_sequence(char [], int, int, RngStream);

extern void print_genotype(Genotype *, int);

extern void print_all_binding_sites(int [NGENES],
                                    AllTFBindingSites *, 
                                    int ,
                                    char [NGENES][TF_ELEMENT_LEN],
                                    char [NGENES][CISREG_LEN]
//                                    int [NGENES][2]
									);

extern void initialize_genotype(Genotype *, RngStream) ;

extern void initialize_genotype_fixed(Genotype *, RngStream);

extern void calc_all_binding_sites_copy(Genotype *, int);

extern void calc_all_binding_sites(Genotype *);

extern void calc_TF_dist_from_all_BS(  AllTFBindingSites *,
                                        int ,
                                        int ,
                                        int ,                                         
                                        int [NPROTEINS],                                    
                                        int [3],
                                        float [NGENES],
                                        int,
                                        float *,
                                        float *,
                                        float *,
                                        float *);

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

extern void initialize_cell(Genotype *,
                            CellState *,
                            int,
                            int,  
                            int [NPROTEINS][2][NGENES],
                            float [NGENES],
                            int [NGENES],
                            float [NPROTEINS],
                            RngStream,float,float);

extern int does_fixed_event_end(CellState*,
                                float);

extern int does_fixed_event_end_plotting(CellState *,float);

extern void update_protein_number_and_fitness(Genotype *,
                                            CellState *,
                                            GillespieRates *,
                                            float, 
//                                            float,
//                                            char,
                                            int *,                                           
                                            int,
                                            Mutation *,
                                            float,
                                            float,
                                            float);

extern void fixed_event_end_transcription(  float *,
                                            CellState *,
                                            GillespieRates *,
                                            Genotype *,
                                            int *,                              
                                            int,
                                            Mutation *,
                                            float, float, float);
extern float calc_tprime(Genotype *, CellState *, float[NGENES], float, float,int[NGENES],int,int);

extern float calc_integral(Genotype *, CellState *, float[MAX_OUTPUT_GENES] , float, float, int[MAX_OUTPUT_GENES],int,int);

extern void calc_instantaneous_fitness(   Genotype *,
                                        CellState *,
                                        float[MAX_OUTPUT_GENES],  
                                        float[MAX_OUTPUT_GENES],
                                        float,
//                                        float,
//                                        char,
                                        int *,                                        
                                        int,
                                        Mutation *,
                                        float,
                                        float,
                                        float);

extern int Gillespie_event_mRNA_decay(GillespieRates *, CellState *, Genotype *, RngStream);

extern void Gillespie_event_repressed_to_intermediate(GillespieRates *, CellState *, Genotype *, RngStream);

extern void Gillespie_event_intermediate_to_repressed(GillespieRates *, CellState *, Genotype *, RngStream); 

extern void Gillespie_event_intermediate_to_active(GillespieRates *, CellState *, Genotype *, RngStream);

extern void Gillespie_event_active_to_intermediate(Genotype *, CellState *, GillespieRates *, RngStream);

extern void Gillespie_event_transcription_init(GillespieRates *, CellState *, Genotype *, float, RngStream);

extern void do_single_timestep( Genotype *, 
                                CellState *,
                                GillespieRates *,  
                                float,
                                float, 
                                float,
                                float,
                                float,
                                float,
                                float,
                                float,        
                                RngStream,                             
                                int,
                                Mutation *,
                                int *,                          
                                int,
                                float,
                                float *) ;

extern void do_single_timestep_plotting(    Genotype *, 
                                            CellState *,
                                            GillespieRates *,                                
                                            float,
                                            float, 
                                            float,
                                            float,
                                            float,
                                            float,
                                            float,
                                            float,  
                                            float (*)[N_TIMEPOINTS],
                                            float [N_TIMEPOINTS],
                                            RngStream,                              
                                            int *,
                                            int *,
                                            int,
                                            float,
                                            float *) ;
							   
extern void free_fixedevent(CellState *);
 
extern void calc_cellular_fitness(   Genotype *,                                    
                                    int [NGENES],
                                    float [NPROTEINS],                                  
                                    RngStream [N_THREADS],                                                                     
                                    int,
                                    float *,
                                    float *,
                                    Mutation *); 
  
extern int init_run_pop(unsigned long int [6], int);

extern void calc_all_rates(Genotype *,
                            CellState *,
                            GillespieRates *,  
                            float,
                            int,
                            int);

extern int fixed_event_end_translation_init(   Genotype *, 
                                                CellState *,    
                                                GillespieRates *, 
                                                float *,   
//						float,
                                                int *,                                    
                                                int,
                                                Mutation *,
//                                                char *,
                                                float, float, float);

extern int do_fixed_event(  Genotype *, 
                            CellState *, 
                            GillespieRates *, 
                            float *,  
//                            float,
                            int , 
                            float,
                            float,
                            float,  
                            float,
                            float,
                            float,
                            float,
                            float,
//                            char *,
                            int *,                          
                            int,
                            Mutation *,
                            float *);

extern int do_fixed_event_plotting( Genotype *, 
                                    CellState *, 
                                    GillespieRates *,                                     
                                    float *,                                
                                    int , 
                                    float,
                                    float,
                                    float,  
                                    float,
                                    float,
                                    float,
                                    float, 
                                    float,
                                    int *,
                                    float *);

extern int do_Gillespie_event(Genotype*, CellState *, GillespieRates *, float, RngStream, int *, int, Mutation *);

extern void initialize_cache(Genotype *);

extern void clone_genotype(Genotype *, Genotype *);

extern float try_fixation(Genotype *, Genotype *, int, int, int *, RngStream);

extern void calc_cellular_fitness_plotting(  Genotype *,                               
                                            int [NGENES],
                                            float [NPROTEINS], 
                                            RngStream [N_THREADS]); 

extern void summarize_binding_sites(Genotype *,int);

extern int check_concurrence(   float , 
                                FixedEvent *, 
                                FixedEvent *, 
                                FixedEvent *, 
                                FixedEvent *,
                                FixedEvent *,
                                float,
                                FixedEvent *);

extern void set_signal(CellState *, float, float, float *, float,float);

extern void output_genotype(Genotype *, int);

extern void release_memory(Genotype *,Genotype *, RngStream *, RngStream [N_THREADS]);

extern void calc_fx_dfx(float, int, float, float*, float*, float*, float*, float*, int);

extern int evolve_N_steps(  Genotype *, 
                            Genotype *,
                            int *, 
                            int, 
                            int *,   
                            float [NPROTEINS],
                            int [NGENES],
                            Mutation *, 
                            RngStream,
                            RngStream [N_THREADS],
                            int );

extern void run_simulation( Genotype *, 
                            Genotype *,                     
                            float [NPROTEINS],
                            int [NGENES],
                            Mutation *,
                            int,
                            int,
                            RngStream,
                            RngStream [N_THREADS]);

extern void continue_simulation(Genotype *, 
                                Genotype *,                                
                                int, 
                                float [NPROTEINS],
                                int [NGENES],
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

extern void calc_genotype_fitness( Genotype *,
                                    float (*)[N_REPLICATES],
                                    float (*)[N_REPLICATES],
                                    int);


extern void evolve_neutrally(   Genotype *,
                                Genotype *,                             
                                Mutation *,
                                int,
                                RngStream);


extern void replay_mutations(   Genotype *,
                                Genotype *,
                                FILE *,
                                Mutation *,
                                int,
                                RngStream);

extern void plot_alternative_fitness(   Genotype *,
                                        Genotype *,
                                        int [NGENES],
                                        float [NGENES],
                                        RngStream [N_THREADS],
                                        Mutation *,
                                        FILE *,
                                        int);

extern void find_i1ffl(Genotype *);

extern void tidy_output_files(char*, char*);

extern void print_core_i1ffls(Genotype *);

extern float calc_replicate_fitness(CellState *, int, float, float, float, float, int);

extern void calc_leaping_interval(Genotype*, CellState*, float *, float, int);

#endif /* !FILE_NETSIM_SEEN */