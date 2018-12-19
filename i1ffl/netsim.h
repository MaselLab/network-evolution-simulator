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
#define PHENOTYPE 0
#define PERTURB 1
#define NEUTRAL 0
#define IRREG_SIGNAL 0

/*Runtime control*/  
#define MAX_MUTATIONS 1600000
#define MAX_TRIALS 2000
#define HI_RESOLUTION_RECALC 5
#define N_THREADS 5
#define N_REPLICATES 100
#define OUTPUT_INTERVAL 10
#define SAVING_INTERVAL 10
#define N_TIMEPOINTS 900 // for plotting
#define OUTPUT_MUTANT_DETAILS 0

/*Miscellaneous settings*/
#define EPSILON 1.0e-6       /* original code used EPSILON 10^-6 */
#define OUTPUT_RNG_SEEDS 1
#define TIME_INFINITY 9.99e10
#define TIME_OFFSET 0.01
#define CAUTIOUS 0

/*Biology and evolution settings*/
#define EFFECTOR_NOT_TF 0
#define POOL_EFFECTORS 1
#define FAST_SS 0
#define SELECT_SENSITIVITY_AND_PRECISION 0
#define SELECT_ON_DURATION 0
#define REALLY_COMPLETECATE 0
#define PEAK_SEARCH 1
#define CHEMOTAXIS 0
#define DIRECT_REG 1
#define NO_PENALTY 0
#define IGNORE_BS 1
#define RM_PF 0
#define RM_I1FFL 0
#define RM_NFBL 0
#define RM_PF_NFBL 1
#define MERGE_PROTEIN 0
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
//#ifndef MAX_NON_OUTPUT_GENES             /* number of genes encoding TFs */
#define MIN_NON_OUTPUT_GENES 2         /* the initial value is set in initiate_genotype*/
#define MAX_COPIES_PER_NON_OUTPUT_GENE 4
//#endif
//#ifndef MAX_NON_OUTPUT_PROTEINS
//#define MIN_NON_OUTPUT_PROTEINS 4
//#endif
//#ifndef MAX_OUTPUT_GENES
#define MAX_OUTPUT_GENES 4  /* this is the upper limit of effector gene copies*/
//#endif
//#ifndef MAX_OUTPUT_PROTEINS
#define MAX_OUTPUT_PROTEINS 4
//#endif
//#ifndef MAX_GENES
#define MAX_GENES 26//MIN_NON_OUTPUT_GENES+MAX_OUTPUT_GENES+1+N_SIGNAL_TF  /* total number of genes: add the (non-TF) selection gene to the total (default case) */
//#endif
//#ifndef MAX_PROTEINS           
#define MAX_PROTEINS 26//MIN_NON_OUTPUT_GENES+MAX_OUTPUT_PROTEINS+1+N_SIGNAL_TF
//#endif
#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 8      /* length of binding element on TF */
#define NMIN 6                /* minimal number of nucleotide that matches the binding sequence of a TF in its binding site*/
                              /* DO NOT MAKE NMIN<TF_ELEMENT_LEN/2, OTHERWISE calc_all_binding_sites_copy will make mistake*/  
//#define NUM_K_DISASSEMBLY 131 /* number of differents for PIC disassembly from data file  */
//#ifndef HIND_LENGTH
#define HIND_LENGTH 3         /* default length of hindrance on each side of the binding site,i.e. a tf occupies TF_ELEMENT_LEN+2*HIND_LENGTH */
                              /* the binding of Lac repressor blockes 12 bp. Record MT 1981*/
//#endif
#define MAX_BINDING 10  /* MAX_MODE is the max number of tf that can bind to a promoter plus 1*/

/* 
 * define macros for logging warning/errors 
 */
#define MAKE_LOG 0
#if MAKE_LOG
  #define LOG(...) { fprintf(fperrors, "%s: ", __func__); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); } 
#endif

/*
 * primary data structures for model
 */
enum PROTEIN_IDENTITY {ACTIVATOR=1, REPRESSOR=0, NON_OUTPUT_PROTEIN=-1,OUTPUT_PROTEIN=-2,NON_TF=-3};
enum BOOLEAN {NA=-1, NO=0, YES=1};

typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         /* transcription factor */
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
typedef struct Environment Environment;
struct Environment
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
//    float benefit;
    float response_amplification;
    float minimal_peak_response;
    float opt_pulse_amplitude_fold;
//    float max_t_in_bias;    
    float opt_pulse_duration;
    float fitness_factor_of_duration;
    float fitness_factor_of_amplitude;
    float opt_ss_response;
    float fitness_factor_of_ss_response;
    float fitness_factor_of_diff_in_ss;
    float opt_peak_response;
    float fitness_factor_of_peak_response;    
    int window_size;
};

struct Selection
{
    Environment env1;
    Environment env2;
    float env1_weight;
    float env2_weight;
    float temporary_miu_ACT_TO_INT_RATE;
    float temporary_miu_protein_syn_rate;
    float temporary_miu_Kd;
    float temporary_DUPLICATION;
    float temporary_SILENCING;
    int temporary_N_effector_genes;
    int temporary_N_tf_genes;
    int MAX_STEPS;
};


typedef struct Genotype Genotype;
struct Genotype {
    int ngenes;                                             /* the number of loci */
    int nproteins;                                          /* because of gene duplication, the number of proteins and mRNAs can be different 
                                                               from the number of loci. nprotein is the number of elements in protein_pool*/
    int n_output_genes;    
    int N_node_families;
    int N_cisreg_clusters;
    int which_protein[MAX_GENES];                              /* in case of gene duplication, this array tells the protein corresponding to a given gene id */ 
    
    char cisreg_seq[MAX_GENES][CISREG_LEN];    
    
    /*these apply to protein, not loci*/
    int N_act;                                              /* number of activators*/
    int N_rep;                                              /* number of repressors*/    
    int protein_identity[MAX_PROTEINS];                  /* entry 1 marks activator (1) or repressor (0)*/        
    char tf_binding_seq[MAX_PROTEINS][TF_ELEMENT_LEN];
    char tf_binding_seq_rc[MAX_PROTEINS][TF_ELEMENT_LEN];       /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
    int protein_pool[MAX_PROTEINS][2][MAX_GENES];               /* element 1 record how many genes/mRNAs producing this protein,ele 2 stores which genes/mRNAs*/
    float Kd[MAX_PROTEINS];
    
    /*group genes into node family*/
    int which_node_family[MAX_GENES];                          /* motifs must be made of genes from different node family.
                                                                * In a node family, genes
                                                                * are either all activators or all repressors
                                                                * are either all output protein or non-output protein
                                                                * have identical TFBS (but not necessarily same Kd)  
                                                                */
    int node_family_pool[MAX_GENES][2][MAX_GENES];                            
    
    
    /*these apply to loci*/
    int output_protein_ids[MAX_OUTPUT_GENES];
    int is_output[MAX_GENES];                                   /*is an output protein (id in output_protein_id) or non-output protein (-1)   */
    int locus_length[MAX_GENES];
    int total_loci_length;
    float mRNA_decay_rate[MAX_GENES];                                /* kinetic rates*/
    float protein_decay_rate[MAX_GENES];                             /* kinetic rates*/
    float protein_syn_rate[MAX_GENES];                              /* kinetic rates*/   
    float active_to_intermediate_rate[MAX_GENES];                          /* kinetic rates*/    
    int min_N_activator_to_transc[MAX_GENES];                          /* 1 for OR GATE, at leat 2 FOR AND GATE */ 
    int locus_specific_TF_behavior[MAX_GENES][MAX_PROTEINS];      /* whether a TF behaves as activator or repressor depends on locus*/
    
    /* binding sites related data, applying to loci*/   
    int cisreg_cluster_pool[MAX_GENES][2][MAX_GENES];                     /* For genes having the same cis-reg, tf distribution can be shared.
                                                             * Genes having the same cis-reg are clustered.
                                                               1st dim stores cluster ids, 2nd dim stores gene_ids in a cluster.
                                                               cisreg_cluster works with which_cluster*/
    int which_cluster[MAX_GENES];                              /* which_cluster stores the cluster id of a gene*/                                                           
    int recalc_TFBS[MAX_GENES];                                /* whether to recalc the TFBS*/       
    int binding_sites_num[MAX_GENES];                          /* total number of binding sites */ 
    int N_allocated_elements;
    int max_unhindered_sites[MAX_GENES][3];                    /* maximal number of binding sites that do not hinder each other. 0 for activator BS, 1 for repressor BS*/  
    int max_hindered_sites[MAX_GENES];                        /* maximal number of BSs a BS can hinder*/ 
    int N_act_BS[MAX_GENES];                                   /* total number of binding sites of activating TF */
    int N_rep_BS[MAX_GENES];                                   /* total number of binding sites of repressing TF */
    AllTFBindingSites *all_binding_sites[MAX_GENES];   

    /*statistics of fitness*/
    float fitness1;
    float fitness2;
    float sq_SE_fitness1;
    float sq_SE_fitness2;
    float avg_fitness;
    float sq_SE_avg_fitness;
    float fitness_measurement[MAX_RECALC_FITNESS*N_REPLICATES];
    
    /*measurement of network topology*/
    int N_motifs[17]; 
    int TF_in_core_C1ffl[MAX_GENES][MAX_PROTEINS];
    int gene_in_core_C1ffl[MAX_GENES];
    int N_act_genes;  
    int N_act_effector_genes;
};

typedef struct Phenotype Phenotype;
struct Phenotype
{ 
    int timepoint;
    int total_time_points;
    float *protein_concentration;
    float *gene_specific_concentration;
    float *instantaneous_fitness;    
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


/*output files. Initialized in main.c*/
extern char mutation_file[32];
extern char setup_summary[32];
extern char evo_summary[32];
extern float signal_profile_matrix[N_THREADS][200][15];

/*Initialized in netsim.c */
extern int MAXELEMENTS;
extern const float MEAN_GENE_LENGTH;
extern int N_EFFECTOR_GENES;
extern int N_TF_GENES;
extern const float MAX_ACT_TO_INT_RATE;
extern const float MIN_ACT_TO_INT_RATE;
extern const float MAX_MRNA_DECAY;
extern const float MIN_MRNA_DECAY;
extern const float MAX_PROTEIN_DECAY;
extern const float MIN_PROTEIN_DECAY;
extern const float MAX_PROTEIN_SYN_RATE;
extern const float MIN_PROTEIN_SYN_RATE;
extern const float MAX_KD;
extern const float MIN_KD;
extern const int MAX_GENE_LENGTH;
extern const int MIN_GENE_LENGTH;
extern const float KD2APP_KD;
extern const float sampling_interval;

/* function prototypes */
char set_base_pair(float);

void initialize_genotype(Genotype *, int, int, int, int, RngStream) ;

void calc_all_binding_sites_copy(Genotype *, int);

void calc_all_binding_sites(Genotype *); 

void initialize_cache(Genotype *);

int evolve_under_selection(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_GENES], RngStream, RngStream [N_THREADS]);

void show_phenotype( Genotype *,
                    Genotype *,
                    Mutation *,
                    Selection *,
                    int [MAX_GENES],
                    float [MAX_GENES],                  
                    RngStream [N_THREADS]);

void evolve_neutrally(  Genotype *,
                        Genotype *,                             
                        Mutation *,
                        Selection *,
                        Selection *,
                        RngStream);

void modify_network(Genotype *,
                    Genotype *,
                    Mutation *,
                    Selection *,
                    int [MAX_GENES],
                    float [MAX_GENES],
                    RngStream [N_THREADS]);

void calc_genotype_fitness( Genotype *,
                            Selection *,
                            float (*)[N_REPLICATES],
                            float (*)[N_REPLICATES],
                            int);

void print_mutatable_parameters(Genotype*,int);

#endif /* !FILE_NETSIM_SEEN */
