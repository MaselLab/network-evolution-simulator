/*
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2018 Arizona Board of Regents on behalf of the University of Arizona
 
 * This file is part of network-evolution-simulator.
 * network-evolution-simulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * network-evolution-simulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * You should have received a copy of the GNU Affero General Public License
 * along with network-evolution-simulator. If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>
#include "RngStream.h"

/*Simulation mode*/
#define PHENOTYPE 0
#define PERTURB 0
#define NEUTRAL 0 //useless

/*Runtime control*/  
#define MAX_MUTATIONS 1600000 // try at most 1600000 mutants during evolution before shut down anyway 
#define MAX_TRIALS 2000 //try at most 2000 mutants at an evolutionary step before shut down anyway 
#define N_REPLICATES 200 //calculate the fitness of a mutant with 200 replicates
#define HI_RESOLUTION_RECALC 5 //calcualte the high-resolution fitness of a resident with 5*N_REPLICATES replicates
#define LOW_RESOLUTION_RECALC 1 //calcualte the low-resolution fitness of a resident with N_REPLICATES replicates
#define N_THREADS 10 //the number of parallel OpenMP threads
#define OUTPUT_INTERVAL 50 //pool results from evolutionary steps before writing to disk
#define OUTPUT_MUTANT_DETAILS 1 //output every mutant genotype and its fitness, whetehr the mutant is accepted. Don't change it if need to characterize mutations

/*simulation mode under PHENOTYPE*/
#if PHENOTYPE
#define SAMPLE_GENE_EXPRESSION 0 //output expression timecourse of all genes at particular evolutionary step 
#define SAMPLE_EFFECTOR_EXPRESSION 0 //output expression timecourse of the effector during the first N evolutionary steps 
#endif

/*simulation mode under PHENOTYPE*/
#if PERTURB
#define CLEAN_UP_NETWORK 0 //identify and exclude non-adaptive 2-mismatch TFBSs 
#define CLASSIFY_MUTATION 0 //classify all mutations into motif-creating and motif-destroying
#endif

/*parameters used by PERTURB to call 2-mismatch TFBSs*/
#define CUT_OFF_NONADAPTIVE_TFBS -0.01 // fitness cutoff used to call non-adaptive 2-mismatch TFBSs when fitness is calculated at high resolution 
#define CUT_OFF_NONADAPTIVE_TFBS_2 -0.02 // fitness cutoff used to call non-adaptive 2-mismatch TFBSs when fitness is calculated at low resolution 

/*or manually remove 2-mismatch TFBSs in networks.txt and N_motifs.txt*/
#define MAX_MISMATCH_SIGNAL2GENES  2 // set to 1 to remove 2-mismatch TFBSs of the signal in all genes
#define MAX_MISMATCH_EFFECTOR2GENES  2 // set to 1 to remove 2-mismatch TFBSs of the effector in all genes
#define MAX_MISMATCH_NONEFFECTOR2GENES  2 // set to 1 to remove 2-mismatch TFBSs of the non-effector (excluding the signal) in all genes

/* Because mutation can change the number of genes, the numbers defined here are used to allocate storage space only.
 * Set the numbers to be 8 folds of the initial ngenes and ntfgenes, so that we can have two whole genome duplications*/
#define MIN_NON_OUTPUT_GENES 2         // the initial value is set in initiate_genotype
#define MAX_COPIES_PER_NON_OUTPUT_GENE 4 // this is the upper limit of non-effector gene copies
#define MAX_OUTPUT_GENES 4  // this is the upper limit of effector gene copies
#define MAX_OUTPUT_PROTEINS 4
#define MAX_GENES 26 // This is the total number of genes = MIN_NON_OUTPUT_GENES+MAX_OUTPUT_GENES+1+N_SIGNAL_TF.        
#define MAX_PROTEINS 26

/* parameters that control TFBSs */
#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 8      /* length of binding element on TF */
#define NMIN 6                /* minimal number of nucleotide that matches the binding sequence of a TF in its binding site*/
                              /* DO NOT MAKE NMIN<TF_ELEMENT_LEN/2, OTHERWISE calc_all_binding_sites_copy will make mistake*/ 
#define HIND_LENGTH 3         /* default length of hindrance on each side of the binding site,i.e. a tf occupies TF_ELEMENT_LEN+2*HIND_LENGTH */
                              /* the binding of Lac repressor blockes 12 bp. Record MT 1981*/
#define MAX_BINDING 10  /* MAX_MODE is the max number of tf that can bind to a promoter plus 1*/

/* 
 * define macros for logging warning/errors 
 */
#define MAKE_LOG 0
#if MAKE_LOG
  #define LOG(...) { fprintf(fperrors, "%s: ", __func__); fprintf (fperrors, __VA_ARGS__) ; fflush(fperrors); } 
#endif

/*Miscellaneous settings*/
#define EPSILON 1.0e-6       
#define OUTPUT_RNG_SEEDS 1
#define TIME_INFINITY 9.99e10
#define TIME_OFFSET 0.01
#define CAUTIOUS 0

/*don't change these parameters*/
#define EVOLVE_I1FFL 1
#define DIRECT_REG 1
#define NO_PENALTY 0
#define N_SIGNAL_TF 1 // the 1st TF enables basal activity in TFN. The 2nd is the actual signal TF. 
#define POOL_EFFECTORS 1 //useless
#define N_OUTPUT 1 //useless
#define MERGE_PROTEIN 0  //useless
#define FORCE_OR_GATE 0 //useless
#define IGNORE_BS 1 //useless

/*
 * primary data structures for model
 */
enum PROTEIN_IDENTITY {ACTIVATOR=1, REPRESSOR=0, NON_OUTPUT_PROTEIN=-1,OUTPUT_PROTEIN=-2, NON_TF=-3};
enum BOOLEAN {NA=-1, NO=0, YES=1};
enum PERTURBATION {RM_NONE=0, RM_T2E=50, RM_E2T=100, RM_S2T=200, RM_ES2T=300, RM_E2E=400, RM_ES2ET=600, RM_EE2ET=500, RM_EES2ETT=700, RM_TE2ET=150, RM_TS2ET=250, RM_TE2E=450, RM_EST2TTE=350, RM_EET2TEE=550, RM_SET2TEE=650, RM_ESET2TTEE=750};

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
    float signal_strength_stage1;
    float signal_strength_stage2;
    float t_stage1;
    float t_stage2;
    int signal_on_aft_burn_in;
    char initial_effect_of_effector;
    char effect_of_effector_aft_burn_in;
    int fixed_effector_effect;  
    float *external_signal;
    float max_duration_of_burn_in_growth_rate;    
    float avg_duration_of_burn_in_growth_rate;
    int is_burn_in;
    
/*selection condition for I1-FFL*/
    float max_t_mid; 
    float opt_peak_response;
    float min_reduction_relative_to_peak;
    float effector_level_before_stage2; 
    float fitness_decay_constant;
    int window_size; // size of windows used to average effector levels
    float w1; //weight of fitness component
    float w2;
    float w3;
    float w4;
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
    int effector_is_TF;
    int signal_ctrl_repressor;
    int aux_is_TF;    
};


typedef struct Genotype Genotype;
struct Genotype {
    int flag_effector_is_TF;
    int flag_rm_which_TFBS;
    int flag_signal_ctrl_repressor;
    int ngenes;                                             /* the number of loci */
    int nproteins;                                          /* because of gene duplication, the number of proteins and mRNAs can be different 
                                                               from the number of loci. nprotein is the number of elements in protein_pool*/
    int n_output_genes;    
    int N_node_families;
    int N_cisreg_clusters;    
    int which_protein[MAX_GENES];                              /* in case of gene duplication, this array tells the protein corresponding to a given gene id */ 
    
    char cisreg_seq[MAX_GENES][CISREG_LEN];    
    
    /*these apply to protein, not loci. Genes encode the same protein if they only differ in expression kinetics*/
    int N_act;                                              /* number of activators*/
    int N_rep;                                              /* number of repressors*/    
    int protein_identity[MAX_PROTEINS];                  /* entry marks activator (1) or repressor (0)*/        
    char tf_binding_seq[MAX_PROTEINS][TF_ELEMENT_LEN];
    char tf_binding_seq_rc[MAX_PROTEINS][TF_ELEMENT_LEN];       /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
    int protein_pool[MAX_PROTEINS][2][MAX_GENES];               /* element 1 record how many genes/mRNAs producing this protein,ele 2 stores which genes/mRNAs*/
    float Kd[MAX_PROTEINS];
    int N_min_match[MAX_PROTEINS];                              /* minimal number of nucleotides that need to match the consensus binding sequence*/
    
    /*group genes into node family*/
    int which_node_family[MAX_GENES];                          /* a node family include genes that are effectively the same TF (node).
                                                                * In a node family, genes
                                                                * are either all activators or all repressors
                                                                * are either all output protein or non-output protein
                                                                * have identical TFBS (but not necessarily same Kd)  
                                                                */
    int node_family_pool[MAX_GENES][2][MAX_GENES];                            
    
    
    /*these apply to loci*/
    int output_gene_ids[MAX_OUTPUT_GENES];
    int is_output[MAX_GENES];                                   /*is an output protein (-2) or non-output protein (-1)   */
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
    float SE_fitness1;
    float SE_fitness2;
    float avg_fitness;
    float SE_avg_fitness;
    float fitness_measurement[HI_RESOLUTION_RECALC*N_REPLICATES];
    
    /*measurement of network topology*/
    int N_motifs[33]; 
    int N_near_AND_gated_motifs[12];
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
    float max_change_in_probability_of_binding;
    float peak_lvl;
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

/*output buffer*/
typedef struct Output_buffer Output_buffer;
struct Output_buffer
{    
    int step;
    int n_tot_mut;
    int n_mut_at_the_step;
    int n_hit_bound;
    float selection_coefficient;
    float avg_f;
    float f1;
    float f2;
    float se_avg_f;
    float se_f1;
    float se_f2;
    int n_gene;
    int n_output_genes;
    int n_node_families;
    int n_act;
    int n_rep;    
    char mut_type;
    int which_gene;
    int which_protein;
    int which_nuc;
    char new_nuc[3];
    int which_kinetic;
    float new_kinetic;   
    int n_motifs[33];
    int n_near_AND_gated_motifs[12];    
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
                    Mutation *,
                    Selection *,
                    int [MAX_GENES],
                    float [MAX_GENES], 
                    int,
                    RngStream [N_THREADS]);

void evolve_neutrally(  Genotype *,
                        Genotype *,                             
                        Mutation *,
                        Selection *,
                        Selection *,
                        RngStream);

void characterize_network(Genotype *,
                            Genotype *,
                            Mutation *,                            
                            Selection *,
                            int,
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
