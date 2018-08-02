/*  
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2007-2018 Arizona Board of Regents (University of Arizona)
 */
#ifndef FILE_NETSIM_SEEN
#define FILE_NETSIM_SEEN

#include <stdio.h>
#include "RngStream.h"

/******************************************************************************/
/*                            Adjustable knobs                                */
/******************************************************************************/


/*1. Simulation mode*/
/******************************************************************************/
/*Set only one of NEUTRAL, PHENOTYPE, and PERTURB to 1 to enable
 *the corresponding mode. If none of the three is enabled, the default mode is 
 *evolving a TRN under the selection condition specified in main.c*/
#define NEUTRAL 0 //run neutral evolution
#define PHENOTYPE 0 //output the expression of genes over time
#define PERTURB 1//run perturbation analysis


/*2. Runtime control*/ 
/******************************************************************************/
#define MAX_TRIALS 2000 //try at most 2000 mutants at an evolutionary step
#define N_THREADS 10 //the number of parallel OpenMP threads
#define N_REPLICATES 200 //calculate the fitness of a mutant with 200 replicates
#define HI_RESOLUTION_RECALC 5 //calcualte the fitness of a resident with 5*N_REPLICATES replicates
#define OUTPUT_INTERVAL 10 //output Summary_BS.txt every 10 evolutionary steps
#define SAVING_INTERVAL 10 //make a saving point every 10 evoluationary steps
#define OUTPUT_MUTANT_DETAILS 0 //output every mutant genotype and its fitness, whetehr the mutant is accepted
#define OUTPUT_RNG_SEEDS 1 //output the state of random number generator every evolutionary step
#define MAKE_LOG 0 //generate error log
#if MAKE_LOG
#define LOG(...) { FILE *fperror; fperror=fopen("error.txt","a+"); fprintf(fperror, "%s: ", __func__); fprintf (fperror, __VA_ARGS__) ; fflush(fperror); fclose(fperror);} 
#endif


/*3. Biology and evolution settings*/
/******************************************************************************/
#define DIRECT_REG 1 //set it to "1" to allow the signal to directly regulate the effector
#define NO_PENALTY 0 //"1" to remove the penalty to fitness when effector is expressed under a wrong environment
#define RANDOM_COOPERATION_LOGIC 0 //"1" to randomly set whether a gene (including TF genes) is AND-gated-capable at initialization 


/*4. Perturbation options*/
/******************************************************************************/
#if PERTURB 
#define DISABLE_AND_GATE 1 //change the logic of an effector gene
#if DISABLE_AND_GATE
#define WHICH_MOTIF 1 //only one type of motif can be disturbed at a time: 0 for C1-FFL, 1 for FFL-in-diamond, 2 for diamond
#define ADD_STRONG_TFBS 1 //by default we add TFBS as strong as the strongest TFBS that already exists in the cis-reg
#define FORCE_MASTER_CONTROLLED 1 // 1 for master TF controlled; 0 for aux. TF controlled 
#endif
#define FORCE_DIAMOND 0  // change an AND-gated FFL-in-diamond to diamond
#define FORCE_SINGLE_FFL 0 // change an AND-gated FFL-in-diamond to isolated FFL
#endif


/*5. Count additional motifs*/
/******************************************************************************/
#define COUNT_NEAR_AND 0 //count near-AND-gated motifs.
#define COUNT_LONG_ARM 0 //count long-arm C1FFLs


/*6. Analyzing weak TFBSs*/
/******************************************************************************/
/*We can exclude different weak TFBSs when scoring motifs. We recommend doing 
 *this analysis under PHENOTYPE mode.
 */
#define CUT_OFF_MISMATCH_SIG2EFFECTOR 2 //the maximum number of mismatches in TFBSs of the signal in effector genes
#define CUT_OFF_MISMATCH_TF2EFFECTOR 2 //the maximum number of mismatches in TFBSs of TFs in effector genes
#define CUT_OFF_MISMATCH_SIGNAL_TO_TF 2 //the maximum number of mismatches in TFBSs of the signal in TF genes
#define CUT_OFF_MISMATCH_TF_TO_TF 2 //the maximum number of mismatches in TFBSs of TFs in TF genes


/*7. Use an irregular signal in selection condition*/
/*****************************************************************************/
/*An irregular signal can be specified with an external file that describe the signal (see main.c)*/
#define IRREG_SIGNAL 0
#if IRREG_SIGNAL
float signal_profile_matrix[N_THREADS][100][90];
#endif


/*8. Other default settings*/    
/******************************************************************************/
#define MAX_TF_GENES 20 
#define MAX_EFFECTOR_GENES 5  
#define MAX_GENES 25  //total number of genes=effector genes + TF genes (including the signal)     
#define MAX_PROTEINS 25
#define CISREG_LEN 150        //length of cis-regulatory region in base-pairs 
#define CONSENSUS_SEQ_LEN 8      //length of binding element on TF */
#define NMIN 6                //minimal number of nucleotide that matches the binding sequence of a TF in its binding site
                              //DO NOT MAKE NMIN<CONSENSUS_SEQ_LEN/2, OTHERWISE calc_all_binding_sites_copy will make mistake  
#define HIND_LENGTH 3         //default length of hindrance on each side of the binding site,i.e. a tf occupies CONSENSUS_SEQ_LEN+2*HIND_LENGTH
#define MAX_BINDING 10  //MAX_BINDING is the max number of tf that can bind to a promoter plus 1*/
/******************************************************************************/
/*                          End of controlling knobs                          */
/******************************************************************************/






#define N_SIGNAL_TF 1 // number of input signals.
#define RANDOMIZE_SIGNAL2 0 //

/*
 * enum for 'CellState'->transcriptional_state indices
 */
enum PROTEIN_IDENTITY {ACTIVATOR=1, REPRESSOR=0, NON_TF=-1};
enum BOOLEAN {NA=-1, NO=0, YES=1};

/*struct of a TF binding site*/
typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites
{
  int tf_id;         /*which TF */
  float Kd;
  int mis_match;     /* number of mismatched nuc */
  int BS_pos;        /* start position of BS on DNA, always with reference to forward strand */                     
  int N_hindered;    /* the number of BS hindered by this TF when it binds to the current BS */  
};

/*
 * Selection contains two envs; each env states the signal whether 
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
    char tf_seq[MAX_PROTEINS][CONSENSUS_SEQ_LEN];                 /* concensus binding sequences*/
    char tf_seq_rc[MAX_PROTEINS][CONSENSUS_SEQ_LEN];              /* reversed complementary sequence of BS. Used to identify BS on the non-template strand*/
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
    float avg_fitness;
    float fitness1;
    float fitness2; 
    float sq_SE_avg_fitness;
    float sq_SE_fitness1;
    float sq_SE_fitness2;
    float fitness_measurement[HI_RESOLUTION_RECALC*N_REPLICATES];
    
    /*Motifs related*/
    int N_motifs[39];  
    int N_near_AND_gated_motifs[9];
    int N_long_arm_c1ffls[9];
    int TF_in_core_C1ffl[MAX_GENES][MAX_PROTEINS];
    int gene_in_core_C1ffl[MAX_GENES];    
};

/*A mutation is specified by its type, target, and mutant value*/
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

/*Expression levels and instantaneous fitness*/
typedef struct Phenotype Phenotype;
struct Phenotype
{ 
    int timepoint;
    int total_time_points;
    float *protein_concentration;
    float *gene_specific_concentration;
    float *instantaneous_fitness;    
};

/*
 * global variables
 */
/*output files. Initialized in main.c*/
extern char mutation_file[32];
extern char setup_summary[32];
extern char evo_summary[32];

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

/* global function prototypes */
char set_base_pair(float);

void initialize_genotype(Genotype *, int, int, int, int, RngStream) ;

int evolve_under_selection(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_GENES], RngStream, RngStream [N_THREADS]);

void initialize_cache(Genotype *);

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

void perturbation_analysis(Genotype *,
                            Genotype *,
                            Mutation *,
                            Selection *,
                            int [MAX_GENES],
                            float [MAX_GENES],
                            RngStream [N_THREADS]);

void print_mutatable_parameters(Genotype*,int);

void calc_all_binding_sites_copy(Genotype *, int);

void calc_all_binding_sites(Genotype *);

#endif /* !FILE_NETSIM_SEEN */
