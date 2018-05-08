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

/*Simulate evolution under selection by default,
 *These are the alternative modes, which can be specified by add -D mode (e.g.
 * -D NEUTRAL) when compile the source code.
 */
#ifndef NEUTRAL
#define NEUTRAL 0
#endif
#ifndef REPLAY
#define REPLAY 0
#endif
#ifndef MODIFY 
#define MODIFY 0
#endif
#ifndef EXTERNAL_SIGNAL
#define EXTERNAL_SIGNAL 0
#endif

/*Runtime control*/ 
#define MAX_MUTATIONS 1000000
#define MAX_TRIALS 2000
#define N_THREADS 10
#define N_REPLICATES 200
#define HI_RESOLUTION_RECALC 5
#define OUTPUT_INTERVAL 10
#define SAVING_INTERVAL 10
#define OUTPUT_MUTANT_DETAILS 0
#define N_TIMP_POINTS 90

/*Miscellaneous settings*/
#define OUTPUT_RNG_SEEDS 1
#define CAUTIOUS 0

/*Biology and evolution settings*/
#ifndef DIRECT_REG  //direct regulation can be enabled by adding -D DIRECT_REG in compilation command
#define DIRECT_REG 0
#endif
#ifndef NO_PENALTY
#define NO_PENALTY 0
#endif
#define minimal_selection_coefficient 1.0e-8
#define NO_REGULATION_COST 0
#define NO_REGULATION 0 // this locks the state of TF genes to Rep
#define RANDOM_COOPERATION_LOGIC 0 
/******************************************************************************/
/*DO NOT MODIFY THE FOLLOWING TWO PARAMETERS!!!
 *I had intended to model more than one signals, but gave it up due to the
 *complexity
 */
/******************************************************************************/
#define N_SIGNAL_TF 1 // number of input signals
#define RANDOMIZE_SIGNAL2 0
/***********************************/

/*Available modifications to network.
 *Enabled under mode MODIFICATION. 
 *Set the desired modification to 1 to enabled it.
 */
#define DISABLE_AND_GATE 0 // change the logic of an effector gene
#if DISABLE_OR_GATE
#define FORCE_MASTER_CONTROLLED 0 // 1 for master TF controlled; 0 for aux. TF controlled 
#endif
#define FORCE_DIAMOND 0  // change an AND-gated FFL-in-diamond to diamond
#define FORCE_SINGLE_FFL 0 // change an AND-gated FFL-in-diamond to isolated FFL

/*Maximum number of genes*/          
#define MAX_TF_GENES 20 
#define MAX_EFFECTOR_GENES 5  
#define MAX_GENES 25  /* total number of genes=effector genes + TF genes (including the signal) */     
#define MAX_PROTEINS 25

/*Cis-reg sequence and TF binding sequence*/
#define CISREG_LEN 150        /* length of cis-regulatory region in base-pairs */
#define TF_ELEMENT_LEN 8      /* length of binding element on TF */
#define NMIN 6                /* minimal number of nucleotide that matches the binding sequence of a TF in its binding site*/
                              /* DO NOT MAKE NMIN<TF_ELEMENT_LEN/2, OTHERWISE calc_all_binding_sites_copy will make mistake*/  
#define HIND_LENGTH 3         /* default length of hindrance on each side of the binding site,i.e. a tf occupies TF_ELEMENT_LEN+2*HIND_LENGTH */
#define MAX_BINDING 10  /* MAX_MODE is the max number of tf that can bind to a promoter plus 1*/



/*
 * enum for 'CellState'->transcriptional_state indices
 */

enum PROTEIN_IDENTITY {ACTIVATOR=1, REPRESSOR=0, NON_TF=-1};
enum BOOLEAN {NA=-1, NO=0, YES=1};

/*struct of a TF binding site*/
typedef struct AllTFBindingSites AllTFBindingSites;
struct AllTFBindingSites {
  int tf_id;         /*which TF */
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
    float fitness_measurement[HI_RESOLUTION_RECALC*N_REPLICATES];
    
    /*Motifs related*/
    int N_motifs[39];  
    int TF_in_core_C1ffl[MAX_GENES][MAX_PROTEINS];
    int gene_in_core_C1ffl[MAX_GENES];
    int N_act_genes; 
    int N_act_genes_reg_by_env;
    int N_act_genes_not_reg_by_env;    
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
extern char RuntimeSumm[32];
extern char output_file[32];
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

/* global function prototypes */

char set_base_pair(float);

void initialize_genotype(Genotype *, int, int, int, int, RngStream) ;

int init_run_pop(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_GENES], RngStream, RngStream [N_THREADS]);

void initialize_cache(Genotype *);

void run_plotting(  Genotype *,
                    Genotype *,
                    Mutation *,
                    Selection *,
                    int [MAX_GENES],
                    float [MAX_GENES],
                    int,
                    RngStream [N_THREADS]);

void evolve_neutrally(  Genotype *,
                        Genotype *,                             
                        Mutation *,                              
                        RngStream);

void plot_alternative_fitness(  Genotype *,
                                Genotype *,
                                Mutation *,
                                Selection *,
                                int [MAX_GENES],
                                float [MAX_GENES]);

void print_mutatable_parameters(Genotype*,int);

void calc_all_binding_sites_copy(Genotype *, int);

void calc_all_binding_sites(Genotype *);

#endif /* !FILE_NETSIM_SEEN */
