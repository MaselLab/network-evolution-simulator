/*  
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2007-2018 Arizona Board of Regents (University of Arizona)
 */

#ifndef CELLULAR_ACTIVITY_H
#define CELLULAR_ACTIVITY_H

#include "netsim.h"
#include "RngStream.h"

#define EPSILON 1.0e-6       /* This is the "0" in changes of TF binding probabilities and of time */
#define TIME_INFINITY 9.99e10 //min
#define TIME_OFFSET 0.01    /* Added to the timing of a fixed event to avoid concurrnce. Unit is minute */

enum TRANSCRIPTIONAL_STATE {REPRESSED, INTERMEDIATE, ACTIVE};

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

/*
 * CellState store the current cellular state, e.g. protein concentration
 */
typedef struct CellState CellState;
struct CellState { 
    float t;
    float cumulative_fitness[MAX_OUTPUT_PROTEINS];                    
    float cumulative_fitness_after_burn_in[MAX_OUTPUT_PROTEINS];          
    float instantaneous_fitness[MAX_OUTPUT_PROTEINS];  
    float cumulative_cost;                    
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
//    float t_burn_in;
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
    float protein_synthesis_index[MAX_GENES];  /*this is N_mRNA*protein_syn_rate/degradation rate.*/
    int transcriptional_state[MAX_GENES];       /*can be REPRESSED, INTERMEDIATE, or ACTIVE */
    
    /*measurement of pulse*/  
    float *sampled_response;
    int N_samples;
};

void initialize_cell(Genotype *, CellState *, Environment *, float, int [MAX_GENES], float [MAX_PROTEINS]);

void do_single_timestep(Genotype *, CellState *, GillespieRates *, Environment *, float, Phenotype *, RngStream) ;

void calc_all_rates(Genotype *, CellState *, GillespieRates *, Environment *, Phenotype *, float, int);

#endif /* EXPRESSION_DYNAMICS_H */

