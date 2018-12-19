/* 
 * Simulator of yeast transcriptional regulatory network evolution
 * 
 * This file contains the functions to simulate gene expression and
 * instantaneous fitness during development. 
 * 
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2007-2018 Arizona Board of Regents (University of Arizona)
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cellular_activity.h"
#include "lib.h"
#include "numerical.h"

#define DO_NOTHING -2
#define INITIALIZATION -1

/*expression rate parameters*/
static const float TRANSCRIPTINIT=6.75; 
static const float TRANSLATION_INITIATION_TIME=0.5; //min
static const float TRANSCRIPTION_TERMINATION_TIME=1.0; //min
static const float TRANSCRIPTION_ELONGATION_RATE=600.0; //codon/min
static const float TRANSLATION_ELONGATION_RATE=330.0; //codon/min
static const float MAX_REP_TO_INT_RATE=0.92;
static const float BASAL_REP_TO_INT_RATE=0.15;
static const float MAX_INT_TO_REP_RATE=4.11;
static const float BASAL_INT_TO_REP_RATE=0.67;
static const float MAX_INT_TO_ACT_RATE=3.3; 
static const float BASAL_INT_TO_ACT_RATE=0.025;
static const float DEFAULT_UPDATE_INTERVAL=10.0; /*min*/
static const float MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING=0.01;

/*fitness*/
static const float Ne_saturate = 10000.0;
static const float c_transl=2.0e-6;
static const float bmax=1.0; 

/******************************************************************************
 * 
 *                          Private function prototypes
 *
 *****************************************************************************/

static float calc_tprime(Genotype*, CellState*, float [MAX_OUTPUT_GENES], float, float, int [MAX_OUTPUT_GENES], int, int);

static float calc_integral(Genotype *, CellState *, float [MAX_OUTPUT_GENES], float, float, int [MAX_OUTPUT_GENES], int, int);

static void calc_fx_dfx(float, int, float, float*, float*, float*, float*, float*, int);

static void calc_leaping_interval(Genotype*, CellState*, float *, float, int);

static void calc_TF_dist_from_all_BS(Genotype *, CellState*, int);

static int Gillespie_event_mRNA_decay(GillespieRates *, CellState *, Genotype *, RngStream);

static void Gillespie_event_repressed_to_intermediate(GillespieRates *, CellState *, Genotype *, RngStream);

static void Gillespie_event_intermediate_to_repressed(GillespieRates *, CellState *, Genotype *, RngStream); 

static void Gillespie_event_intermediate_to_active(GillespieRates *, CellState *, Genotype *, RngStream);

static void Gillespie_event_active_to_intermediate(Genotype *, CellState *, GillespieRates *, RngStream);

static void Gillespie_event_transcription_init(GillespieRates *, CellState *, Genotype *, float, RngStream);

static int fixed_event_end_translation_init(Genotype *, CellState *, GillespieRates *, Environment *, float *);

static void fixed_event_end_transcription(float *, CellState *, GillespieRates *, Genotype *, Environment *);

static int does_fixed_event_end(CellState*, float);

static int do_fixed_event(Genotype *, CellState *, GillespieRates *, Environment *, Phenotype *, float *, int);

static void calc_fitness(float *, Genotype *, CellState *, Environment *, float [MAX_OUTPUT_GENES], float [MAX_OUTPUT_GENES], float);

static void update_protein_number_and_fitness(Genotype *, CellState *, GillespieRates *, Environment *, float);

static int do_Gillespie_event(Genotype*, CellState *, GillespieRates *, float, RngStream);



/******************************************************************************
 * 
 *                              Global functions
 *
 *****************************************************************************/
/*
 * initialize the cell state with the specified initial protein
 * concentration, mean mRNA number and mRNA decay
 */
void initialize_cell(   Genotype *genotype,
                        CellState *state, 
                        Environment *env,
                        int init_mRNA_number[MAX_GENES],
                        float init_protein_number[MAX_PROTEINS])
{
    int i, j;
    state->t=0.0;
    state->sensitivity[0] = 0.0;	
	state->cumulative_cost = 0.0;
    state->Duration_pulse=0.0;
    state->T_pulse_on=env->t_development;
    state->T_pulse_off=0.0;
    state->Pulse_is_on=0;
    state->first_pulse=1;
    state->N_samples=0;
    state->cumulative_advanced_benefit=0.0;
    state->cumulative_basal_benefit=0.0;
    state->cumulative_cost=0.0;
    state->cumulative_damage=0.0;  
    state->recording_basal_response=0;
    state->found_gradient=0;
    state->position=0;
    state->cumulative_t_in_bias=0;    
    state->moving=0;
    state->cumulative_benefit=0.0;
    
    for(i=0;i<MAX_OUTPUT_PROTEINS;i++)
    {
        state->cumulative_fitness[i] = 0.0;     
        state->cumulative_fitness_after_burn_in[i] = 0.0;   
        state->instantaneous_fitness[i] = 0.0; 
        state->threshold_response_to_bias[i]=0.0;
        state->basal_response[i]=0.0;
    }
    state->effect_of_effector=env->initial_effect_of_effector;

    /* initialize linked tables*/
    state->mRNA_transcr_time_end_head = NULL;
    state->mRNA_transcr_time_end_tail = NULL;
    state->mRNA_transl_init_time_end_head = NULL;
    state->mRNA_transl_init_time_end_tail = NULL;
    state->signal_off_head = NULL;
    state->signal_off_tail = NULL;
    state->signal_on_head = NULL;
    state->signal_on_tail = NULL;
    state->burn_in_growth_rate_head =NULL;
    state->burn_in_growth_rate_tail=NULL;
    state->sampling_point_end_head=NULL;
    state->sampling_point_end_tail=NULL;
    state->last_event_t=0.0;  
    state->change_signal_strength_head=NULL;
    state->change_signal_strength_tail=NULL;
    state->t_to_update_probability_of_binding=TIME_INFINITY;
    state->cell_activated=0;
    /*initialize gene state, mRNA number*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {
        state->transcriptional_state[i]=REPRESSED;       
        state->mRNA_aft_transl_delay_num[i]=init_mRNA_number[i];
        state->mRNA_under_transl_delay_num[i]=0;
        state->mRNA_under_transc_num[i]=0;
        state->last_P_A[i]=0.0;
        state->last_P_R[i]=0.0;
        state->last_P_A_no_R[i]=0.0;;
        state->last_P_NotA_no_R[i]=0.0;
        state->protein_synthesis_index[i]=(float)state->mRNA_aft_transl_delay_num[i]*genotype->protein_syn_rate[i]/genotype->protein_decay_rate[i];
    }       
    /* initiate protein concentration*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
        state->gene_specific_protein_number[i] = init_protein_number[i];    
    for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
    {
        state->protein_number[i]=0.0;        
        for(j=0;j<genotype->protein_pool[i][0][0];j++)
            state->protein_number[i]+=state->gene_specific_protein_number[genotype->protein_pool[i][1][j]];
    }    
    /* deal with the sensor tf*/
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        state->mRNA_aft_transl_delay_num[i]=0;        
        state->mRNA_under_transc_num[i]=0;
        state->mRNA_under_transl_delay_num[i]=0;
        state->gene_specific_protein_number[i]=0.0;
    }        
    /*mark when to start calculating average fitness*/
    if(env->duration_of_burn_in_growth_rate!=0.0)
        add_fixed_event(-1,env->duration_of_burn_in_growth_rate,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail)); 
    else
        add_fixed_event(-1,(float)TIME_INFINITY,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));                
    /*plot protein concentration and fitness vs time*/
    
#if !PHENOTYPE
#if PEAK_SEARCH
    float t;
    int N_data_points; 
    t=env->t_signal_off+TIME_OFFSET;
    N_data_points=(int)((env->t_development-env->t_signal_off)/sampling_interval);
    for(i=0;i<N_data_points;i++)
    {
        add_fixed_event(-1,t,&(state->sampling_point_end_head),&(state->sampling_point_end_tail)); //get a timepoint each minute
        t+=sampling_interval;            
    } 
    state->sampled_response=(float*)calloc(N_data_points,sizeof(float)); 
#elif FAST_SS
    float t;
    int N_data_points; 
    t=env->t_signal_off+TIME_OFFSET;
    N_data_points=(int)((env->t_development-env->t_signal_off)/sampling_interval);
    for(i=0;i<N_data_points;i++)
    {
        add_fixed_event(-1,t,&(state->sampling_point_end_head),&(state->sampling_point_end_tail)); //get a timepoint each minute
        t+=sampling_interval;            
    } 
    state->sampled_response=(float*)calloc(N_data_points,sizeof(float)); 
#elif CHEMOTAXIS
    add_fixed_event(-1,env->t_signal_on-1.0,&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
#endif 
#else
    float t;
    int N_data_points; 
    t=0.0+TIME_OFFSET;
    N_data_points=(int)(env->t_development/sampling_interval);
    for(i=0;i<N_data_points;i++)
    {
        add_fixed_event(-1,t,&(state->sampling_point_end_head),&(state->sampling_point_end_tail)); //get a timepoint each minute
        t+=sampling_interval;            
    } 
#endif    
}

/*
 * Calculate the rates of all Gillespie events
 */
void calc_all_rates(Genotype *genotype,
                    CellState *state,
                    GillespieRates *rates,
                    Environment *env,
                    int UPDATE_WHAT)
{
    int i,cluster_id,gene_id;
    int concurrent;
    float t_to_update_probability_of_binding,interval_to_update_probability_of_binding;  
    float diff_PA,diff_PR,diff_PnotAnoR,diff_PAnoR,diff_max;  
    /* reset rates */
    rates->total_mRNA_decay_rate=0.0;
    rates->total_active_to_intermediate_rate=0.0;
    rates->total_repressed_to_intermediate_rate=0.0;
    rates->total_intermediate_to_repressed_rate=0.0;
    rates->total_N_gene_transcript_initiated=0;
    rates->total_intermediate_to_active_rate=0.0;
    rates->total_Gillespie_rate=0.0;    
    for(i=0;i<genotype->ngenes;i++)
    {
        rates->repressed_to_intermediate_rate[i]=0.0;
        rates->intermediate_to_repressed_rate[i]=0.0;
        rates->intermediate_to_active_rate[i]=0.0;
        rates->active_to_intermediate_rate[i]=0.0;
        rates->mRNA_decay_rate[i]=0.0;
        rates->transcript_initiation_state[i]=0;
        state->P_A[i]=0.0;
        state->P_R[i]=0.0;
        state->P_A_no_R[i]=0.0;
        state->P_NotA_no_R[i]=0.0;
    } 
    /* update mRNA decay rates*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        rates->mRNA_decay_rate[i] = genotype->mRNA_decay_rate[i] * (state->mRNA_aft_transl_delay_num[i] + state->mRNA_under_transl_delay_num[i]);
        rates->total_mRNA_decay_rate += rates->mRNA_decay_rate[i];        
    }
    /*update probability of binding configurations that activates expression
     * and use it to update other rates*/
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {    
        cluster_id=genotype->which_cluster[i];        
        if(genotype->cisreg_cluster_pool[cluster_id][1][0]!=i)  /*if this gene does not have a unique cis-reg sequence*/
        {                
            state->P_A[i]=state->P_A[genotype->cisreg_cluster_pool[cluster_id][1][0]]; /* copy TF distribution from elsewhere*/
            state->P_R[i]=state->P_R[genotype->cisreg_cluster_pool[cluster_id][1][0]];
            state->P_A_no_R[i]=state->P_A_no_R[genotype->cisreg_cluster_pool[cluster_id][1][0]];
            state->P_NotA_no_R[i]=state->P_NotA_no_R[genotype->cisreg_cluster_pool[cluster_id][1][0]];
        }
        else /* otherwise, we need to calc the ratio*/
        {
            if(genotype->N_act_BS[i]!=0 || genotype->N_rep_BS[i]!=0)                
                calc_TF_dist_from_all_BS(genotype, state, i);
            else
            {
                state->P_A[i]=0.0;     
                state->P_R[i]=0.0;
                state->P_A_no_R[i] = 0.0;
                state->P_NotA_no_R[i]=0.0;
            }
        }
        
        /* calc other rates*/
        switch (state->transcriptional_state[i])
        {
            case REPRESSED:
                rates->repressed_to_intermediate_rate[i]=state->P_A[i]*(MAX_REP_TO_INT_RATE-BASAL_REP_TO_INT_RATE)+BASAL_REP_TO_INT_RATE;             
                rates->total_repressed_to_intermediate_rate+=rates->repressed_to_intermediate_rate[i];
                rates->intermediate_to_repressed_rate[i]=0.0;
                rates->intermediate_to_active_rate[i]=0.0;
                rates->active_to_intermediate_rate[i]=0.0;
                break;                
            case INTERMEDIATE:
                rates->intermediate_to_repressed_rate[i]=state->P_R[i]*(MAX_INT_TO_REP_RATE-BASAL_INT_TO_REP_RATE)+BASAL_INT_TO_REP_RATE;
                rates->total_intermediate_to_repressed_rate+=rates->intermediate_to_repressed_rate[i];
                rates->intermediate_to_active_rate[i]=MAX_INT_TO_ACT_RATE*state->P_A_no_R[i]+BASAL_INT_TO_ACT_RATE*state->P_NotA_no_R[i];
                rates->total_intermediate_to_active_rate+=rates->intermediate_to_active_rate[i]; 
                rates->active_to_intermediate_rate[i]=0.0;
                rates->repressed_to_intermediate_rate[i]=0.0;
                break;                
            case ACTIVE: 
                rates->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i];
                rates->total_active_to_intermediate_rate+=rates->active_to_intermediate_rate[i]; 
                rates->transcript_initiation_state[i]= 1;
                rates->total_N_gene_transcript_initiated+=1;
                rates->intermediate_to_repressed_rate[i]=0.0; 
                rates->repressed_to_intermediate_rate[i]=0.0;
                rates->intermediate_to_active_rate[i]=0.0;
                break;
        }
    }
    rates->total_Gillespie_rate+=rates->total_intermediate_to_repressed_rate;
    rates->total_Gillespie_rate+=rates->total_intermediate_to_active_rate;
    rates->total_Gillespie_rate+=rates->total_repressed_to_intermediate_rate;
    rates->total_Gillespie_rate+=rates->total_mRNA_decay_rate;
    rates->total_Gillespie_rate+=rates->total_active_to_intermediate_rate;
    rates->total_Gillespie_rate+=(float)rates->total_N_gene_transcript_initiated*TRANSCRIPTINIT;  
    
    /*Check if Pact needs to be updated more or less often*/ 
    if(UPDATE_WHAT!=INITIALIZATION && state->cell_activated==1)
    {  
        diff_max=0.0;          
        for(cluster_id=N_SIGNAL_TF;i<genotype->N_cisreg_clusters;cluster_id++) //check if Pact changes too much
        {
            gene_id=genotype->cisreg_cluster_pool[cluster_id][1][0];
            diff_PA=fabs(state->P_A[gene_id]-state->last_P_A[gene_id]);            
            diff_PAnoR=fabs(state->P_A_no_R[gene_id]-state->last_P_A_no_R[gene_id]); 
            diff_PnotAnoR=fabs(state->P_NotA_no_R[gene_id]-state->last_P_NotA_no_R[gene_id]);  
            diff_PR=fabs(state->P_R[gene_id]-state->last_P_R[gene_id]);
            diff_max=(diff_max>diff_PA)?diff_max:diff_PA;
            diff_max=(diff_max>diff_PR)?diff_max:diff_PR;
            diff_max=(diff_max>diff_PAnoR)?diff_max:diff_PAnoR;
            diff_max=(diff_max>diff_PnotAnoR)?diff_max:diff_PnotAnoR;           
        }
        if(diff_max<EPSILON)
            interval_to_update_probability_of_binding=DEFAULT_UPDATE_INTERVAL;
        else
            interval_to_update_probability_of_binding=MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING/diff_max*(state->t-state->last_event_t);
        if(UPDATE_WHAT!=DO_NOTHING)          
            calc_leaping_interval(genotype,state,&interval_to_update_probability_of_binding,env->t_development,UPDATE_WHAT);  
    
        /*Update the next time that Pact will be updated mandatorily*/
        t_to_update_probability_of_binding=state->t+interval_to_update_probability_of_binding;
        concurrent=check_concurrence(state, t_to_update_probability_of_binding);
        while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
        {
            t_to_update_probability_of_binding+=TIME_OFFSET;
            concurrent=check_concurrence(state, t_to_update_probability_of_binding);        
        }
        state->t_to_update_probability_of_binding=t_to_update_probability_of_binding;
    }
    /*Keep a copy of Pact and time for comparison at next time Pact is updated*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        state->last_P_A[i]=state->P_A[i];     
        state->last_P_R[i]=state->P_R[i];
        state->last_P_A_no_R[i] = state->P_A_no_R[i];
        state->last_P_NotA_no_R[i]=state->P_NotA_no_R[i];
    }
    state->last_event_t=state->t;  
}

/*
 * run the model for a specified cell for a single timestep:
 */
void do_single_timestep(Genotype *genotype, 
                        CellState *state,                         
                        GillespieRates *rates,                        
                        Environment *env,   
                        Phenotype *timecourse,
                        RngStream RS) 
{    
    int event, UPDATE_WHAT;     
    float fixed_time; 
    float dt;
    float x; 
    
    /* draw random number */
    x = expdev(RS);       
    dt = x/rates->total_Gillespie_rate;    
    /* check if a fixed event occurs during dt, or in tdevelopment if running for a fixed development time */
    fixed_time = (state->t+dt<env->t_development)?(state->t+dt):env->t_development;
    event = does_fixed_event_end(state, fixed_time);
    while(event!=0)
    {           
        /*after doing fixed event, return a flag to indicate whether mandatorily update Pact or Prep*/
        UPDATE_WHAT=do_fixed_event(genotype, state, rates, env, timecourse, &dt, event);       
        /* advance time by the dt */  
        state->t += dt;    
        /* we've been running with rates->total_Gillespie_rate for dt, so substract it from x*/   
        x -= dt*rates->total_Gillespie_rate;  
        /* update rates->total_Gillespie_rate and compute a new dt */  
        calc_all_rates(genotype, state, rates, env, UPDATE_WHAT);      
        dt = x/rates->total_Gillespie_rate;
        /*deal with rounding error*/
        if(dt<0.0)
        {  	
#if KEEP_LOG // if enabled, error file can be huge, because this rounding error can happen very often            
            LOG("Negative dt!\n");    
#endif 
            dt=TIME_OFFSET; 
        }
        fixed_time = (state->t+dt<env->t_development)?(state->t+dt):env->t_development;
        /* check to see there aren't more fixed events to do */
        event = does_fixed_event_end(state, fixed_time);                                    
    } 
    /* no remaining fixed events to do in dt, now do stochastic events */  
    /* if we haven't already reached end of development with last delta-t*/          
    if (state->t+dt < env->t_development)
    {        
        /*update protein concentration and fitness after dt*/
        update_protein_number_and_fitness(genotype, state, rates, env, dt); 
        /*do Gillespie event*/
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, RS);       
        /* Gillespie step: advance time to next event at dt */
        state->t += dt;
        calc_all_rates(genotype,state,rates,env,UPDATE_WHAT);        
    } 
    else 
    { 
        /* do remaining dt */
        dt = env->t_development - state->t;
        /* the final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, env, dt);
        /* advance to end of development (this exits the outer while loop) */
        state->t = env->t_development;
#if SELECT_SENSITIVITY_AND_PRECISION
        state->t = env->tdevelopment;
		state->sensitivity[2] = state->sensitivity[0];
		state->precision[2] = state->precision[0];
#endif
    }
}



/*****************************************************************************
 * 
 *                             Private functions
 * 
 *****************************************************************************/

/*
 * compute t' factor used in the integration of fitness
 * t' is the time the effector protein increases or decreases to a given amount
 */
static float calc_tprime(Genotype *genotype, CellState* state, float number_of_selection_protein_bf_dt[MAX_GENES], float dt, float given_amount, int gene_ids[MAX_GENES], int start, int n_copies) 
{    
    int i;      
    float protein_synthesis_rate[n_copies],protein_decay_rate[n_copies];
    
    for(i=start;i<start+n_copies;i++)
    {
        protein_decay_rate[i]=genotype->protein_decay_rate[gene_ids[i]];
        protein_synthesis_rate[i]=state->protein_synthesis_index[gene_ids[i]]*protein_decay_rate[i];
    }       
    return rtsafe(&calc_fx_dfx, n_copies, given_amount, number_of_selection_protein_bf_dt, protein_synthesis_rate, protein_decay_rate, 0.0, dt,start); 
}

/*
 * calculate f(x)-Pp_s and f'(x),
 * f(x) is the number of effector protein molecules at time x
 */
static void calc_fx_dfx(float x, int n_copies, float given_amount, float *intial_protein_number, float *protein_synthesis_rate, float *protein_decay_rate,float *fx, float *dfx, int start)
{
    int i;    
    *fx=0;
    *dfx=0;    
    for(i=0;i<n_copies;i++)
    {
        *fx+=(intial_protein_number[i+start]-protein_synthesis_rate[i]/protein_decay_rate[i])*exp(-protein_decay_rate[i]*x)+protein_synthesis_rate[i]/protein_decay_rate[i];
        *dfx+=(protein_synthesis_rate[i]-protein_decay_rate[i]*intial_protein_number[i+start])*exp(-protein_decay_rate[i]*x);
    }
    *fx-=given_amount;    
}

/*
 * calculate F(delta_t)/Ne_sat. F(x) is the integral of f(x) over delta_t.
 * f(x) is the number of effector protein molecules at time x
 */
static float calc_integral(Genotype *genotype, CellState *state, float initial_protein_number[MAX_OUTPUT_GENES], float dt, float saturate_protein_number, int gene_ids[MAX_OUTPUT_GENES], int start, int n_copies)
{
    int i;
    float integral=0.0,ect_minus_one;     
        
    for(i=start;i<start+n_copies;i++)
    {
        ect_minus_one=exp(-genotype->protein_decay_rate[gene_ids[i]]*dt)-1.0;    
        integral+=(state->protein_synthesis_index[gene_ids[i]]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]-
                initial_protein_number[i]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]+ 
                state->protein_synthesis_index[gene_ids[i]]*dt);
    }
    return integral/saturate_protein_number;
}

/*
 * return the instantaneous fitness given the current cell state and environment,
 * also return the integrated fitness
 */
static void calc_fitness(float *integrated_fitness,
                        Genotype *genotype,
                        CellState *state,
                        Environment *env,
                        float number_of_selection_protein_bf_dt[MAX_OUTPUT_GENES],
                        float number_of_selection_protein_aft_dt[MAX_OUTPUT_GENES],
                        float dt)
{
    int i,j,protein_id,counter;
    int N_output;
    int gene_copy_id[MAX_OUTPUT_GENES];
    int N_gene_copies[MAX_OUTPUT_PROTEINS];
//    float dt_prime, dt_rest, t_on, t_off, t_on_copy,t_off_copy;   
//    float penalty=bmax/Ne_saturate;   
    float cost_of_expression;     
//    float Ne, Ne_next;   
//    int concurrent;
//    float endtime;  
    
//    float instantaneous_fitness=0.0;  /* this is returned from the function */
    float total_protein_syn_rate = 0.0;    
//    float dt_prime;   
//    float cost_of_expression;     
//    float Ne_next=state->protein_number[genotype->nproteins-1];   
//    float Ne=0.0;
//    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++) 
//       Ne+=number_of_selection_protein_bf_dt[i];
                      
    /* compute the total cost of translation across all genes  */
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {     
        total_protein_syn_rate += (genotype->protein_syn_rate[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->protein_syn_rate[i]*(float)state->mRNA_under_transl_delay_num[i])*(float)genotype->locus_length[i]/236.0;
    }
    cost_of_expression=total_protein_syn_rate*c_transl;
    state->cumulative_cost+=cost_of_expression*dt;
    
    /*list the genes that encode effectors*/
    counter=0;
//    for(i=0;i<genotype->n_output_proteins;i++)
//    {        
//        protein_id=genotype->output_protein_ids[i];            
//        for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)         
//        {
//            gene_copy_id[counter]=genotype->protein_pool[protein_id][1][j];
//            counter++;
//        }
//    }    
#if POOL_EFFECTORS 
    N_output=1;
    N_gene_copies[0]=genotype->n_output_genes;    
#else
    N_output=genotype->n_output_proteins;
    for(i=0;i<genotype->n_output_proteins;i++)     
        N_gene_copies[i]=genotype->protein_pool[genotype->output_protein_ids[i]][0][0];    
#endif
    
#if SELECT_SENSITIVITY_AND_PRECISION
	state->sensitivity[0] = (Ne > Ne_next) ? Ne : Ne_next;

	if (effect_of_effector=='d') //calculate precision		
		state->precision[0] = (Ne > Ne_next) ? Ne_next : Ne;	
		
#elif SELECT_ON_DURATION
    /*Ne increases over Ne_sat*/
    t_on=state->T_pulse_on;
    t_off=state->T_pulse_off;
    if (Ne <= Ne_saturate && Ne_next>=Ne_saturate)
    {
        t_on= calc_tprime(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate);
        t_on+=state->t;
        state->Pulse_is_on = 1;
        if(state->first_pulse)
        {
            state->T_pulse_on=t_on;
            state->first_pulse=0;
        }
    }
    /*Ne decreases over Ne_sat*/
    if (Ne >= Ne_saturate && Ne_next <= Ne_saturate && state->Pulse_is_on)
    {
        t_off = calc_tprime(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate);
        t_off+=state->t;
        state->Pulse_is_on = 0;
        if (state->Duration_pulse < t_off - t_on)
        {
            state->T_pulse_on = t_on;
            state->T_pulse_off = t_off;
            state->Duration_pulse = t_off - t_on;
        }
    }
    /*Ne never decreases over Ne_sat*/
    if(state->Pulse_is_on && state->t+dt>=tdevelopment)
    {
        t_off=tdevelopment;
        state->Pulse_is_on=0;
        if (state->Duration_pulse < t_off - t_on)
        {
            state->T_pulse_on = t_on;
            state->T_pulse_off = t_off;
            state->Duration_pulse = t_off - t_on;
        }
    }
    
#elif REALLY_COMPLETECATE
    switch (effect_of_effector)/* effector is beneficial!*/
    {
        case 'b': /* effector is beneficial!*/
            if(Ne>Ne_next)//decrease in effector protein
            {
                if(Ne_next>=Ne_saturate) //too many effector throughout
                {   
                    state->cumulative_advanced_benefit +=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0);
                }
                else if(Ne<=Ne_saturate) // not enough effector throughout
                {
                    state->cumulative_basal_benefit += calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0);
                }
                else // bf dt_prime, the benefit saturates
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
                    state->cumulative_advanced_benefit+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0);                
                    state->cumulative_basal_benefit+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)-
                                                      calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0);
                }                    
            }
            else // increase in effector protein
            {
                if(Ne>=Ne_saturate) //too many effector throughout
                {
                    state->cumulative_advanced_benefit+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0);
                }   
                else if(Ne_next<=Ne_saturate)// not enough effector throughout
                {
                    state->cumulative_basal_benefit+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0);
                }
                else //Aft dt_prime, the benefit saturates
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
                    state->cumulative_basal_benefit+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0);                
                    state->cumulative_advanced_benefit+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)-
                                                      calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0);
                    if(state->first_pulse==1) //start the onset of the pulse
                    {
                        endtime=state->t+dt_prime+opt_pulse_duration;
                        concurrent=check_concurrence(   endtime,
                                                        state->mRNA_transl_init_time_end_head,
                                                        state->mRNA_transcr_time_end_head,
                                                        state->signal_on_head,
                                                        state->signal_off_head,
                                                        state->burn_in_growth_rate_head,
                                                        state->->last_P_NotA_no_R,
                                                        state->change_signal_strength_head);
                        while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
                        {
                            endtime+=TIME_OFFSET;
                            concurrent=check_concurrence(   endtime,
                                                            state->mRNA_transl_init_time_end_head,
                                                            state->mRNA_transcr_time_end_head,
                                                            state->signal_on_head,
                                                            state->signal_off_head,
                                                            state->burn_in_growth_rate_head,
                                                            state->->last_P_NotA_no_R,
                                                            state->change_signal_strength_head);        
                        }  
                        add_fixed_event(-1,endtime,&(state->signal_off_head),&(state->signal_off_tail)); // mark when the pulse should be turned off
                        state->first_pulse=0;
                    }
                }                
            } 
            break;
    
        case 'd': /* effector is deleterious! */            
            if(Ne>Ne_next)//decrease in effector protein
            {
                if(Ne_next>=Ne_saturate) //too many effector throughout
                {
                    state->cumulative_damage+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0);
                }                
                else if(Ne>=Ne_saturate) // 
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
                    state->cumulative_damage+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)-
                                              calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0);                    
                }                    
            }
            else // increase in effector protein
            {
                if(Ne>=Ne_saturate) //too many effector throughout
                {
                    state->cumulative_damage+=calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0);
                }   
                else if(Ne_next>=Ne_saturate)// not enough effector throughout
                {                    
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
                    state->cumulative_damage+= calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0);
                }                
            }   
            break;
    } 
#elif CHEMOTAXIS    
    counter=0;
    for(i=0;i<N_output;i++)        
    {    
        Ne=0.0;
        Ne_next=0.0;
        for(j=counter;j<counter+N_gene_copies[i];j++)
        {            
            Ne_next+=number_of_selection_protein_aft_dt[j];
            Ne+=number_of_selection_protein_bf_dt[j];
        } 
        /*calculate the cumulative level of the effector protein during the 
         *1 minute before the gradient is found*/        
        if(state->recording_basal_response)
            state->basal_response[i]+=calc_integral(state,number_of_selection_protein_bf_dt,dt,1.0,gene_copy_id,counter,N_gene_copies[i]);
        
        /*calculate if moving upstream of the gradient*/
        if(state->found_gradient)
        {    
            switch (state->position)
            {
                case 0: //still in poor media
                    if(Ne_next>=state->threshold_response_to_bias[i]) //expression increases
                    {
                        dt_prime=calc_tprime(state,number_of_selection_protein_bf_dt,dt,state->threshold_response_to_bias[i],gene_copy_id,counter,N_gene_copies[i]); 
                        state->cumulative_benefit+=dt_prime*benefit1;
                        state->position=1; //moving to rich media
                        state->cumulative_benefit+=(dt-dt_prime)*benefit2;
                    }
                    else
                        state->cumulative_benefit+=dt*benefit1;                        
                    break;
                case 1: //in rich media
//                    switch (state->moving)
//                    {
//                        case 0:
//                            break;
//                        case 1:
//                            break;
//                    }
                    
                    if(Ne>Ne_next)//decrease in effector protein
                    {
                        if(Ne_next>=state->threshold_response_to_bias[i]) //too many effector throughout                   
                            state->cumulative_t_in_bias+=dt; 
                        else if(Ne>state->threshold_response_to_bias[i]) // bf dt_prime, the benefit saturates 
                            state->cumulative_t_in_bias+=calc_tprime(state,number_of_selection_protein_bf_dt,dt,state->threshold_response_to_bias[i],gene_copy_id,counter,N_gene_copies[i]);  
                    }
                    else // increase in effector protein
                    {
                        if(Ne>=state->threshold_response_to_bias[i]) //too many effector throughout                
                            state->cumulative_t_in_bias+=dt;
                        else if(Ne_next>state->threshold_response_to_bias[i]) //Aft dt_prime, the benefit saturates
                            state->cumulative_t_in_bias+=dt-calc_tprime(state,number_of_selection_protein_bf_dt,dt,state->threshold_response_to_bias[i],gene_copy_id,counter,N_gene_copies[i]);  
                    }                     
                    
                    /*swimming in bias for too long?*/
                    if(state->cumulative_t_in_bias>max_t_bias)//Yes. 
                    {
                        state->cumulative_benefit+=benefit2*(max_t_bias-(state->cumulative_t_in_bias-dt));
                        state->position=2; //Then back to poor media
                        state->cumulative_benefit+=benefit1*(state->cumulative_t_in_bias-max_t_bias);
                    }
                    else                    
                        state->cumulative_benefit+=dt*benefit2;
                    
                    break;
                case 2: //back to poor media
                    state->cumulative_benefit+=dt*benefit1;
                    break;
            }      
        }      
        else        
            state->cumulative_benefit+=dt*benefit1;        
        counter+=N_gene_copies[i];
    }     
#elif PEAK_SEARCH  
#elif !FAST_SS
    counter=0;
    for(i=0;i<N_output;i++)        
    {    
        Ne=0.0;
        Ne_next=0.0;
        for(j=counter;j<counter+N_gene_copies[i];j++)
        {            
            Ne_next+=number_of_selection_protein_aft_dt[j];
            Ne+=number_of_selection_protein_bf_dt[j];
        }        
        switch (effect_of_effector)
        {
            case 'b': /* effector is beneficial!*/
                if(Ne>Ne_next)//decrease in effector protein
                {
                    if(Ne_next>=Ne_saturate) //too many effector throughout
                    {
                        state->cumulative_fitness[i]+=dt*bmax;
                    }
                    else if(Ne<=Ne_saturate) // not enough effector throughout
                    {
                        state->cumulative_fitness[i]+= bmax*calc_integral(state,number_of_selection_protein_bf_dt,dt,Ne_saturate,gene_copy_id,counter,N_gene_copies[i]);                                            
                    }
                    else // bf dt_prime, the benefit saturates
                    {
                        dt_prime=calc_tprime(state,number_of_selection_protein_bf_dt,dt,Ne_saturate,gene_copy_id,counter,N_gene_copies[i]); 
                        state->cumulative_fitness[i]+=bmax*dt_prime+bmax*(calc_integral(state, number_of_selection_protein_bf_dt, dt, Ne_saturate,gene_copy_id,counter,N_gene_copies[i])-
                                                      calc_integral(state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate,gene_copy_id,counter,N_gene_copies[i]));
                    }                    
                }
                else // increase in effector protein
                {
                    if(Ne>=Ne_saturate) //too many effector throughout
                    {
                        state->cumulative_fitness[i]+=dt*bmax;
                    }   
                    else if(Ne_next<=Ne_saturate)// not enough effector throughout
                    {
                        state->cumulative_fitness[i]+=bmax*calc_integral(state, number_of_selection_protein_bf_dt, dt, Ne_saturate,gene_copy_id,counter,N_gene_copies[i]);                                                   
                    }
                    else //Aft dt_prime, the benefit saturates
                    {
                        dt_prime=calc_tprime(state,number_of_selection_protein_bf_dt,dt,Ne_saturate,gene_copy_id,counter,N_gene_copies[i]); 
                        state->cumulative_fitness[i]+=bmax*(dt-dt_prime)+bmax*calc_integral(state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate,gene_copy_id,counter,N_gene_copies[i]);                                                   
                    }                
                } 
                /* compute instantaneous growth rate at t */
                if (Ne_next < Ne_saturate)
                    state->instantaneous_fitness[i] = bmax*Ne_next/Ne_saturate;
                else
                    state->instantaneous_fitness[i] = bmax;            
                break;

            case 'd': /* effector is deleterious! */      
        #if !NO_PENALTY
                if(Ne>Ne_next)//decrease in effector protein
                {
                    if(Ne_next>=Ne_saturate) //too many effector throughout
                    {
                        state->cumulative_fitness[i]+=0.0;
                    }
                    else if(Ne<=Ne_saturate) // not enough effector throughout
                    {
                        state->cumulative_fitness[i]+=bmax*dt-penalty*calc_integral(state, number_of_selection_protein_bf_dt, dt, 1.0,gene_copy_id,counter,N_gene_copies[i]);                                                    
                    }
                    else // aft dt_prime, the benefit becomes positive
                    {
                        dt_prime=calc_tprime(state,number_of_selection_protein_bf_dt,dt,Ne_saturate,gene_copy_id,counter,N_gene_copies[i]); 
                        state->cumulative_fitness[i]+=bmax*(dt-dt_prime)-penalty*(calc_integral(state, number_of_selection_protein_bf_dt, dt, 1.0,gene_copy_id,counter,N_gene_copies[i])-
                                                      calc_integral(state, number_of_selection_protein_bf_dt, dt_prime, 1.0,gene_copy_id,counter,N_gene_copies[i]));                                                                     
                    }                    
                }
                else // increase in effector protein
                {
                    if(Ne>=Ne_saturate) //too many effector throughout
                    {
                        state->cumulative_fitness[i]+=0.0;
                    }   
                    else if(Ne_next<=Ne_saturate)// not enough effector throughout
                    {
                        state->cumulative_fitness[i]+=bmax*dt-penalty*calc_integral(state, number_of_selection_protein_bf_dt, dt, 1.0,gene_copy_id,counter,N_gene_copies[i]);
                    }
                    else //Aft dt_prime, the benefit becomes zero
                    {
                        dt_prime=calc_tprime(state,number_of_selection_protein_bf_dt,dt,Ne_saturate,gene_copy_id,counter,N_gene_copies[i]); 
                        state->cumulative_fitness[i]+=bmax*dt_prime-penalty*calc_integral(state, number_of_selection_protein_bf_dt, dt_prime, 1.0,gene_copy_id,counter,N_gene_copies[i]);
                    }                
                } 
                if(Ne_next<Ne_saturate)
                    state->instantaneous_fitness[i] = bmax - penalty*Ne_next;
                else
                    state->instantaneous_fitness[i] = 0.0;
                break; 
        #else
                state->cumulative_fitness[i] +=bmax*dt;
                state->instantaneous_fitness[i]=bmax;            
        #endif
        }
        counter+=N_gene_copies[i];
    } 
#endif     
    
//    switch (state->effect_of_effector)
//    {
//        case 'b': /* effector is beneficial!*/
//            if(Ne>Ne_next)//decrease in effector protein
//            {
//                if(Ne_next>=Ne_saturate) //too many effector throughout
//                {
//                    *integrated_fitness =dt*(bmax-cost_of_expression);
//                }
//                else if(Ne<=Ne_saturate) // not enough effector throughout
//                {
//                    *integrated_fitness = bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
//                                            -cost_of_expression*dt;
//                }
//                else // bf dt_prime, the benefit saturates
//                {
//                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
//                    *integrated_fitness = bmax*dt_prime+bmax*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)-
//                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate))-
//                                                    cost_of_expression*dt;                    
//                }                    
//            }
//            else // increase in effector protein
//            {
//                if(Ne>=Ne_saturate) //too many effector throughout
//                {
//                    *integrated_fitness =dt*(bmax-cost_of_expression);
//                }   
//                else if(Ne_next<=Ne_saturate)// not enough effector throughout
//                {
//                    *integrated_fitness = bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
//                                                    -cost_of_expression*dt;
//                }
//                else //Aft dt_prime, the benefit saturates
//                {
//                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
//                    *integrated_fitness = bmax*(dt-dt_prime)+bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate)-
//                                                    cost_of_expression*dt;
//                }                
//            } 
//            /* compute instantaneous fitness at t */
//            if (Ne_next < Ne_saturate)
//                instantaneous_fitness = bmax*Ne_next/Ne_saturate;
//            else
//                instantaneous_fitness = bmax;
//            break;
//    
//        case 'd': /* effector is deleterious! */      
//#if !NO_PENALTY
//            if(Ne>Ne_next)//decrease in effector protein
//            {
//                if(Ne_next>=Ne_saturate) //too many effector throughout
//                {
//                    *integrated_fitness =0.0-dt*cost_of_expression;
//                }
//                else if(Ne<=Ne_saturate) // not enough effector throughout
//                {
//                    *integrated_fitness = bmax*dt-bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
//                                                    -cost_of_expression*dt;
//                }
//                else // aft dt_prime, the benefit becomes positive
//                {
//                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
//                    *integrated_fitness = bmax*(dt-dt_prime)-bmax*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)-
//                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate))-
//                                                    cost_of_expression*dt;                    
//                }                    
//            }
//            else // increase in effector protein
//            {
//                if(Ne>=Ne_saturate) //too many effector throughout
//                {
//                    *integrated_fitness =0.0-dt*cost_of_expression;
//                }   
//                else if(Ne_next<=Ne_saturate)// not enough effector throughout
//                {
//                    *integrated_fitness = bmax*dt-bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
//                                                    -cost_of_expression*dt;
//                }
//                else //Aft dt_prime, the benefit becomes zero
//                {
//                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
//                    *integrated_fitness = bmax*dt_prime-bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate)-
//                                                    cost_of_expression*dt;
//                }                
//            } 
//            if(Ne_next<Ne_saturate)
//                instantaneous_fitness = bmax - bmax/Ne_saturate*Ne_next;
//            else
//                instantaneous_fitness = 0.0;            
//#else
//            *integrated_fitness=(bmax-cost_of_expression)*dt;
//            instantaneous_fitness=bmax;            
//#endif
//            break; 
//    } 
//    /* and instantaneous integrated rate */
//    instantaneous_fitness -= cost_of_expression;
//    /* return the instantaneous fitness */
//    return instantaneous_fitness;
}

/*Calculate probability of binding configurations*/
static void calc_TF_dist_from_all_BS(Genotype *genotype, CellState *state, int gene_id) 
{
    int max_N_binding_act=genotype->max_unhindered_sites[gene_id][1]+1; //Binding configurations can contain at most x activators, plus 1 type of configurations that don't have activators at all. 
    int max_N_binding_rep=genotype->max_unhindered_sites[gene_id][2]+1; //Binding configurations can contain at most y repressors, plus 1 type of configurations that don't have repressors at all. 
    double ratio_matrices[genotype->binding_sites_num[gene_id]+1][max_N_binding_rep][max_N_binding_act]; 
    double sum;    
    register double product_of_freq; 
    register float cache_Kd;
    int pos_of_last_record;  
    int pos_of_mat_nH;
    int pos_next_record;
    int i,j,k,m,n;
    double temp;
    AllTFBindingSites *BS_info;
    float *protein_number;
    protein_number=&(state->protein_number[0]);
    
    BS_info=genotype->all_binding_sites[gene_id];
    
    /* initializing matrices to all zeros */
    for(i=0;i<genotype->max_hindered_sites[gene_id]+1;i++)
    {
        for(j=0;j<max_N_binding_rep;j++) 
        {
            for(k=0;k<max_N_binding_act;k++)
                ratio_matrices[i][j][k]=0.0;
        }
    }   
    /* body of the forward algorithm*/    
    pos_next_record=0; //where in the ratio_matrices to put the next record
    ratio_matrices[pos_next_record][0][0]=BS_info[0].Kd;   
    /*calculate distribution based on the first BS*/
    if(genotype->locus_specific_TF_behavior[gene_id][BS_info[0].tf_id]==ACTIVATOR)// if a activator binds to BS 0  
        ratio_matrices[pos_next_record][0][1]=protein_number[BS_info[0].tf_id];
    else    
        ratio_matrices[pos_next_record][1][0]=protein_number[BS_info[0].tf_id]; 
    /*keep calculating distribution from the remaining BS*/
    for(m=1;m<genotype->binding_sites_num[gene_id];m++)
    {
        pos_next_record++;
        /*If binding of site m blocks other binding sites*/
        product_of_freq = protein_number[BS_info[m].tf_id]; 
        if(BS_info[m].N_hindered!=0) 
        {
            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
                product_of_freq*=BS_info[n].Kd;            
        }
        cache_Kd=BS_info[m].Kd;
        /*Check whether m is a site of activator or repressor*/
        switch(genotype->locus_specific_TF_behavior[gene_id][BS_info[m].tf_id])
        {
            case ACTIVATOR: // a BS of activators              
                if(m-BS_info[m].N_hindered!=0)//if binding of m does not block all of the BS evaluated before
                {      
                    /*find matrix(n-H)*/
                    pos_of_mat_nH=pos_next_record-BS_info[m].N_hindered-1; 
                   /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;
                    for(i=0;i<max_N_binding_rep;i++)
                        for(j=1;j<max_N_binding_act;j++) 
                            ratio_matrices[pos_next_record][i][j]=cache_Kd*ratio_matrices[pos_of_last_record][i][j]+product_of_freq*ratio_matrices[pos_of_mat_nH][i][j-1];
                   
                    for(i=0;i<max_N_binding_rep;i++)
                        ratio_matrices[pos_next_record][i][0]=cache_Kd*ratio_matrices[pos_of_last_record][i][0];                    
                }
                else
                {
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;  //find last record     
                    for(i=0;i<max_N_binding_rep;i++)                  
                        for(j=0;j<max_N_binding_act;j++) 
                            ratio_matrices[pos_next_record][i][j]=cache_Kd*ratio_matrices[pos_of_last_record][i][j];
               
                    ratio_matrices[pos_next_record][0][1]+=product_of_freq;                    
                }
                break;

            case REPRESSOR: // a BS of repressors            
                if(m-BS_info[m].N_hindered!=0)
                {
                    /*find matrix(n-H)*/
                    pos_of_mat_nH=pos_next_record-BS_info[m].N_hindered-1; 
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1; 
                    /*Similar problem when pos_of_mat_nH=pos_next_record*/
                    for(i=1;i<max_N_binding_rep;i++)                  
                        for(j=0;j<max_N_binding_act;j++)
                            ratio_matrices[pos_next_record][i][j]=cache_Kd*ratio_matrices[pos_of_last_record][i][j]+product_of_freq*ratio_matrices[pos_of_mat_nH][i-1][j];
                   
                    for(j=0;j<max_N_binding_act;j++)
                        ratio_matrices[pos_next_record][0][j]=cache_Kd*ratio_matrices[pos_of_last_record][0][j];     
                }
                else
                {
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;  
                    for(i=0;i<max_N_binding_rep;i++)                  
                        for(j=0;j<max_N_binding_act;j++)
                            ratio_matrices[pos_next_record][i][j]=cache_Kd*ratio_matrices[pos_of_last_record][i][j];                                                   
                   
                    ratio_matrices[pos_next_record][1][0]+=product_of_freq;
                } 
                break;
        }
    }

    sum=0.0;
    for(i=0;i<max_N_binding_rep;i++)    
        for(j=0;j<max_N_binding_act;j++)     
            sum+=ratio_matrices[pos_next_record][i][j];
    
    temp=0.0;
    for(i=0;i<max_N_binding_rep;i++)  
        for(j=genotype->min_N_activator_to_transc[gene_id];j<max_N_binding_act;j++)   
            temp+=ratio_matrices[pos_next_record][i][j];
    state->P_A[gene_id]=(float)(temp/sum);
    
    temp=0.0;
    for(i=1;i<max_N_binding_rep;i++)  
        for(j=0;j<max_N_binding_act;j++)
            temp+=ratio_matrices[pos_next_record][i][j];  
    state->P_R[gene_id]=(float)(temp/sum);
   
	temp=0.0;
	for(j=genotype->min_N_activator_to_transc[gene_id];j<max_N_binding_act;j++)
		temp+=ratio_matrices[pos_next_record][0][j];
	state->P_A_no_R[gene_id]=(float)(temp / sum);
    
    temp=0.0;
    for(j=0;j<genotype->min_N_activator_to_transc[gene_id];j++)
        temp+=ratio_matrices[pos_next_record][0][j];
    state->P_NotA_no_R[gene_id]=(float)(temp/sum);    
}


/* 
 * update both the protein concentration and current cell size *
 * 
 */
static void update_protein_number_and_fitness( Genotype *genotype,
                                                CellState *state,                                   
                                                GillespieRates *rates,
                                                Environment *env,                                                
                                                float dt)
{
    int i,j,counter;
    float ct, ect, one_minus_ect;  
    float integrated_fitness;
    float N_output_molecules_bf_dt[MAX_OUTPUT_GENES];
    float N_output_molecules_aft_dt[MAX_OUTPUT_GENES];
  
    /* store the numbers of the effector proteins encoded by each copy of gene before updating*/   
    counter=0;
    for(i=0;i<genotype->n_output_genes;i++)
    {      
        N_output_molecules_bf_dt[counter]=state->gene_specific_protein_number[genotype->output_protein_ids[i]];
        counter++;        
    }
    /* update protein numbers*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {     
        ct=genotype->protein_decay_rate[i]*dt;
        ect = exp(-ct);
        if (fabs(ct)<EPSILON) one_minus_ect=ct;
        else one_minus_ect = 1.0-ect;      
        /* get the new protein concentration for this gene */
        state->gene_specific_protein_number[i]=ect*state->gene_specific_protein_number[i]+state->protein_synthesis_index[i]*one_minus_ect;        
    }    
    /* now, use protein_pool to pool gene specific protein number*/
    for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
    {
        state->protein_number[i]=0.0;        
        for(j=0;j<genotype->protein_pool[i][0][0];j++)
            state->protein_number[i]+=state->gene_specific_protein_number[genotype->protein_pool[i][1][j]];
    } 
    /*store the numbers of the effector proteins encoded by each copy of gene after updating*/
    counter=0;
    for(i=0;i<genotype->n_output_genes;i++)
    { 
        N_output_molecules_aft_dt[counter]=state->gene_specific_protein_number[genotype->output_protein_ids[i]];
        counter++;        
    }    
    /* now find out the protein numbers at end of dt interval and compute instantaneous and cumulative fitness */   
    calc_fitness(&integrated_fitness, 
                    genotype, 
                    state, 
                    env,
                    N_output_molecules_bf_dt, 
                    N_output_molecules_aft_dt,
                    dt);  
//    /* update cumulative fitness at the end of dt*/
//    state->cumulative_fitness += integrated_fitness;    
//    /* update the instantaneous fitness at the end of dt */
//    state->instantaneous_fitness = instantaneous_fitness;      
}

/*
 * 
 * Functions that handle each possible Gillespie event 
 *
 */
static int Gillespie_event_mRNA_decay(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id;
    float x;
    int mRNA_id;
    while(1)/*in case of numerical error*/    
    {           
        x=RngStream_RandU01(RS)*rates->total_mRNA_decay_rate;
        gene_id=N_SIGNAL_TF-1;
        /* loop through mRNA products, to choose the mRNA with the
        proportionally higher decay rate */
        while (gene_id < genotype->ngenes-1 && x > 0.0) 
        {
            gene_id++;
            x-= rates->mRNA_decay_rate[gene_id];
        }
        /*rarely, numerical error picks up a gene that has no mRNA at all*/
        if((state->mRNA_aft_transl_delay_num[gene_id]+state->mRNA_under_transl_delay_num[gene_id])>=1)
            break;
    }    
    /* assume mRNAs are equally likely to be degraded */
    x = RngStream_RandInt(RS,1,state->mRNA_aft_transl_delay_num[gene_id] + state->mRNA_under_transl_delay_num[gene_id]);    
    /* decay mRNA in cytoplasm */
    if (x <= state->mRNA_aft_transl_delay_num[gene_id])
    {
        /* remove the mRNA from the cytoplasm count */
        (state->mRNA_aft_transl_delay_num[gene_id])--;  
        /*update protein synthesis rate*/
        state->protein_synthesis_index[gene_id] = (float)state->mRNA_aft_transl_delay_num[gene_id]*genotype->protein_syn_rate[gene_id]/genotype->protein_decay_rate[gene_id];
//        if(genotype->which_protein[gene_id]==genotype->nproteins-1)
//        	return DO_NOTHING;
//    	else // an mRNA of transcription factor is degraded, which can cause fluctuation in transcription factor concentrations.
        	return gene_id;    
    } 
    else 
    {
        /* decay mRNA in process of translation initialization */       
        mRNA_id = RngStream_RandInt(RS,0,state->mRNA_under_transl_delay_num[gene_id]-1);
        /* delete this fixed event: this mRNA will never be translated */     
        delete_fixed_event(gene_id, mRNA_id, &(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));       
        /* remove the mRNA from the count */
        (state->mRNA_under_transl_delay_num[gene_id])--; 
        return DO_NOTHING;
    }
}

static void Gillespie_event_repressed_to_intermediate(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x= RngStream_RandU01(RS)*rates->total_repressed_to_intermediate_rate;
        gene_id=N_SIGNAL_TF-1;
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->repressed_to_intermediate_rate[gene_id];
        }
        if(rates->repressed_to_intermediate_rate[gene_id]>0.0)
            break;
    }
    /* set state */
    state->transcriptional_state[gene_id]=INTERMEDIATE;
}

static void Gillespie_event_intermediate_to_repressed(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id; 
    float x; 
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->total_intermediate_to_repressed_rate;
        gene_id=N_SIGNAL_TF-1;
        /* choose a particular gene copy to change state */
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->intermediate_to_repressed_rate[gene_id];
        }
        if(rates->intermediate_to_repressed_rate[gene_id]>0.0)
            break;
    }
    /* set state */
    state->transcriptional_state[gene_id]=REPRESSED;
}

static void Gillespie_event_intermediate_to_active(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    float x; 
    int gene_id;
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->total_intermediate_to_active_rate;
        gene_id=N_SIGNAL_TF-1;
        /* choose a particular gene copy to change state */
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->intermediate_to_active_rate[gene_id];
        }
        if(rates->intermediate_to_active_rate[gene_id]>0.0)
            break;
    }
    /* set state */
    state->transcriptional_state[gene_id] =ACTIVE;  
}

static void Gillespie_event_active_to_intermediate(Genotype *genotype, CellState *state,GillespieRates *rates, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x=RngStream_RandU01(RS)*rates->total_active_to_intermediate_rate;   
        gene_id=N_SIGNAL_TF-1;       
        while(gene_id < genotype->ngenes-1 && x>0.0) 
        {
            gene_id++;       
            x -= rates->active_to_intermediate_rate[gene_id];
        }
        if(rates->active_to_intermediate_rate[gene_id]>0.0)
            break;
    }
    state->transcriptional_state[gene_id]=INTERMEDIATE;   
}

static void Gillespie_event_transcription_init(GillespieRates *rates, CellState *state, Genotype *genotype, float dt, RngStream RS)
{
    int gene_id;  
    int x;
    float candidate_t;
    int concurrent;    
    gene_id=N_SIGNAL_TF-1;    
    x=RngStream_RandInt(RS,1,rates->total_N_gene_transcript_initiated);
    while(gene_id<genotype->ngenes-1 && x>0)
    {
        gene_id++;
        x-=rates->transcript_initiation_state[gene_id];
    }     
    /* now that transcription of gene has been initiated, 
     * we add the timepoint at which the transcription ends, 
     * which is dt+time-of-transcription from now */
    candidate_t=state->t+dt+(float)genotype->locus_length[gene_id]/TRANSCRIPTION_ELONGATION_RATE+TRANSCRIPTION_TERMINATION_TIME;
    concurrent=check_concurrence(state, candidate_t);
    while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
    {
        candidate_t+=TIME_OFFSET;
        concurrent=check_concurrence(state, candidate_t);        
    }    
    add_fixed_event(gene_id, candidate_t,&(state->mRNA_transcr_time_end_head), &(state->mRNA_transcr_time_end_tail));
    /* increase the number mRNAs being transcribed */
    (state->mRNA_under_transc_num[gene_id])++;                      
}
/*
 * END
 * Functions that handle each possible Gillespie event 
 */

/* do a fixed event that occurs in current t->dt window */
static int do_fixed_event(Genotype *genotype, 
                            CellState *state, 
                            GillespieRates *rates, 
                            Environment *env,
                            Phenotype *timecourse,
                            float *dt,                           
                            int event)
{     
    int i, j, protein_id, return_value;
    return_value=DO_NOTHING;
    switch (event) 
    {
        case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, state, rates, genotype, env); 
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, env, dt);
            state->cell_activated=1;
            break;
        case 3:     /* turn signal off*/ 
            *dt = state->signal_off_head->time - state->t;     
            update_protein_number_and_fitness(genotype, state, rates, env, *dt); 
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));
#if SELECT_SENSITIVITY_AND_PRECISION
			state->precision[0] = state->protein_number[genotype->nproteins - 1]; //reset running precision
#endif
//            if(env->fixed_effector_effect)
//                state->effect_of_effector=env->initial_effect_of_effector;
//            else
//                state->effect_of_effector='d';
//            state->protein_number[N_SIGNAL_TF-1]=env->signal_off_strength;         
            break;
        case 4:     /*turn signal on*/
            *dt = state->signal_on_head->time - state->t;   
            update_protein_number_and_fitness(genotype, state, rates, env, *dt);  
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
            state->protein_number[N_SIGNAL_TF-1]=env->signal_on_strength;
            state->gene_specific_protein_number[N_SIGNAL_TF-1]=env->signal_on_strength;
            state->recording_basal_response=0;
            state->found_gradient=1;
#if POOL_EFFECTORS
            state->threshold_response_to_bias[0]=(env->minimal_peak_response>state->basal_response[0]*env->response_amplification)?env->minimal_peak_response:state->basal_response[0]*env->response_amplification;
#endif 
#if SELECT_SENSITIVITY_AND_PRECISION
			state->precision[1] = state->precision[0]; //record the max precision to signal change 1
			state->sensitivity[1] = state->sensitivity[0]; //record the max sensitivity to signal change 1
			state->sensitivity[0] = 0.0;	// reset running sensitivity
#endif
//            if(env->fixed_effector_effect)                               
//                state->effect_of_effector=env->initial_effect_of_effector;            
//            else                 
//                state->effect_of_effector='b';            
            break;	
        case 5: /* finishing burn-in developmental simulation*/
            *dt=env->duration_of_burn_in_growth_rate-state->t;     
            update_protein_number_and_fitness(genotype, state, rates, env, *dt);
            for(i=0;i<MAX_OUTPUT_PROTEINS;i++)
                state->cumulative_fitness_after_burn_in[i]=state->cumulative_fitness[i];           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6: /* mandatorily updating Pact and Prep*/
            *dt=state->t_to_update_probability_of_binding-state->t;
            update_protein_number_and_fitness(genotype, state, rates, env, *dt);          
            break;
        case 7: /* update signal strength */
            *dt=state->change_signal_strength_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, env, *dt);
            state->protein_number[N_SIGNAL_TF-1]=env->external_signal[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            break;     
        case 8: /* record expression levels*/
            *dt=state->sampling_point_end_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, env, *dt);
            delete_fixed_event_from_head(&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
#if PHENOTYPE
            for(i=0;i<genotype->N_node_families;i++)
            {                
                for(j=0;j<genotype->node_family_pool[i][0][0];j++)
                    timecourse->protein_concentration[i*timecourse->total_time_points+timecourse->timepoint]+=state->gene_specific_protein_number[genotype->node_family_pool[i][1][j]];
            }
//            for(i=0;i<genotype->ngenes;i++)
//                timecourse->gene_specific_concentration[i*timecourse->total_time_points+timecourse->timepoint]=state->gene_specific_protein_number[i];                
//            timecourse->instantaneous_fitness[timecourse->timepoint]=state->instantaneous_fitness;
            timecourse->timepoint++;   
#else
#if CHEMOTAXIS
            state->recording_basal_response=1;
#else
            state->sampled_response[state->N_samples]=0.0;
            for(i=0;i<genotype->n_output_genes;i++)
                state->sampled_response[state->N_samples]+=state->gene_specific_protein_number[genotype->output_protein_ids[i]];
            state->N_samples++;   
#endif
#endif
            break;
    }  
    return return_value;
}

/*
 *Gillespie events
 */
static int do_Gillespie_event(Genotype *genotype,
                                CellState *state,
                                GillespieRates *rates,
                                float dt,                        
                                RngStream RS)
{
    float x; 
    int return_value;
    return_value=DO_NOTHING;
    while(1)
    {    
        x=RngStream_RandU01(RS)*(rates->total_Gillespie_rate);          
        if (x <= rates->total_mRNA_decay_rate)  /*STOCHASTIC EVENT: an mRNA decay event */
        {    
            return_value=Gillespie_event_mRNA_decay(rates, state, genotype, RS);                
            break;
        } 
        else 
        {               
            x -= rates->total_mRNA_decay_rate;             
            if (x <= rates->total_active_to_intermediate_rate) /* STOCHASTIC EVENT: active transcriptional state to intermediate*/
            {                   
                Gillespie_event_active_to_intermediate(genotype, state, rates, RS);
                break;
            } 
            else 
            {  
                x -= rates->total_active_to_intermediate_rate;  
                if (x <= rates->total_repressed_to_intermediate_rate)  /* repressed to intermediate*/
                {    
                    Gillespie_event_repressed_to_intermediate(rates, state, genotype, RS);
                    break;
                } 
                else 
                { 
                    x-= rates->total_repressed_to_intermediate_rate;                  
                    if (x <= rates->total_intermediate_to_repressed_rate)/* intermediate to repressed */ 
                    {   
                        Gillespie_event_intermediate_to_repressed(rates, state, genotype, RS);
                        break;
                    } 
                    else 
                    {
                        x -= rates->total_intermediate_to_repressed_rate;                        
                        if (x <= rates->total_intermediate_to_active_rate)/* intermediate to active*/
                        {  
                            Gillespie_event_intermediate_to_active(rates, state, genotype, RS); 
                            break;
                        } 
                        else 
                        {
                            x -= rates->total_intermediate_to_active_rate;                            
                            if (x <= rates->total_N_gene_transcript_initiated * TRANSCRIPTINIT) /* STOCHASTIC EVENT: transcription initiation */
                            {   
                                Gillespie_event_transcription_init(rates, state, genotype, dt, RS);
                                break;
                            }   
                        }
                    }
                }
            }       
        }
    }
    return return_value;
}

/* 
 * check to see if a fixed event ends within dt
 *
 * returns:
 *  0 if there is no fixed event occuring before time t+dt
 *  1 a transcription event happens first before time t+dt
 *  2 a translation initiation event starts 
 *  3 time to turn on the environmental signal 
 *  4 time to turn off the signal
 *  5 burn-in developmental simulation ends
 *  6 time to update TF binding probabilities mandatorily
 *  7 time to change signal strength (used when an external signal profile is provided)
 *  8 time to sample expression levels (used under PHENOTYPE mode)
 */
 
static int does_fixed_event_end(CellState *state, float t) 
{
    int retval=0;
    float t1;
    float t2;
    float t3;
    float t4;
    float t5;  
    float t6;    
    float t7;
    float t8;
    t1 = state->mRNA_transcr_time_end_head ? state->mRNA_transcr_time_end_head->time : TIME_INFINITY;
    t2 = state->mRNA_transl_init_time_end_head ? state->mRNA_transl_init_time_end_head->time : TIME_INFINITY;
    t3 = state->signal_off_head? state->signal_off_head->time : TIME_INFINITY;
    t4 = state->signal_on_head ? state->signal_on_head->time : TIME_INFINITY;
    t5 = state->burn_in_growth_rate_head ? state->burn_in_growth_rate_head->time : TIME_INFINITY;
    t6 = state->t_to_update_probability_of_binding;
    t7 = state->change_signal_strength_head ? state->change_signal_strength_head->time : TIME_INFINITY;
    t8 = state->sampling_point_end_head?state->sampling_point_end_head->time : TIME_INFINITY;
    if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5) &&(t1<=t6) && (t1<=t7) && (t1<=t8))
    {
        retval = 1;	
    }
    else if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5)&&(t2<=t6)&&(t2<=t7)&& (t2<=t8))
    { 
        retval = 2;
    }  
    else if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5)&&(t3<=t6) &&(t3<=t7)&& (t3<=t8)) 
    {
        retval = 3;
    }
    else if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5)&&(t4<=t6)&&(t4<=t7)&& (t4<=t8)) 
    {
        retval = 4;
    }               
    else if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4)&&(t5<=t6)&&(t5<=t7)&& (t5<=t8))
    {
        retval = 5;
    }             
    else if((t6 <= t1) && (t6 <= t) && (t6 <= t2) && (t6 <= t3) && (t6<=t4)&&(t6<=t5)&&(t6<=t7)&& (t6<=t8))
    {
        retval=6;
    }
    else if((t7 <= t1) && (t7 <= t) && (t7 <= t2) && (t7 <= t3) && (t7<=t4)&&(t7<=t5)&&(t7<=t6)&& (t7<=t8))
    {
        retval=7;
    }
    else if((t8 <= t1) && (t8 <= t) && (t8 <= t2) && (t8 <= t3) && (t8<=t4)&&(t8<=t5)&&(t8<=t6)&& (t8<=t7))
    {
        retval = 8;
    }
    else
    {
        retval=0;
    }
    return retval;
}

/*
 * end transcription: update the mRNAs ready for translation initiation
 * etc. accordingly and delete the event from the queue
 */
static void fixed_event_end_transcription( float *dt,                                     
                                            CellState *state,
                                            GillespieRates *rates,
                                            Genotype *genotype,
                                            Environment *env)
{
    int gene_id;
    float concurrent;
    float endtime;
    /* recompute the delta-t based on difference between now and the time of transcription end */
    *dt = state->mRNA_transcr_time_end_head->time - state->t;   
    /* update fitness and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, env, *dt);
    /* get the gene which is ending transcription */
    gene_id = state->mRNA_transcr_time_end_head->event_id;    
    /* increase number of mRNAs that are initializing translation*/
    (state->mRNA_under_transl_delay_num[gene_id])++;
    /* decrease the number of mRNAs undergoing transcription */
    (state->mRNA_under_transc_num[gene_id])--;
    /* delete the fixed even which has just occurred */
    if(state->mRNA_transcr_time_end_head==NULL)
    {
#if MAKE_LOG
        LOG("No transcription is going on\n");
#endif
        exit(-3);
    }
    delete_fixed_event_from_head(&(state->mRNA_transcr_time_end_head), &(state->mRNA_transcr_time_end_tail));   
    /*add transcription initialization event*/ 
    endtime=state->t+*dt+(float)genotype->locus_length[gene_id]/TRANSLATION_ELONGATION_RATE+TRANSLATION_INITIATION_TIME;

    concurrent=check_concurrence(state, endtime);
    while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
    {
        endtime+=TIME_OFFSET;
        concurrent=check_concurrence(state, endtime);        
    }  
    /*add to translation initiation event*/
    add_fixed_event(gene_id, endtime, &(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));
}

/*
 * end translation initiation
 */
static int fixed_event_end_translation_init(Genotype *genotype, 
                                            CellState *state, 
                                            GillespieRates *rates, 
                                            Environment *env,
                                            float *dt)
{
    int gene_id;    
    /* calc the remaining time till translation initiation ends */
    *dt = state->mRNA_transl_init_time_end_head->time - state->t;         
    /* update fitness and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, env, *dt);
    /* get identity of gene that has just finished translating */
    gene_id=state->mRNA_transl_init_time_end_head->event_id; 
    /* there is one less mRNA that is initializing translation */
    (state->mRNA_under_transl_delay_num[gene_id])--;  
    /* delete the event that just happened */
    if(state->mRNA_transl_init_time_end_head==NULL)
    {
#if MAKE_LOG
        LOG("No translation initiation going on!\n");
#endif
        exit(-3);
    }
    delete_fixed_event_from_head(&(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));    
    /* there is one more mRNA that produces protein */
    (state->mRNA_aft_transl_delay_num[gene_id])++;   
    /* update protein synthesis rate*/
    state->protein_synthesis_index[gene_id]= (float)state->mRNA_aft_transl_delay_num[gene_id]*genotype->protein_syn_rate[gene_id]/genotype->protein_decay_rate[gene_id];
    
//    if(genotype->which_protein[gene_id]==genotype->nproteins-1)//if the mRNA encodes a non-sensor TF, there could be a huge change in TF concentration
//        return DO_NOTHING;
//    else
        return gene_id; //so update Pact and reset updating interval to MIN_INTERVAL_TO_UPDATE_PACT
}

static void calc_leaping_interval(Genotype *genotype, CellState *state, float *minimal_interval, float t_unreachable, int which_gene)
{
    int protein_id,j;
    float dt;
    float N_proteins_cause_change,Kd,N_at_end_of_simulation;
    float overall_rate;
    float P_binding;
    float t_remaining;    
    int gene_ids[genotype->ngenes]; 
    float ct, ect, one_minus_ect;
    
    t_remaining=t_unreachable-state->t;
    dt=t_unreachable;  
 	protein_id=genotype->which_protein[which_gene];
    Kd=KD2APP_KD*genotype->Kd[protein_id];
    P_binding=state->protein_number[protein_id]/(state->protein_number[protein_id]+Kd);
    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
        gene_ids[j]=genotype->protein_pool[protein_id][1][j];

    /*determine whether the protein tends to increase or decrease concentration*/
    overall_rate=0.0;
    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
        overall_rate+=(state->protein_synthesis_index[gene_ids[j]]-state->gene_specific_protein_number[gene_ids[j]])*genotype->protein_decay_rate[gene_ids[j]];
    
    if(overall_rate>0.0) //tend to increase
    {            
        if(P_binding>=(1.0-MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING)) //concentration already very high
            dt=t_unreachable;
        else
        {
            /* check if much change is possible within the duration of simulation*/
            N_proteins_cause_change=Kd*(P_binding+MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING)/(1.0-P_binding-MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING); 
            /* calc N_protein at the end of simulation*/
            N_at_end_of_simulation=0.0;
            for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
            {                  
                ct=genotype->protein_decay_rate[gene_ids[j]]*t_remaining;
                ect = exp(-ct);
                if (fabs(ct)<EPSILON) one_minus_ect=ct;
                else one_minus_ect = 1.0-ect;   
                N_at_end_of_simulation+=ect*state->gene_specific_protein_number[gene_ids[j]]+state->protein_synthesis_index[gene_ids[j]]*one_minus_ect;        
            }               
            if(N_at_end_of_simulation<N_proteins_cause_change)
               dt=t_unreachable; 
            else /*We need to solve an equation*/
                dt=calc_tprime(genotype,state,state->gene_specific_protein_number,t_remaining,N_proteins_cause_change,gene_ids,0,genotype->protein_pool[protein_id][0][0]); 
        }
        *minimal_interval=(*minimal_interval<dt)?*minimal_interval:dt; 
    }
    else //tends to decrease
    {
        if(P_binding<=MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING) //concentration already very low
            dt=t_unreachable;
        else
        {
            /* first check if much change is possible within the duration of simulation*/
            N_proteins_cause_change=Kd*(P_binding-MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING)/(1.0-P_binding+MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING); 
            /* calc N_protein at the end of simulation*/
            N_at_end_of_simulation=0.0;
            for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
            {  
                ct=genotype->protein_decay_rate[gene_ids[j]]*t_remaining;
                ect = exp(-ct);
                if (fabs(ct)<EPSILON) one_minus_ect=ct;
                else one_minus_ect = 1.0-ect;   
                N_at_end_of_simulation+=ect*state->gene_specific_protein_number[gene_ids[j]]+state->protein_synthesis_index[gene_ids[j]]*one_minus_ect;        
            }               
            if(N_at_end_of_simulation>N_proteins_cause_change)
               dt=t_unreachable; 
            else
               dt=calc_tprime(genotype,state,state->gene_specific_protein_number,t_remaining,N_proteins_cause_change,gene_ids,0,genotype->protein_pool[protein_id][0][0]);
        }
        *minimal_interval=(*minimal_interval<dt)?*minimal_interval:dt;  
    }
}



