/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
static const float NS_Kd=1.0e-5;
static const float KD2APP_KD=1.8e10;
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
static float calc_tprime(Genotype*, CellState*, float*, float, float, int);

static float calc_integral(Genotype *, CellState *, float *, float, float);

static void calc_fx_dfx(float, int, float, float*, float*, float*, float*, float*);

static void calc_leaping_interval(Genotype*, CellState*, float *, float, int);

static void calc_TF_dist_from_all_BS(Genotype *, CellState*, int);

static int Gillespie_event_mRNA_decay(GillespieRates *, CellState *, Genotype *, RngStream);

static void Gillespie_event_repressed_to_intermediate(GillespieRates *, CellState *, Genotype *, RngStream);

static void Gillespie_event_intermediate_to_repressed(GillespieRates *, CellState *, Genotype *, RngStream); 

static void Gillespie_event_intermediate_to_active(GillespieRates *, CellState *, Genotype *, RngStream);

static void Gillespie_event_active_to_intermediate(Genotype *, CellState *, GillespieRates *, RngStream);

static void Gillespie_event_transcription_init(GillespieRates *, CellState *, Genotype *, float, RngStream);

static int fixed_event_end_translation_init(Genotype *, CellState *, GillespieRates *, float *);

static void fixed_event_end_transcription(float *, CellState *, GillespieRates *, Genotype *);

static int does_fixed_event_end(CellState*, float);

static int do_fixed_event(Genotype *, CellState *, GillespieRates *, Test *, Phenotype *, float *, int);

static float calc_fitness(float *, Genotype *, CellState *, float*, float);

static void update_protein_number_and_fitness(Genotype *, CellState *, GillespieRates *, float);

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
                        Test *test,
                        int init_mRNA_number[MAX_GENES],
                        float init_protein_number[MAX_PROTEINS])
{
    int i, j;
    state->t=0.0;
    state->cumulative_fitness = 0.0;     
    state->cumulative_fitness_after_burn_in = 0.0;   
    state->instantaneous_fitness = 0.0; 
    state->effect_of_effector=test->initial_effect_of_effector;
    state->error=0; //no error in simulation
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
        state->protein_synthesis_index[i]=(float)state->mRNA_aft_transl_delay_num[i]*genotype->translation_rate[i]/genotype->protein_decay_rate[i];
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
    /*mark when to start calculating average growth rate*/
    if(test->duration_of_burn_in_growth_rate!=0.0)
        add_fixed_event(-1,test->duration_of_burn_in_growth_rate,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail)); 
    else
        add_fixed_event(-1,(float)TIME_INFINITY,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));                
    /*plot protein concentration and fitness vs time*/
    #if RUN_PLOTTING 
        float t;
        int N_data_points; 
        t=0.0+TIME_OFFSET;
        N_data_points=(int)test->t_development;
        for(i=0;i<N_data_points;i++)
        {
            add_fixed_event(-1,t,&(state->sampling_point_end_head),&(state->sampling_point_end_tail)); //get a timepoint each minute
            t+=1.0;            
        } 
    #endif    
}

/*
 * Calculate the rates of all Gillespie events
 */
void calc_all_rates(Genotype *genotype,
                    CellState *state,
                    GillespieRates *rates,
                    Test *test,
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
        #if NO_REGULATION //if we manually turn off the expression of non-sensor TFs
            if(genotype->protein_identity[genotype->which_protein[i]]==NON_TF) // we only calculate the binding configurations for effector genes
            {
                cluster_id=genotype->which_cluster[i];        
                if(genotype->cisreg_cluster[cluster_id][0]!=i)  /*if this gene does not have a unique cis-reg sequence*/
                {                
                    state->P_A[i]=state->P_A[genotype->cisreg_cluster[cluster_id][0]]; /* copy TF distribution from elsewhere*/
                    state->P_R[i]=state->P_R[genotype->cisreg_cluster[cluster_id][0]];
                    state->P_A_no_R[i]=state->P_A_no_R[genotype->cisreg_cluster[cluster_id][0]];
                    state->P_NotA_no_R[i]=state->P_NotA_no_R[genotype->cisreg_cluster[cluster_id][0]];
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
            }
            else
            {
               state->P_A[i]=0.0; 
               state->P_R[i]=0.0;
			   state->P_A_no_R[i] = 0.0;
               state->P_NotA_no_R[i]=0.0;
            }
        #else
            cluster_id=genotype->which_cluster[i];        
            if(genotype->cisreg_cluster[cluster_id][0]!=i)  /*if this gene does not have a unique cis-reg sequence*/
            {                
                state->P_A[i]=state->P_A[genotype->cisreg_cluster[cluster_id][0]]; /* copy TF distribution from elsewhere*/
                state->P_R[i]=state->P_R[genotype->cisreg_cluster[cluster_id][0]];
                state->P_A_no_R[i]=state->P_A_no_R[genotype->cisreg_cluster[cluster_id][0]];
                state->P_NotA_no_R[i]=state->P_NotA_no_R[genotype->cisreg_cluster[cluster_id][0]];
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
        #endif
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
        cluster_id=1;    
        while(genotype->cisreg_cluster[cluster_id][0]!=NA) //check if Pact changes too much
        {
            gene_id=genotype->cisreg_cluster[cluster_id][0];
            diff_PA=fabs(state->P_A[gene_id]-state->last_P_A[gene_id]);            
            diff_PAnoR=fabs(state->P_A_no_R[gene_id]-state->last_P_A_no_R[gene_id]); 
            diff_PnotAnoR=fabs(state->P_NotA_no_R[gene_id]-state->last_P_NotA_no_R[gene_id]);  
            diff_PR=fabs(state->P_R[gene_id]-state->last_P_R[gene_id]);
            diff_max=(diff_max>diff_PA)?diff_max:diff_PA;
            diff_max=(diff_max>diff_PR)?diff_max:diff_PR;
            diff_max=(diff_max>diff_PAnoR)?diff_max:diff_PAnoR;
            diff_max=(diff_max>diff_PnotAnoR)?diff_max:diff_PnotAnoR;
            cluster_id++;
        }
        if(diff_max<EPSILON)
            interval_to_update_probability_of_binding=DEFAULT_UPDATE_INTERVAL;
        else
            interval_to_update_probability_of_binding=MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING/diff_max*(state->t-state->last_event_t);
        if(UPDATE_WHAT!=DO_NOTHING)          
            calc_leaping_interval(genotype,state,&interval_to_update_probability_of_binding,test->t_development,UPDATE_WHAT);  
    
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
                        Test *test,   
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
    fixed_time = (state->t+dt<test->t_development)?(state->t+dt):test->t_development;
    event = does_fixed_event_end(state, fixed_time);
    while(event!=0)
    {           
        /*after doing fixed event, return a flag to indicate whether mandatorily update Pact or Prep*/
        UPDATE_WHAT=do_fixed_event(genotype, state, rates, test, timecourse, &dt, event);        
        if(state->error==1) 
            return;
        state->t += dt;  /* advance time by the dt */       
        x -= dt*rates->total_Gillespie_rate;  /* we've been running with rates->total_Gillespie_rate for dt, so substract it from x*/        
        calc_all_rates(genotype, state, rates, test, UPDATE_WHAT);  /* update rates->total_Gillespie_rate and re-compute a new dt */      
        dt = x/rates->total_Gillespie_rate;
        /*deal with rounding error*/
        if(dt<0.0)
        {  	
#if CAUTIOUS // this rounding error can happen very often, therefore the error_log can be huge            
            LOG("rounding error in dt at mut_step %d\n",mut_step);    
#endif 
            dt=TIME_OFFSET; 
        }
        fixed_time = (state->t+dt<test->t_development)?(state->t+dt):test->t_development;
        /* check to see there aren't more fixed events to do */
        event = does_fixed_event_end(state, fixed_time);                                    
    } 
  /* no remaining fixed events to do in dt, now do stochastic events */  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (state->t+dt < test->t_development)
    {
        update_protein_number_and_fitness(genotype, state, rates, dt);        
        if(state->error==1)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, RS);
        if(state->error==1)
            return;  
        /* Gillespie step: advance time to next event at dt */
        state->t += dt;
        calc_all_rates(genotype,state,rates,test,UPDATE_WHAT);        
    } 
    else 
    { 
        /* do remaining dt */
        dt = test->t_development - state->t;
        /* final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, dt); 
        if(state->error==1)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        state->t = test->t_development;
    }
}



/*****************************************************************************
 * 
 *                             Private functions
 * 
 *****************************************************************************/

/*
 * compute t' factor used in the integration of growth rate
 * t' is the time the effector protein increases or decreases to a given amount
 */
static float calc_tprime(Genotype *genotype, CellState* state, float *number_of_selection_protein_bf_dt, float dt, float given_amount, int protein_id) 
{
    int n_copies;
    int i;          
    n_copies=genotype->protein_pool[protein_id][0][0];
    float protein_synthesis_rate[n_copies],protein_decay_rate[n_copies];
    for(i=0;i<n_copies;i++)
    {
        protein_decay_rate[i]=genotype->protein_decay_rate[genotype->protein_pool[protein_id][1][i]];
        protein_synthesis_rate[i]=state->protein_synthesis_index[genotype->protein_pool[protein_id][1][i]]*protein_decay_rate[i];
    }       
    return rtsafe(&calc_fx_dfx, n_copies, given_amount, number_of_selection_protein_bf_dt, protein_synthesis_rate, protein_decay_rate, 0.0, dt); 
}

/*
 * calculate f(x)-Pp_s and f'(x),
 * f(x) is the number of effector protein molecules at time x
 */
static void calc_fx_dfx(float x, int n_copies, float given_amount, float *intial_protein_number, float *protein_synthesis_rate, float *protein_decay_rate,float *fx, float *dfx)
{
    int i;    
    *fx=0;
    *dfx=0;    
    for(i=0;i<n_copies;i++)
    {
        *fx+=(intial_protein_number[i]-protein_synthesis_rate[i]/protein_decay_rate[i])*exp(-protein_decay_rate[i]*x)+protein_synthesis_rate[i]/protein_decay_rate[i];
        *dfx+=(protein_synthesis_rate[i]-protein_decay_rate[i]*intial_protein_number[i])*exp(-protein_decay_rate[i]*x);
    }
    *fx-=given_amount;    
}

/*
 * calculate F(delta_t)/Ne_sat. F(x) is the integral of f(x) over delta_t.
 * f(x) is the number of effector protein molecules at time x
 */
static float calc_integral(Genotype *genotype, CellState *state, float *initial_protein_number, float dt, float saturate_protein_number)
{
    int i,n_copies,gene_ids[MAX_PROTEINS];
    float integral=0.0,ect_minus_one;    
   
    n_copies=genotype->protein_pool[genotype->nproteins-1][0][0];
    for(i=0;i<n_copies;i++)
        gene_ids[i]=genotype->protein_pool[genotype->nproteins-1][1][i];
        
    for(i=0;i<n_copies;i++)
    {
        ect_minus_one=exp(-genotype->protein_decay_rate[gene_ids[i]]*dt)-1.0;    
        integral+=(state->protein_synthesis_index[gene_ids[i]]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]-
                initial_protein_number[i]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]+ 
                state->protein_synthesis_index[gene_ids[i]]*dt);
    }
    return integral/saturate_protein_number;
}

/*
 * return the instantaneous growth rate given the current cell state and environment,
 * also return the integrated growth rate as a pointer
 */
float calc_fitness(float *integrated_fitness,
                            Genotype *genotype,
                            CellState *state,
                            float* number_of_selection_protein_bf_dt,
                            float dt)
{
    int i;
    float instantaneous_fitness=0.0;  /* this is returned from the function */
    float total_translation_rate = 0.0;    
    float dt_prime;   
    float cost_of_expression;     
    float Ne_next=state->protein_number[genotype->nproteins-1];   
    float Ne=0.0;
    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++) 
       Ne+=number_of_selection_protein_bf_dt[i];
                      
    /* compute the total cost of translation across all genes  */
#if NO_REGULATION_COST
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {    
        if(genotype->which_protein[i]==genotype->nproteins-1)           
            total_translation_rate += genotype->translation_rate[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation_rate[i]*(float)state->mRNA_under_transl_delay_num[i];
    } 
#else
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {     
        total_translation_rate += (genotype->translation_rate[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation_rate[i]*(float)state->mRNA_under_transl_delay_num[i])*(float)genotype->locus_length[i]/pow(10.0,MEAN_GENE_LENGTH);
    } 
#endif  
    cost_of_expression=total_translation_rate*c_transl;

    switch (state->effect_of_effector)
    {
        case 'b': /* effector is beneficial!*/
            if(Ne>Ne_next)//decrease in effector protein
            {
                if(Ne_next>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =dt*(bmax-cost_of_expression);
                }
                else if(Ne<=Ne_saturate) // not enough effector throughout
                {
                    *integrated_fitness = bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                            -cost_of_expression*dt;
                }
                else // bf dt_prime, the benefit saturates
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*dt_prime+bmax*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)-
                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate))-
                                                    cost_of_expression*dt;                    
                }                    
            }
            else // increase in effector protein
            {
                if(Ne>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =dt*(bmax-cost_of_expression);
                }   
                else if(Ne_next<=Ne_saturate)// not enough effector throughout
                {
                    *integrated_fitness = bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                                    -cost_of_expression*dt;
                }
                else //Aft dt_prime, the benefit saturates
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*(dt-dt_prime)+bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate)-
                                                    cost_of_expression*dt;
                }                
            } 
            /* compute instantaneous growth rate at t */
            if (Ne_next < Ne_saturate)
                instantaneous_fitness = bmax*Ne_next/Ne_saturate;
            else
                instantaneous_fitness = bmax;
            break;
    
        case 'd': /* effector is deleterious! */      
#if !NO_PENALTY
            if(Ne>Ne_next)//decrease in effector protein
            {
                if(Ne_next>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =0.0-dt*cost_of_expression;
                }
                else if(Ne<=Ne_saturate) // not enough effector throughout
                {
                    *integrated_fitness = bmax*dt-bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                                    -cost_of_expression*dt;
                }
                else // aft dt_prime, the benefit becomes positive
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*(dt-dt_prime)-bmax*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)-
                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate))-
                                                    cost_of_expression*dt;                    
                }                    
            }
            else // increase in effector protein
            {
                if(Ne>=Ne_saturate) //too many effector throughout
                {
                    *integrated_fitness =0.0-dt*cost_of_expression;
                }   
                else if(Ne_next<=Ne_saturate)// not enough effector throughout
                {
                    *integrated_fitness = bmax*dt-bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                                    -cost_of_expression*dt;
                }
                else //Aft dt_prime, the benefit becomes zero
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate,genotype->nproteins-1); 
                    *integrated_fitness = bmax*dt_prime-bmax*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate)-
                                                    cost_of_expression*dt;
                }                
            } 
            if(Ne_next<Ne_saturate)
                instantaneous_fitness = bmax - bmax/Ne_saturate*Ne_next;
            else
                instantaneous_fitness = 0.0;            
#else
            *integrated_fitness=(bmax-cost_of_expression)*dt;
            instantaneous_fitness=bmax;            
#endif
            break; 
    } 
    /* and instantaneous integrated rate */
    instantaneous_fitness -= cost_of_expression;
    /* return the instantaneous growth rate */
    return instantaneous_fitness;
}

/*Calculate probability of binding configurations*/
void calc_TF_dist_from_all_BS(Genotype *genotype, CellState *state, int gene_id) 
{
    int max_N_binding_act=genotype->max_unhindered_sites[gene_id][1]+1; //Binding configurations can contain at most x activators, plus 1 type of configurations that don't have activators at all. 
    int max_N_binding_rep=genotype->max_unhindered_sites[gene_id][2]+1; //Binding configurations can contain at most y repressors, plus 1 type of configurations that don't have repressors at all. 
    double ratio_matrices[genotype->binding_sites_num[gene_id]+1][max_N_binding_rep][max_N_binding_act]; 
    double sum;    
    double product_of_freq; 
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
    if(genotype->protein_identity[BS_info[0].tf_id]==1) // if a activator binds to BS 0   
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
        /*Check whether m is a site of activator or repressor*/
        switch(genotype->protein_identity[BS_info[m].tf_id])
        {
            case ACTIVATOR: // a BS of activators              
               if(m-BS_info[m].N_hindered!=0)//if binding of m does not block all of the BS evaluated before
                {      
                    /*find matrix(n-H)*/
                    pos_of_mat_nH=pos_next_record-BS_info[m].N_hindered-1; 
                   /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;
                    for(i=0;i<max_N_binding_rep;i++)
                    {                
                        /*Note that it is possible for pos_of_mat_nH=pos_next_record, 
                         *which means we will be reading and writing to the same
                         *matrix. Updating the xth column of matrix(n+1) uses
                         *the (x-1)th column of matrix(n-H). To avoiding changing the values
                         *of matrix(n-H) before using it to update matrix(n+1), we
                         *need to update matrix(n+1) from its last column.*/
                        for(j=max_N_binding_act-1;j>0;j--) 
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j]+product_of_freq*ratio_matrices[pos_of_mat_nH][i][j-1];
                    }  
                    for(i=0;i<max_N_binding_rep;i++)
                        ratio_matrices[pos_next_record][i][0]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][0];                    
                }
                else
                {
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;  //find last record     
                    for(i=0;i<max_N_binding_rep;i++)
                    {
                        for(j=0;j<max_N_binding_act;j++) 
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j];
                    }
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
                    for(i=max_N_binding_rep-1;i>0;i--)
                    {
                        for(j=0;j<max_N_binding_act;j++)
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j]+product_of_freq*ratio_matrices[pos_of_mat_nH][i-1][j];
                            
                    }  
                    for(j=0;j<max_N_binding_act;j++)
                        ratio_matrices[pos_next_record][0][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][0][j];     
                }
                else
                {
                    /*find matrix(n)*/
                    pos_of_last_record=pos_next_record-1;  
                    for(i=0;i<max_N_binding_rep;i++)
                    {
                        for(j=0;j<max_N_binding_act;j++)
                            ratio_matrices[pos_next_record][i][j]=BS_info[m].Kd*ratio_matrices[pos_of_last_record][i][j];                                                   
                    }
                    ratio_matrices[pos_next_record][1][0]+=product_of_freq;
                } 
                break;
        }
    }

    sum=0.0;

    for(i=0;i<max_N_binding_rep;i++)
    {
        for(j=0;j<max_N_binding_act;j++)
//        for(j=0;j<((max_N_binding_act<max_unhindered_sites[0]+1-i)?max_N_binding_act:max_unhindered_sites[0]+1-i);j++)        
            sum+=ratio_matrices[pos_next_record][i][j];
    }   

    temp=0.0;
    for(i=0;i<max_N_binding_rep;i++)
    {  
        for(j=genotype->min_N_activator_to_transc[gene_id];j<max_N_binding_act;j++)
//        for(;j<((max_N_binding_act<max_unhindered_sites[0]+1-i)?max_N_binding_act:max_unhindered_sites[0]+1-i);j++)        
            temp+=ratio_matrices[pos_next_record][i][j];        
    } 
    state->P_A[gene_id]=(float)(temp/sum);
    
    temp=0.0;
    for(i=1;i<max_N_binding_rep;i++)
    {
        for(j=0;j<max_N_binding_act;j++)
            temp+=ratio_matrices[pos_next_record][i][j];
    }    
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
                                        float dt)
{
    int i,j;
    float ct, ect, one_minus_ect;
    float N_effector_molecules_bf_dt[genotype->protein_pool[genotype->nproteins-1][0][0]];
    float instantaneous_fitness = 0.0;
    float integrated_fitness = 0.0;
  
    /* store the numbers of the effector proteins encoded by each copy of gene before updating*/   
    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++)
        N_effector_molecules_bf_dt[i]=state->gene_specific_protein_number[genotype->protein_pool[genotype->nproteins-1][1][i]];
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
    /* now find out the protein numbers at end of dt interval and compute growth rate */   
    instantaneous_fitness = calc_fitness(&integrated_fitness, 
                                            genotype, 
                                            state, 
                                            N_effector_molecules_bf_dt, 
                                            dt);  
//    if(*end_state==0)
//    /*use 0 to indicate abnormal behavior of the program. I expect rounding error to raise this flag.
//    *In case this flag is raised, quit the current replicate of growth and rerun a replicate.*/
//    {
//        fperror=fopen(error_file,"a+");
//        LOG("at mut step %d",mut_step);
//        fclose(fperror);
//        return;
//    } 
    /* use the integrated growth rate to compute the cell size in the next timestep */
    state->cumulative_fitness += integrated_fitness;    
    /* update the instantaneous growth rate for the beginning of the next timestep */
    state->instantaneous_fitness = instantaneous_fitness;      
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
        state->protein_synthesis_index[gene_id] = (float)state->mRNA_aft_transl_delay_num[gene_id]*genotype->translation_rate[gene_id]/genotype->protein_decay_rate[gene_id];
        if(genotype->which_protein[gene_id]==genotype->nproteins-1)
        	return DO_NOTHING;
    	else // an mRNA of transcription factor is degraded, which can cause fluctuation in transcription factor concentrations.
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

/* while there are fixed events
     occuring in current t->dt window */
static int do_fixed_event(Genotype *genotype, 
                    CellState *state, 
                    GillespieRates *rates, 
                    Test *test,
                    Phenotype *timecourse,
                    float *dt,                           
                    int event)
{     
    int i, return_value;
    return_value=DO_NOTHING;
    switch (event) 
    {
        case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, state, rates, genotype); 
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, dt);
            state->cell_activated=1;
            break;
        case 3:     /* turn signal off*/ 
            *dt = state->signal_off_head->time - state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt); 
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));
            if(test->fixed_effector_effect)
                state->effect_of_effector=test->initial_effect_of_effector;
            else
                state->effect_of_effector='d';
            state->protein_number[N_SIGNAL_TF-1]=test->signal_off_strength;         
            break;
        case 4:     /*turn signal on*/
            *dt = state->signal_on_head->time - state->t;   
            update_protein_number_and_fitness(genotype, state, rates, *dt);  
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
            state->protein_number[N_SIGNAL_TF-1]=test->signal_on_strength;
            if(test->fixed_effector_effect)                               
                state->effect_of_effector=test->initial_effect_of_effector;            
            else                 
                state->effect_of_effector='b';            
            break;	
        case 5: /* finishing burn-in growth rate*/
            *dt=test->duration_of_burn_in_growth_rate-state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt);
            state->cumulative_fitness_after_burn_in=state->cumulative_fitness;           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6: /* mandatorily updating Pact and Prep*/
            *dt=state->t_to_update_probability_of_binding-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt);          
            break;
        case 7: /* update signal strength */
            *dt=state->change_signal_strength_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt);
            state->protein_number[N_SIGNAL_TF-1]=test->external_signal[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            break;     
        case 8:
            *dt=state->sampling_point_end_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt);
            delete_fixed_event_from_head(&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
            for(i=0;i<genotype->nproteins;i++)
                timecourse->protein_concentration[i*timecourse->total_time_points+timecourse->timepoint]=state->protein_number[i];
            for(i=0;i<genotype->ngenes;i++)
                timecourse->gene_specific_concentration[i*timecourse->total_time_points+timecourse->timepoint]=state->gene_specific_protein_number[i];                
            timecourse->instantaneous_fitness[timecourse->timepoint]=state->instantaneous_fitness;
            timecourse->timepoint++;   
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
 *  1 if a transcription event happens first before time t+dt
 *  2 if a translation initiation event  
 *  3 if the environmental signal is turned on 
 *  4 if the signal is turned off
 *  5 if burn-in growth ends
 *  6 if time to update Pact mandatorily
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
                                    Genotype *genotype)
{
    int gene_id;
    float concurrent;
    float endtime;
    /* recompute the delta-t based on difference between now and the time of transcription end */
    *dt = state->mRNA_transcr_time_end_head->time - state->t;   
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt);
    /* get the gene which is ending transcription */
    gene_id = state->mRNA_transcr_time_end_head->event_id;    
    /* increase number of mRNAs that are initializing translation*/
    (state->mRNA_under_transl_delay_num[gene_id])++;
    /* decrease the number of mRNAs undergoing transcription */
    (state->mRNA_under_transc_num[gene_id])--;
    /* delete the fixed even which has just occurred */
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
static int fixed_event_end_translation_init(   Genotype *genotype, 
                                        CellState *state, 
                                        GillespieRates *rates, 
                                        float *dt)
{
    int gene_id;    
    /* calc the remaining time till translation initiation ends */
    *dt = state->mRNA_transl_init_time_end_head->time - state->t;         
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt);
    /* get identity of gene that has just finished translating */
    gene_id=state->mRNA_transl_init_time_end_head->event_id; 
    /* there is one less mRNA that is initializing translation */
    (state->mRNA_under_transl_delay_num[gene_id])--;  
    /* delete the event that just happened */
    delete_fixed_event_from_head(&(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));    
    /* there is one more mRNA that produces protein */
    (state->mRNA_aft_transl_delay_num[gene_id])++;   
    /* update protein synthesis rate*/
    state->protein_synthesis_index[gene_id]= (float)state->mRNA_aft_transl_delay_num[gene_id]*genotype->translation_rate[gene_id]/genotype->protein_decay_rate[gene_id];
    
    if(genotype->which_protein[gene_id]==genotype->nproteins-1)//if the mRNA encodes a non-sensor TF, there could be a huge change in TF concentration
        return DO_NOTHING;
    else
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
    float N_proteins_at_now[genotype->ntfgenes];
    int gene_id; 
    float ct, ect, one_minus_ect;
    
    t_remaining=t_unreachable-state->t;
    dt=t_unreachable;  
 	protein_id=genotype->which_protein[which_gene];
    Kd=KD2APP_KD*genotype->Kd[protein_id];
    P_binding=state->protein_number[protein_id]/(state->protein_number[protein_id]+Kd);        

    /*determine whether the protein tends to increase or decrease concentration*/
    overall_rate=0.0;
    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
    {
        gene_id=genotype->protein_pool[protein_id][1][j];
        overall_rate+=(state->protein_synthesis_index[gene_id]-state->gene_specific_protein_number[gene_id])*genotype->protein_decay_rate[gene_id];
    }

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
                gene_id=genotype->protein_pool[protein_id][1][j];
                ct=genotype->protein_decay_rate[gene_id]*t_remaining;
                ect = exp(-ct);
                if (fabs(ct)<EPSILON) one_minus_ect=ct;
                else one_minus_ect = 1.0-ect;   
                N_at_end_of_simulation+=ect*state->gene_specific_protein_number[gene_id]+state->protein_synthesis_index[gene_id]*one_minus_ect;        
            }               
            if(N_at_end_of_simulation<N_proteins_cause_change)
               dt=t_unreachable; 
            else /*We need to solve an equation*/
            {                    
                for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                {
                    gene_id=genotype->protein_pool[protein_id][1][j];
                    N_proteins_at_now[j]=state->gene_specific_protein_number[gene_id];
                }
                dt=calc_tprime(genotype,state,N_proteins_at_now,t_remaining,N_proteins_cause_change,protein_id); 
            }
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
                gene_id=genotype->protein_pool[protein_id][1][j];
                ct=genotype->protein_decay_rate[gene_id]*t_remaining;
                ect = exp(-ct);
                if (fabs(ct)<EPSILON) one_minus_ect=ct;
                else one_minus_ect = 1.0-ect;   
                N_at_end_of_simulation+=ect*state->gene_specific_protein_number[gene_id]+state->protein_synthesis_index[gene_id]*one_minus_ect;        
            }               
            if(N_at_end_of_simulation>N_proteins_cause_change)
               dt=t_unreachable; 
            else
            {
                for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                {
                    gene_id=genotype->protein_pool[protein_id][1][j];
                    N_proteins_at_now[j]=state->gene_specific_protein_number[gene_id];
                }
                dt=calc_tprime(genotype,state,N_proteins_at_now,t_remaining,N_proteins_cause_change,protein_id); 
            }
        }
        *minimal_interval=(*minimal_interval<dt)?*minimal_interval:dt;  
    }
}



