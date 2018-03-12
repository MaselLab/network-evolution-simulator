/* -*- Mode: C; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe, Kun Xiong
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>
#include "numerical.h"
#include "lib.h"
#include "netsim.h"
#include "RngStream.h"

#define INITIALIZATION -1
#define DO_NOTHING -2

int MAXELEMENTS=100; 
//int min_act_to_transcribe[N_THREADS][MAX_BINDING];
const float TRANSLATION_TIME=2.0; 
const float TRANSCRIPTION_TIME=2.0;
const float PROB_ACTIVATING=0.62;
const float TRANSCRIPTINIT=6.75; 
const float MAX_REP_TO_INT_RATE=0.92;
const float BASAL_REP_TO_INT_RATE=0.15;
const float MAX_INT_TO_REP_RATE=4.11;
const float BASAL_INT_TO_REP_RATE=0.67;
const float MAX_INT_TO_ACT_RATE=3.3; 
const float BASAL_INT_TO_ACT_RATE=0.025;
const float MEAN_PROTEIN_DECAY_RATE=-2.83;
const float SD_PROTEIN_DECAY_RATE=0.286;
const float MEAN_ACT_TO_INT_RATE=1.27;
const float SD_ACT_TO_INT_RATE=0.226;
const float MEAN_MRNA_DECAY_RATE=-1.49;
const float SD_MRNA_DECAY_RATE=0.267;
const float MEAN_TRANSLATION_INIT_RATE=0.408;
const float SD_TRANSLATION_INIT_RATE=0.384;
const float MIN_Kd=1.0e-9;
const float MAX_Kd=1.0e-6;
const float log_MIN_Kd=-9.0;
const float log_MAX_Kd=-6.0;
const float NS_Kd=1.0e-5;
const float KD2APP_KD=1.8e10;
const float DEFAULT_UPDATE_INTERVAL=10.0; /*min*/
const float MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING=0.01;

/*Mutations*/
float SUBSTITUTION = 3.5e-10; /* susbstitution rate per site per cell division. Lynch 2008*/
/*the following mutations are enabled after burn-in*/
float DUPLICATION = 0.0;   
float SILENCING = 0.0;      
float MUTKINETIC = 4.0e-9; /* 1% subs in a gene (~1kb, including introns) will change kinetic rates and binding seq */        
float proportion_mut_binding_seq = 1.0;
float proportion_mut_identity = 1.0;
float proportion_mut_koff=1.0;
float proportion_mut_kdis=1.0;
float proportion_mut_mRNA_decay=1.0;
float proportion_mut_protein_decay=1.0;
float proportion_mut_translation_rate=1.0;
float proportion_mut_locus_specific_effect=1.0;
float proportion_mut_cooperation=0.0;
float proportion_effector2TF=0.1;
float mutational_regression_rate=0.5;
float sigma_ACT_TO_INT_RATE=0.39; 
float sigma_mRNA_decay=0.39; 
float sigma_protein_decay=0.39; 
float sigma_translation_init=0.39; 
float sigma_Kd=0.39;
float miu_ACT_TO_INT_RATE=1.75;
float miu_mRNA_decay=-0.49;
float miu_protein_decay=-2.83;
float miu_translation_init=-0.592;
float miu_Kd=-5.0;
const float MAX_ACT_TO_INT_RATE=64.7;
const float MIN_ACT_TO_INT_RATE=0.59;
const float MAX_MRNA_DECAY=0.54;
const float MIN_MRNA_DECAY=7.5e-4;
const float MAX_PROTEIN_DECAY=0.69;
const float MIN_PROTEIN_DECAY=3.9e-6;
const float MAX_TRANSLATION_RATE=33.0;
const float MIN_TRANSLATION_RATE=2.7e-3;
const float MAX_KD=1.0e-5;
const float MIN_KD=0.0;

/*fitness*/
float sampling_interval;
float saturate_cumulative_response_from_pulse;
float saturate_pulse_amplitude;
float opt_pulse_duration;
float sd_opt_pulse_duration;
float tolerable_delay_bf_pulse;
float Ne_saturate = 15000.0;
float c_transl=2.0e-6;
float bmax=1.0; 
float duration_of_burn_in_growth_rate; /* allow cells to reach (possiblly) steady growth*/
float growth_rate_scaling = 1.0; /* set default growth rate scaling factor */
int recalc_new_fitness; /*calculate the growth rate of the new genotype four more times to increase accuracy*/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   


/*initial conditions*/
int init_N_non_output_act=3;
int init_N_non_output_rep=3;
int init_N_output_act=1;
int init_N_output_rep=0;

/*Environments*/
float background_signal_strength = 0.0; /* number of TF0 protein */
float signal_off_strength=0.0; /* number of signal protein at basal level */
float signal_profile_matrix[N_THREADS][200][15];
float env1_minimal_peak_response;
float env2_minimal_peak_response;
float env1_response_amplification;
float env2_response_amplification;
float env1_benefit1;
float env1_benefit2;
float env1_max_duration_bias;            
float env2_benefit1;
float env2_benefit2;
float env2_max_duration_bias;
float env1_signal1_strength;
float env1_signal2_strength;
float env2_signal1_strength;
float env2_signal2_strength;
float env1_t_development;
float env2_t_development;
float env1_t_signal_on;    
float env1_t_signal_off;     
float env2_t_signal_on;
float env2_t_signal_off;
char env1_initial_effect_of_effector;
char env2_initial_effect_of_effector;
char env1_effect_of_effector;
char env2_effect_of_effector;
int env1_fixed_effector_effect;    
int env2_fixed_effector_effect;
float env1_occurence;             /* one environment can occur more frequently than*/                                
float env2_occurence;             /* the other*/ 

char set_base_pair(float x) 
{
    char base;
    if (x<0.25)
        base = 'a';
    else if (x<0.5)
        base = 'c';
    else if (x<0.75)
        base = 'g';
    else 
        base = 't';  
    return base;
}

void initialize_sequence(char *Seq, 
                         int len,                         
                         int num_elements,
                         RngStream RS)
{
    float x;
    int i;
    int current_element = len/num_elements;
    int pos_n;

    for (i=0; i<len; i++) 
    {
        pos_n = (i / current_element)*current_element + i % current_element;  
        x = RngStream_RandU01(RS);     
        Seq[pos_n] = set_base_pair(x);
    }
}

/*This function initialize kinetic constants for gene expression, as well the identities of TFs*/
void initialize_genotype_fixed(Genotype *genotype, RngStream RS)
{
    int i, j, counter1,counter2,counter3;
    /* the first N_SIGNAL_TF genes encode the sensor TFs. The concentration of a sensor TF
     * is determined by certain environmental signal*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {  
        #if RANDOM_COOPERATION_LOGIC        
            genotype->min_act_to_transc[i]=RngStream_RandInt(RS,1,2); //if one activator is sufficient to induce expression, the gene is regualted by OR gate.
        #else
            genotype->min_N_activator_to_transc[i]=1; 
            genotype->min_N_activator_to_transc[genotype->ngenes-1]=1;
        #endif             
        /* tf affinity */
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
        while(genotype->Kd[i]>=MAX_Kd || genotype->Kd[i]<=MIN_Kd)
            genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
        /* mRNA decay */
        genotype->mRNA_decay_rate[i] = pow(10.0,SD_MRNA_DECAY_RATE*gasdev(RS)+MEAN_MRNA_DECAY_RATE);
        while(genotype->mRNA_decay_rate[i]>=MAX_MRNA_DECAY || genotype->mRNA_decay_rate[i]<=MIN_MRNA_DECAY)
            genotype->mRNA_decay_rate[i] = pow(10.0,SD_MRNA_DECAY_RATE*gasdev(RS)+MEAN_MRNA_DECAY_RATE);
        /* protein decay */
        genotype->protein_decay_rate[i] = pow(10.0,SD_PROTEIN_DECAY_RATE*gasdev(RS)+MEAN_PROTEIN_DECAY_RATE);  
        while(genotype->protein_decay_rate[i]>=MAX_PROTEIN_DECAY || genotype->protein_decay_rate[i]<=MIN_PROTEIN_DECAY)
            genotype->protein_decay_rate[i] = pow(10.0,SD_PROTEIN_DECAY_RATE*gasdev(RS)+MEAN_PROTEIN_DECAY_RATE);
        /* translation rate */
        genotype->translation_rate[i] = pow(10.0,SD_TRANSLATION_INIT_RATE*gasdev(RS)+MEAN_TRANSLATION_INIT_RATE);  
        while(genotype->translation_rate[i]>=MAX_TRANSLATION_RATE || genotype->translation_rate[i]<=MIN_TRANSLATION_RATE)
            genotype->translation_rate[i] = pow(10.0,SD_TRANSLATION_INIT_RATE*gasdev(RS)+MEAN_TRANSLATION_INIT_RATE); 
        /*PIC disassembly rate*/
        genotype->active_to_intermediate_rate[i]=pow(10.0,SD_ACT_TO_INT_RATE*gasdev(RS)+MEAN_ACT_TO_INT_RATE);  
        while(genotype->active_to_intermediate_rate[i]>=MAX_ACT_TO_INT_RATE || genotype->active_to_intermediate_rate[i]<=MIN_ACT_TO_INT_RATE)
            genotype->active_to_intermediate_rate[i]=pow(10.0,SD_ACT_TO_INT_RATE*gasdev(RS)+MEAN_ACT_TO_INT_RATE);  
    }          
    /* assign tf identity*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        /*initialize the signal as an activator of every gene*/
        genotype->locus_specific_TF_behavior[i][0]=ACTIVATOR;
        /*for each gene, the first init_N_non_output_act TF are activators, and init_N_non_output_rep are repressors*/
        for(j=1;j<init_N_non_output_act+1;j++)
            genotype->locus_specific_TF_behavior[i][j]=ACTIVATOR;
        for(;j<init_N_non_output_act+1+init_N_non_output_rep;j++)
            genotype->locus_specific_TF_behavior[i][j]=REPRESSOR;
        /*do the same for activators and repressors that are also the output*/
        for(;j<init_N_non_output_act+1+init_N_non_output_rep+init_N_output_act;j++)
            genotype->locus_specific_TF_behavior[i][j]=ACTIVATOR;
        for(;j<init_N_non_output_act+1+init_N_non_output_rep+init_N_output_act+init_N_output_rep;j++)
            genotype->locus_specific_TF_behavior[i][j]=REPRESSOR;        
    }
    for(i=0;i<init_N_non_output_act+N_SIGNAL_TF;i++)
    {
        genotype->protein_identity[i][0]=ACTIVATOR;
        genotype->protein_identity[i][1]=NON_OUTPUT_PROTEIN;
    }
    for(;i<init_N_non_output_act+N_SIGNAL_TF+init_N_non_output_rep;i++)
    {
        genotype->protein_identity[i][0]=REPRESSOR;
        genotype->protein_identity[i][1]=NON_OUTPUT_PROTEIN;
    }
    j=0;
    for(;i<init_N_non_output_rep+init_N_non_output_act+init_N_output_act+N_SIGNAL_TF;i++)
    {
        genotype->protein_identity[i][0]=ACTIVATOR;
        genotype->protein_identity[i][1]=i;
        genotype->output_protein_id[j]=i;
        j++;
    }
    for(;i<init_N_non_output_rep+init_N_non_output_act+init_N_output_act+init_N_output_rep+N_SIGNAL_TF;i++)
    {
        genotype->protein_identity[i][0]=REPRESSOR;
        genotype->protein_identity[i][1]=i;
        genotype->output_protein_id[j]=i;
        j++;
    }
    genotype->N_act=init_N_non_output_act+init_N_output_act;
    genotype->N_rep=init_N_non_output_rep+init_N_output_rep;
    /* parameterize sensor TF*/
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        genotype->mRNA_decay_rate[i]=0.0; // we assume environmental signal toggles the state of sensor TF between active and inactive 
        genotype->protein_decay_rate[i]=0.0; // the concentration of sensor TF is constant.
        genotype->translation_rate[i]=0.0;
        genotype->active_to_intermediate_rate[i]=0.0; 
        genotype->protein_identity[i][0]=ACTIVATOR; /*make sensor TF an activator*/
        genotype->protein_identity[i][1]=NON_OUTPUT_PROTEIN;
        genotype->N_act++;
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
    }
#if RANDOMIZE_SIGNAL2
    #if N_SIGNAL_TF==2
        if(RngStream_RandU01(RS)<=0.5) // we assume there is a background "on" signal, which is sensor TF 0, in the network.
            genotype->protein_identity[1]=1; // Other sensor TFs can be either activators or repressors.
        else
        {
            genotype->protein_identity[1]=0;
            genotype->N_act--;
            genotype->N_rep++;
        }
    #endif
#endif   
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype, RngStream RS)
{ 
    int i,k;

    genotype->ngenes=N_SIGNAL_TF+init_N_output_act+init_N_output_rep+init_N_non_output_act+init_N_non_output_rep; /*including the signal genes and 1 selection gene*/  
    genotype->nproteins=genotype->ngenes;  /*at initialization, each protein is encoded by one copy of gene*/   
    genotype->n_output_genes=init_N_output_act+init_N_output_rep;
    genotype->n_output_proteins=genotype->n_output_genes;
    /*at initialization, each copy of gene should have a unique cis-regulatory sequence*/
    for(i=0;i<genotype->ngenes;i++)
    {    
        genotype->which_cluster[i]=i; 
        genotype->cisreg_cluster[i][0]=i;
    }     
    /* initially, each protein has only one copy of gene*/    
    for(i=0;i<genotype->nproteins;i++)
    {
        genotype->protein_pool[i][0][0]=1;
        genotype->protein_pool[i][1][0]=i;
        genotype->which_protein[i]=i;
    }  
    /* put output genes to the end*/
    for(i=0;i<genotype->n_output_proteins;i++)
        genotype->output_protein_id[i]=genotype->nproteins-genotype->n_output_proteins+i;
    
    initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*NGENES, genotype->ngenes, RS);  // initialize cis-reg sequence
    initialize_sequence((char *)genotype->tf_binding_seq, TF_ELEMENT_LEN*NGENES, genotype->ngenes, RS);    //initialize binding sequence of TFs    
    /* We now generate the complementary sequence of BS that are on the non-template strand.
     * The complementary sequence is used to search for BS that on the non-template strand.  
     * We also assume that all the TFs can work on both strands, but cen induce expression in one direction.*/  
    for(i=0;i< genotype->ngenes;i++)
    {        
        for(k=0;k<TF_ELEMENT_LEN;k++)
        {
            switch (genotype->tf_binding_seq[i][TF_ELEMENT_LEN-k-1])
            {
                case 'a': genotype->tf_binding_seq_rc[i][k]='t'; break;
                case 't': genotype->tf_binding_seq_rc[i][k]='a'; break;
                case 'c': genotype->tf_binding_seq_rc[i][k]='g'; break;
                case 'g': genotype->tf_binding_seq_rc[i][k]='c'; break;
            }
        }        
    }
    initialize_genotype_fixed(genotype, RS);     
    #if !SET_BS_MANUALLY    
        calc_all_binding_sites(genotype);
    #endif
}

/*
 * Manually set binding sites. 
 * Check calc_all_binding_sites to understand this function.
 *
 */
void set_binding_sites(Genotype *genotype,int *N_BS_per_gene, int (*tf_id_per_site)[5], int (*hindrance_per_site)[5], float (*Kd_per_site)[5])
{
    int i,j;
    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        genotype->binding_sites_num[i]=N_BS_per_gene[i];
        genotype->max_hindered_sites[i]=0;
        for(j=0;j<N_BS_per_gene[i];j++)
        {
            genotype->all_binding_sites[i][j].tf_id=tf_id_per_site[i][j];
            genotype->all_binding_sites[i][j].N_hindered=hindrance_per_site[i][j];
            genotype->all_binding_sites[i][j].Kd=Kd_per_site[i][j];
            genotype->max_hindered_sites[i]=(genotype->max_hindered_sites[i]>hindrance_per_site[i][j])?
                                            genotype->max_hindered_sites[i]:hindrance_per_site[i][j];
            genotype->all_binding_sites[i][j].mis_match=0;
        }       
    }
}

/*
 * compute the list binding sites for specified gene and gene copy
 */
void calc_all_binding_sites_copy(Genotype *genotype, int gene_id)
{
    int i, j, k;
    int match,match_rc; // number of nucleotide that matches the binding sequence of TF, in a binding site in the coding and in the non-coding strand.    
    int N_hindered_BS=0;   
    int N_binding_sites=0;
    int start_TF;
    FILE *fperror;
    genotype->N_act_BS[gene_id]=0;
    genotype->N_rep_BS[gene_id]=0;
    genotype->max_hindered_sites[gene_id]=0;  
    //some helper pointer 
    char *tf_binding_seq;
    char *cis_seq;
    char *tf_binding_seq_rc; 
    cis_seq=&(genotype->cisreg_seq[gene_id][0]); 
  
    for(i=3; i < CISREG_LEN-TF_ELEMENT_LEN-3; i++) /* scan promoter */
    {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(genotype->all_binding_sites[gene_id][j].BS_pos> i-TF_ELEMENT_LEN-2*HIND_LENGTH)
                    N_hindered_BS++;
            }
        }  
        /* loop through TF proteins */        
        #if !DIRECT_REG 
            if(genotype->which_protein[gene_id]==genotype->nproteins-1) // if the gene is an effector gene
                start_TF=N_SIGNAL_TF;// the environmental signals cannot directly regulate the selection gene
            else
                start_TF=0;
        #else
            start_TF=0;
        #endif
        
        for (k=start_TF; k < genotype->nproteins; k++) 
        {    
            tf_binding_seq=&(genotype->tf_binding_seq[k][0]);
            tf_binding_seq_rc=&(genotype->tf_binding_seq_rc[k][0]);            
            /*find BS on the template strand*/
            match=0;
            for (j=i; j < i+TF_ELEMENT_LEN; j++) /*calculate the number of nucleotides that match in each [i,i+TF_ELEMENT_LEN] window. The window slides by 1 each time when scanning the promoter*/
                if (cis_seq[j] == tf_binding_seq[j-i]) match++; 
            if (match >= NMIN)
            {  
                if (N_binding_sites + 1 >= genotype->N_allocated_elements) 
                {  
                    while(genotype->N_allocated_elements<=N_binding_sites+1)
                        genotype->N_allocated_elements+=100;
                   
                    for(j=0;j<NGENES;j++)
                    {
                        genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                        if (!genotype->all_binding_sites[j]) 
                        {
                            fperror=fopen(error_file,"a+");
                            LOG("error in calc_all_binding_sites_copy\n");
                            fclose(fperror);
                            exit(-1);                                       
                        }     
                    }                    
                }                
                genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                      
                genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match)/(TF_ELEMENT_LEN-NMIN+1));
                genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ; 
                genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match;             
                genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                N_hindered_BS++;              
                N_binding_sites++;
                if(genotype->locus_specific_TF_behavior[gene_id][k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand.*/
            {
                match_rc=0;
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_binding_seq_rc[j-i]) match_rc++;
                if (match_rc >= NMIN)
                {
                    /**********************************************************************/     
                    if (N_binding_sites + 1 >= genotype->N_allocated_elements) 
                    {  
                        while(genotype->N_allocated_elements<=N_binding_sites+1)
                            genotype->N_allocated_elements+=100;

                        for(j=0;j<NGENES;j++)
                        {
                            genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                            if (!genotype->all_binding_sites[j]) 
                            {
                                fperror=fopen(error_file,"a+");
                                LOG("error in calc_all_binding_sites_copy\n");
                                fclose(fperror);
                                exit(-1);                                       
                            }     
                        }                    
                    } 
                    /************************************************************************************************************/
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                                     
                    genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match_rc)/(TF_ELEMENT_LEN-NMIN+1));
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
                    genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match_rc;
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                    N_hindered_BS++;                  
                    N_binding_sites++;  //two binding sites on different strands can also hinder each other                  
                    if(genotype->locus_specific_TF_behavior[gene_id][k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                }
            } 
        }/* looping through TFs ends */
    }/*end of promoter scanning*/ 
    
    genotype->binding_sites_num[gene_id]=N_binding_sites;  
    genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);
    /* calculate max_hindered_sites */    
    for(i=0;i<N_binding_sites;i++)
    {
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id] > genotype->all_binding_sites[gene_id][i].N_hindered)?
                                      genotype->max_hindered_sites[gene_id] : genotype->all_binding_sites[gene_id][i].N_hindered;           
    }
    
#if FORCE_OR_GATE
    float temp;    
    if(genotype->N_motifs[0]!=0 && genotype->gene_in_core_C1ffl[gene_id]==1) //
    {
        temp=1.0;
        for(i=0;i<genotype->binding_sites_num[gene_id];i++)
        {
            if(genotype->all_binding_sites[gene_id][i].tf_id==N_SIGNAL_TF-1) 
                temp=(temp>genotype->all_binding_sites[gene_id][i].Kd)?genotype->all_binding_sites[gene_id][i].Kd:temp;  // find a strong binding site                      
        }
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].tf_id = N_SIGNAL_TF-1;
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd = temp;                   
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].BS_pos = 2*CISREG_LEN;
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].mis_match = 0;
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].N_hindered = 0;
        genotype->binding_sites_num[gene_id]++;
        genotype->N_act_BS[gene_id]++;

        for(i=1;i<genotype->nproteins;i++)
        {
            if(genotype->TF_in_core_C1ffl[gene_id][i]==1)
            {
                temp=1.0;
                for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                {
                    if(genotype->all_binding_sites[gene_id][j].tf_id==i) 
                        temp=(temp>genotype->all_binding_sites[gene_id][j].Kd)?genotype->all_binding_sites[gene_id][j].Kd:temp;  // find a strong binding site                      
                }
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].tf_id = i;
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd = temp;                   
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].BS_pos = (2+i)*CISREG_LEN;
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].mis_match = 0;
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].N_hindered = 0;
                genotype->binding_sites_num[gene_id]++;
                genotype->N_act_BS[gene_id]++;                
            }
        }
    }    
#endif
    
    /* calculate max_unhindered_sites */
    /* max_unhindered_sites is maximum number of TFs that can bind to a cis-reg sequence at the same time*/
    /* We use it to faciliate the calculation of Pact and Prep. See calc_TF_dist_from_all_BS for its usage.*/
    int act_BS[MAXELEMENTS][2],rep_BS[MAXELEMENTS][2];
    int N_act_BS,N_rep_BS;    
    N_act_BS=1;
    N_rep_BS=1;
    for(i=0;i<N_binding_sites;i++) /* make lists BS by their types*/    
    {
        if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id][0]==1)
        {
            act_BS[N_act_BS][0]=i;
            N_act_BS++;
        } 
        else
        {
            rep_BS[N_rep_BS][0]=i;
            N_rep_BS++;
        }
    }
    /*calculate the maximum number of activator binding sites that do not hinder each other*/
    /*Assuming that when site n is bound, at most x sites among site 1 to site n can be bound at the same time 
    *If site n+1 does not hinder n, then when site n+1 is bound, at most x+1 sites among site 1 to n+1 can be bound
    *If site n+1 hinders n, we check if it hinders site n-1,n-2.. until one site n-m which is not hindered. 
    *Obviously, site n-m+1 must be hindered by site n, therefore at most x-1 or x site among site 1 to n-m can be bound,
    *depending on whether n-m is hindered by n. This means at most x or x+1 sites among site 1 to n+1 can be bound. 
    *This means as n increases the maximum number of sites that can be bound at the same time won't decrease. 
    *We will see the maximum number of binding sites that won't hinder each other when n=N_act_BS.*/
    act_BS[0][0]=-1;
    act_BS[0][1]=0; 
    for(i=1;i<N_act_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][act_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][act_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        act_BS[i][1]=act_BS[j][1]+1;
    }
    /*calculate the maximum number of repressor binding sites that do not hinder each other*/
    rep_BS[0][0]=-1;
    rep_BS[0][1]=0;
    for(i=1;i<N_rep_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][rep_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][rep_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        rep_BS[i][1]=rep_BS[j][1]+1;
    }
    genotype->max_unhindered_sites[gene_id][1]=act_BS[N_act_BS-1][1];
    genotype->max_unhindered_sites[gene_id][2]=rep_BS[N_rep_BS-1][1];
}

/*
 * compute the list of binding sites for the specified number of gene
 * copies
 */
void calc_all_binding_sites(Genotype *genotype)
{    
    int gene_id;
    if(genotype->N_allocated_elements<MAXELEMENTS)
    {
        for(gene_id=0;gene_id<NGENES;gene_id++)
            genotype->all_binding_sites[gene_id]=realloc(genotype->all_binding_sites[gene_id], MAXELEMENTS*sizeof(AllTFBindingSites));
        genotype->N_allocated_elements=MAXELEMENTS;
    }
    
    for(gene_id=N_SIGNAL_TF;gene_id < genotype->ngenes;gene_id++)
    {        
        if(genotype->recalc_TFBS[gene_id]) /* do not calculate the binding sites if there's no mutation in the promoter or in TF binding seq*/
        {            
            calc_all_binding_sites_copy(genotype,gene_id);          
            genotype->recalc_TFBS[gene_id]=0;
        }
    }
}

int mod(int a, int b) // b is the base
{
    if(a>=0)
        return (a%b);
    else
        return (abs(a+b)%b);
}

/*Calculate probability of binding configurations*/
void calc_TF_dist_from_all_BS( AllTFBindingSites *BS_info,
                                int nproteins,
                                int max_N_hindered_BS,
                                int N_BS,                                                                     
                                int protein_identity[NPROTEINS],
                                int max_unhindered_sites[3],
                                float protein_number[NGENES],
                                int min_act_to_transcr,
                                float *Pact,
                                float *Prep,
								float *PaNr,
								float *Pno_TF)
{
    int max_N_binding_act=max_unhindered_sites[1]+1; //Binding configurations can contain at most x activators, plus 1 type of configurations that don't have activators at all. 
    int max_N_binding_rep=max_unhindered_sites[2]+1; //Binding configurations can contain at most y repressors, plus 1 type of configurations that don't have repressors at all. 
    double ratio_matrices[N_BS+1][max_N_binding_rep][max_N_binding_act]; 
    double sum;    
    double product_of_freq; 
    int pos_of_last_record;  
    int pos_of_mat_nH;
    int pos_next_record;
    int i,j,k,m,n;
    double temp;
    
    /* initializing matrices to all zeros */
    for(i=0;i<max_N_hindered_BS+1;i++)
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
    if(protein_identity[BS_info[0].tf_id]==1) // if a activator binds to BS 0   
        ratio_matrices[pos_next_record][0][1]=protein_number[BS_info[0].tf_id];
    else    
        ratio_matrices[pos_next_record][1][0]=protein_number[BS_info[0].tf_id]; 
    /*keep calculating distribution from the remaining BS*/
    for(m=1;m<N_BS;m++)
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
        switch(protein_identity[BS_info[m].tf_id])
        {
            case 1: // a BS of activators              
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
                         *matrix. Updating the x column of matrix(n+1) uses
                         *the x-1 column of matrix(n-H). To avoiding changing the values
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

            case 0: // a BS of repressors            
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
        for(j=min_act_to_transcr;j<max_N_binding_act;j++)
//        for(;j<((max_N_binding_act<max_unhindered_sites[0]+1-i)?max_N_binding_act:max_unhindered_sites[0]+1-i);j++)        
            temp+=ratio_matrices[pos_next_record][i][j];        
    } 
    *Pact=(float)(temp/sum);
    
    temp=0.0;
    for(i=1;i<max_N_binding_rep;i++)
    {
        for(j=0;j<max_N_binding_act;j++)
            temp+=ratio_matrices[pos_next_record][i][j];
    }    
    *Prep=(float)(temp/sum);
    
    temp=0.0;
	for(j= min_act_to_transcr;j<max_N_binding_act;j++)
		temp+=ratio_matrices[pos_next_record][0][j];       
    for(i=1;i<max_N_binding_rep;i++)
    {  
        for(j=2*i;j<max_N_binding_act;j++)
//        for(;j<((max_N_binding_act<max_unhindered_sites[0]+1-i)?max_N_binding_act:max_unhindered_sites[0]+1-i);j++)        
            temp+=ratio_matrices[pos_next_record][i][j];        
    }     
	*PaNr = (float)(temp / sum);
   
	*Pno_TF = (float)(ratio_matrices[pos_next_record][0][0] / sum);

    // end of the forward algorithm   
}

/*Add fixed event to queue*/
int add_fixed_event(int i,
                    float t,
                    FixedEvent **start,
                    FixedEvent **last)
{
    FixedEvent *newtime;
    int pos;
    FILE *fp;

    newtime = malloc(sizeof*newtime);
    if (!newtime) 
    {
        fp=fopen("output.txt","a+");
        fprintf(fp,"error in add_fixed_event\n");
        fclose(fp); 
        exit(1);
    }
    newtime->event_id = i;  
    newtime->time = t;
    pos = sls_store(newtime, start, last);
    return pos;
}

/**returns 0 if new fixed event won't happen concurrently with any exisiting event*/
int check_concurrence(  float t, //new fixed event is scheduled to happen at t
                        FixedEvent *translation, 
                        FixedEvent *transcription, 
                        FixedEvent *signal_on, 
                        FixedEvent *signal_off,
                        FixedEvent *burn_in,
                        float t_to_update_probability_of_binding,
                        FixedEvent *change_signal)
{   
    while(translation!=NULL)
    {
        if(t==translation->time)            
            return 1;
        translation=translation->next;
    }
    while(transcription!=NULL)
    {
        if(t==transcription->time)        
            return 1;
        transcription=transcription->next;
    }
    while(signal_on!=NULL)
    {
        if(t==signal_on->time)        
            return 1; 
        signal_on=signal_on->next;
    }
    while(signal_off!=NULL)
    {
        if(t==signal_off->time)            
            return 1;
        signal_off=signal_off->next;
    }
    while(burn_in!=NULL)
    {
        if(t==burn_in->time)        
            return 1; 
        burn_in=burn_in->next;
    }
    while(change_signal!=NULL)
    {
        if(t==change_signal->time)        
            return 1; 
        change_signal=change_signal->next;        
    }
    if(t==t_to_update_probability_of_binding)
        return 1;        
    return 0;
}

/*delete linked table from anywhere*/
/*This function is used only to pick up mRNA that is under translation
 *initiation*/
void delete_fixed_event(int gene_x,                      
                        int mRNA_y_of_gene_x,
                        FixedEvent **head,
                        FixedEvent **tail)
{
    FixedEvent *info, *lastinfo=NULL;
    int j, done;
    FILE *fp;
  
    j = -1;
    done = 0;
    info = *head;
    while (info) 
    {
        if (info->event_id==gene_x) 
        {
            j++; //we are looking for the yth mRNA of gene x.
            if (j == mRNA_y_of_gene_x) 
            {
                if (info == *head) //if we found mRNA y in the head of the queue
                {
                    *head = info->next; //the original 2nd in queue becomes the new head
                    if (info == *tail) //if there is only one event in queue, deleting the event leaves an empty queue. 
                        *tail = NULL; // Therefore the tail of the queue is NULL (the head is set to null in the upper line),                       
                } 
                else 
                {
                    lastinfo->next = info->next;
                    if (info == *tail) //if we are going to delete the tail of the queue
                        *tail = lastinfo; //the original second to last event becomes the new tail
                }                
                done = 1; //found mRNA y!
                break;
            } 
            else 
            {
                lastinfo = info;
                info = info->next;
            }
        } 
        else 
        {
            lastinfo = info;
            info = info->next;
        }
    }
    if (done == 0) //if could not find mRNA y
    {
        fp=fopen("output.txt","a+");
        fprintf(fp,"error in delete_fixed_event\n");
        fclose(fp);        
        exit(1);
    }
    free(info);
}

void delete_fixed_event_from_head(FixedEvent **head,FixedEvent **tail)
{
    FixedEvent *info;

    info = *head;
    *head = info->next;
    if (*tail == info) 
        *tail = NULL;
    free(info);
}

/*
 * initialize the cell state with the specified initial protein
 * concentration, mean mRNA number and mRNA decay
 */
void initialize_cell(   Genotype *genotype,
                        CellState *state,                        
                        int ngenes,
                        int nproteins,
                        int protein_pool[NPROTEINS][2][NGENES], 
                        float mRNAdecay[NGENES],
                        int init_mRNA_number[NGENES],
                        float init_protein_number[NPROTEINS], 
                        RngStream RS,
                        float tdevelopment,
                        float t_signal_on)
{
    int i, j;
    /* start fitness at 1.0 */
	state->sensitivity[0] = 0.0;	
	state->cumulative_cost = 0.0;
    state->Duration_pulse=0.0;
    state->T_pulse_on=tdevelopment;
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
    state->cell_activated=0;
    state->t_to_update_probability_of_binding=TIME_INFINITY;
    for(i=0;i<MAX_OUTPUT_PROTEINS;i++)
    {
        state->cumulative_fitness[i] = 0.0;     
        state->cumulative_fitness_after_burn_in[i] = 0.0;   
        state->instantaneous_fitness[i] = 0.0; 
        state->threshold_response_to_bias[i]=0.0;
        state->basal_response[i]=0.0;
    }
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
    /*initialize gene state, mRNA number*/
    for (i=N_SIGNAL_TF; i < ngenes; i++) 
    {
        state->gene_state[i]= REPRESSED;       
        state->mRNA_aft_transl_delay_num[i]=init_mRNA_number[i];
        state->mRNA_under_transl_delay_num[i]=0;
        state->mRNA_under_transc_num[i]=0;
        state->last_P_A[i]=0.0;
        state->last_P_R[i]=0.0;
        state->last_P_A_no_R[i]=0.0;
        state->last_P_NotA_no_R[i]=0.0;
        state->protein_synthesis_index[i]=(float)state->mRNA_aft_transl_delay_num[i]*genotype->translation_rate[i]/genotype->protein_decay_rate[i];
    }       
    /* initiate protein concentration*/
    for (i=N_SIGNAL_TF; i < ngenes; i++) 
        state->gene_specific_protein_number[i] = init_protein_number[i];    
    for(i=N_SIGNAL_TF;i<nproteins;i++)
    {
        state->protein_number[i]=0.0;        
        for(j=0;j<protein_pool[i][0][0];j++)
            state->protein_number[i]+=state->gene_specific_protein_number[protein_pool[i][1][j]];
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
    if(duration_of_burn_in_growth_rate!=0.0)
        add_fixed_event(-1,duration_of_burn_in_growth_rate,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail)); 
    else
        add_fixed_event(-1,(float)TIME_INFINITY,&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));                
    /*plot protein concentration and fitness vs time*/ 
#if !JUST_PLOTTING
#if PEAK_SEARCH
    float t;
    int N_data_points; 
    t=t_signal_on+TIME_OFFSET;
    N_data_points=(int)((tdevelopment-t_signal_on)/sampling_interval);
    for(i=0;i<N_data_points;i++)
    {
        add_fixed_event(-1,t,&(state->sampling_point_end_head),&(state->sampling_point_end_tail)); //get a timepoint each minute
        t+=sampling_interval;            
    } 
    state->sampled_response=(float*)calloc(N_data_points,sizeof(float)); 
#elif CHEMOTAXIS
    add_fixed_event(-1,t_signal_on-1.0,&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
#endif    
#else
    float t;    
    int N_data_points;
    N_data_points=(int)((tdevelopment-t_signal_on)/sampling_interval);
    state->sampled_response=(float*)calloc(N_data_points,sizeof(float)); 
    t=TIME_OFFSET;   
    for(i=0;i<N_TIMEPOINTS;i++)
    {
        add_fixed_event(-1,t,&(state->sampling_point_end_head),&(state->sampling_point_end_tail)); //get a timepoint each minute
        t+=1.0;            
    }      
#endif
}

/*
 * Set how the environmental signal should change
 */
void set_signal(CellState *state, 
                float duration_signal_on, 
                float duration_signal_off, 
                float *signal_profile,                 
                float tdevelopment,
                float signal_on_strength)
{
    float t=0.0;     
    char flag;    
    if(signal_profile==NULL)   
    {
        flag='f'; 
        state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;    //always start with signal on
        #if N_SIGNAL_TF==2
            state->protein_number[0]=background_signal_strength;
        #endif      
        while(t<tdevelopment)
        {
            if(flag=='o')
            {
                add_fixed_event(-1,t+duration_signal_on,&(state->signal_off_head),&(state->signal_off_tail));
                flag='f';
                t=t+duration_signal_on;                    
            }    
            else
            {
                add_fixed_event(-1,t+duration_signal_off,&(state->signal_on_head),&(state->signal_on_tail));
                flag='o';
                t=t+duration_signal_off;
            }
        } 
    }
    else
    {
        int time_point=1;
        state->protein_number[N_SIGNAL_TF-1]=signal_profile[0];
        t=10.0;
        while(t<tdevelopment)
        {
            add_fixed_event(time_point,t,&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            time_point++;
            t+=10.0;
        } 
    }
}

/*
 * Calculate the rates of all Gillespie events
 */
void calc_all_rates(Genotype *genotype,
                    CellState *state,
                    GillespieRates *rates, 
                    float tdevelopment,
                    int UPDATE_WHAT,
                    int thread_id)
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
                        calc_TF_dist_from_all_BS(   genotype->all_binding_sites[i],
                                                    genotype->nproteins,
                                                    genotype->max_hindered_sites[i],
                                                    genotype->binding_sites_num[i],                                                                
                                                    genotype->locus_specific_TF_behavior[i], 
                                                    genotype->max_unhindered_sites[i],
                                                    state->protein_number,
                                                    genotype->min_N_activator_to_transc[i],
                                                    &(state->P_A[i]),
                                                    &(state->P_R[i]),
													&(state->P_A_no_R[i]),
                                                    &(state->P_NotA_no_R[i]));
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
                    calc_TF_dist_from_all_BS(   genotype->all_binding_sites[i],
                                                genotype->nproteins,
                                                genotype->max_hindered_sites[i],
                                                genotype->binding_sites_num[i],                                                                
                                                genotype->locus_specific_TF_behavior[i], 
                                                genotype->max_unhindered_sites[i],
                                                state->protein_number,
                                                genotype->min_N_activator_to_transc[i],
                                                &(state->P_A[i]),
                                                &(state->P_R[i]),
												&(state->P_A_no_R[i]),
                                                &(state->P_NotA_no_R[i]));
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
        switch (state->gene_state[i])
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
        while(genotype->cisreg_cluster[cluster_id][0]!=-1) //check if Pact changes too much
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
            calc_leaping_interval(genotype,state,&interval_to_update_probability_of_binding,tdevelopment,UPDATE_WHAT);  
    
        /*Update the next time that Pact will be updated mandatorily*/
        t_to_update_probability_of_binding=state->t+interval_to_update_probability_of_binding;
        concurrent=check_concurrence(   t_to_update_probability_of_binding,
                                        state->mRNA_transl_init_time_end_head,
                                        state->mRNA_transcr_time_end_head,
                                        state->signal_on_head,
                                        state->signal_off_head,
                                        state->burn_in_growth_rate_head,
                                        state->t_to_update_probability_of_binding,
                                        state->change_signal_strength_head);
        while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
        {
            t_to_update_probability_of_binding+=TIME_OFFSET;
            concurrent=check_concurrence(   t_to_update_probability_of_binding,
                                            state->mRNA_transl_init_time_end_head,
                                            state->mRNA_transcr_time_end_head,
                                            state->signal_on_head,
                                            state->signal_off_head,
                                            state->burn_in_growth_rate_head,
                                            state->t_to_update_probability_of_binding,
                                            state->change_signal_strength_head);        
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
 
int does_fixed_event_end(CellState *state, float t) 
{
    int retval;
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
    t7= state->change_signal_strength_head ? state->change_signal_strength_head->time : TIME_INFINITY;
    t8=state->sampling_point_end_head?state->sampling_point_end_head->time:TIME_INFINITY;
    if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5) &&(t1<=t6) && (t1<=t7)&&(t1<=t8))
    {
        retval = 1;	
    }
    else if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5)&&(t2<=t6)&&(t2<=t7)&&(t2<=t8))
    { 
        retval = 2;
    }  
    else if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5)&&(t3<=t6) &&(t3<=t7)&&(t3<=t8)) 
    {
        retval = 3;
    }
    else if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5)&&(t4<=t6)&&(t4<=t7)&&(t4<=t8)) 
    {
        retval = 4;
    }               
    else if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4)&&(t5<=t6)&&(t5<=t7)&&(t5<=t8))
    {
        retval = 5;
    }             
    else if((t6 <= t1) && (t6 <= t) && (t6 <= t2) && (t6 <= t3) && (t6<=t4)&&(t6<=t5)&&(t6<=t7)&&(t6<=t8))
    {
        retval=6;
    }
    else if((t7 <= t1) && (t7 <= t) && (t7 <= t2) && (t7 <= t3) && (t7<=t4)&&(t7<=t5)&&(t7<=t6)&&(t7<=t8))
    {
        retval=7;
    }
    else if((t8 <= t1) && (t8 <= t) && (t8 <= t2) && (t8 <= t3) && (t8<=t4)&&(t8<=t5)&&(t8<=t6)&&(t8<=t7))
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
void fixed_event_end_transcription( float *dt,   
//									float tdevelopment,
                                    CellState *state,
                                    GillespieRates *rates,
                                    Genotype *genotype,
                                    int *end_state,                        
                                    int mut_step,
                                    Mutation *mut_record,
//                                    char *effect_of_effector,
                                    float benefit1,
                                    float benefit2,
                                    float max_t_bias)
{
    int gene_id;
    float concurrent;
    float endtime;
    /* recompute the delta-t based on difference between now and the time of transcription end */
    *dt = state->mRNA_transcr_time_end_head->time - state->t;   
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, mut_step, mut_record, benefit1,benefit2,max_t_bias);
    /* get the gene which is ending transcription */
    gene_id = state->mRNA_transcr_time_end_head->event_id;    
    /* increase number of mRNAs that are initializing translation*/
    (state->mRNA_under_transl_delay_num[gene_id])++;
    /* decrease the number of mRNAs undergoing transcription */
    (state->mRNA_under_transc_num[gene_id])--;
    /* delete the fixed even which has just occurred */
    delete_fixed_event_from_head(&(state->mRNA_transcr_time_end_head), &(state->mRNA_transcr_time_end_tail));   
    /*add transcription initialization event*/  
    endtime=state->t+*dt+TRANSLATION_TIME;
    concurrent=check_concurrence(   endtime,
                                    state->mRNA_transl_init_time_end_head,
                                    state->mRNA_transcr_time_end_head,
                                    state->signal_on_head,
                                    state->signal_off_head,
                                    state->burn_in_growth_rate_head,
                                    state->t_to_update_probability_of_binding,
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
                                        state->t_to_update_probability_of_binding,
                                        state->change_signal_strength_head);        
    }  
    /*add to translation initiation event*/
    add_fixed_event(gene_id, endtime, &(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));
}

/*
 * end translation initiation
 */
int fixed_event_end_translation_init(   Genotype *genotype, 
                                        CellState *state, 
                                        GillespieRates *rates, 
                                        float *dt,    
//										float tdevelopment,
                                        int *end_state,                          
                                        int mut_step,
                                        Mutation *mut_record,
//                                        char *effect_of_effector,
                                        float benefit1,
                                        float benefit2,
                                        float max_t_in_bias)
{
    int gene_id;    
    /* calc the remaining time till translation initiation ends */
    *dt = state->mRNA_transl_init_time_end_head->time - state->t;         
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, mut_step, mut_record, benefit1,benefit2, max_t_in_bias);
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

/*
 * compute t' factor used in the integration of growth rate
 * t' is the time the effector protein increases or decreases to saturation level
 */
float calc_tprime(Genotype *genotype, CellState *state, float number_of_selection_protein_bf_dt[NGENES], float dt, float RHS,int gene_ids[NGENES],int start, int n_copies) 
{    
    int i;   
    float protein_synthesis_rate[n_copies],protein_decay_rate[n_copies];

    for(i=start;i<start+n_copies;i++)
    {
        protein_decay_rate[i]=genotype->protein_decay_rate[gene_ids[i]];
        protein_synthesis_rate[i]=protein_decay_rate[i]*state->protein_synthesis_index[gene_ids[i]];
    }       
    return rtsafe(  &calc_fx_dfx,
                    n_copies,
                    RHS,
                    number_of_selection_protein_bf_dt,
                    protein_synthesis_rate,
                    protein_decay_rate,
                    0.0,
                    dt,
                    0.01,
                    start); //rtsafe is in numerical.c
}


/*
 * calculate f(x)-Pp_s and f'(x),
 * f(x) is the number of effector protein molecules at time x
 */
void calc_fx_dfx(float x, int n_copies, float RHS, float *p, float *alpha_s, float *c,float *fx, float *dfx, int start)
{
    int i;    
    *fx=0;
    *dfx=0;    
    for(i=0;i<n_copies;i++)
    {
        *fx+=(p[i+start]-alpha_s[i]/c[i])*exp(-c[i]*x)+alpha_s[i]/c[i];
        *dfx+=(alpha_s[i]-c[i]*p[i+start])*exp(-c[i]*x);
    }
    *fx-=RHS;    
}

/*
 * calculate F(delta_t)/Pp_s. F(x) is the integral of f(x) over delta_t.
 * f(x) is the number of effector protein molecules at time x
 */
float calc_integral(Genotype *genotype, CellState *state, float protein[MAX_OUTPUT_GENES], float dt, float saturate_protein_number, int gene_ids[MAX_OUTPUT_GENES],int start, int n_copies)
{
    int i;
    float integral=0.0,ect_minus_one; 
    
    for(i=start;i<start+n_copies;i++)
    {
        ect_minus_one=exp(-genotype->protein_decay_rate[gene_ids[i]]*dt)-1.0;    
        integral+=1.0/saturate_protein_number*(state->protein_synthesis_index[gene_ids[i]]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]-
                protein[i]*ect_minus_one/genotype->protein_decay_rate[gene_ids[i]]+ 
                state->protein_synthesis_index[gene_ids[i]]*dt);
    }
    return integral;
}

/*
 * return the instantaneous growth rate given the current cell state and environment,
 * also return the integrated growth rate as a pointer
 */
void calc_instantaneous_fitness( Genotype *genotype,
                                CellState *state,
                                float number_of_selection_protein_bf_dt[MAX_OUTPUT_GENES],
                                float number_of_selection_protein_aft_dt[MAX_OUTPUT_GENES],
                                float dt,
//                                float tdevelopment,                               
//                                char effect_of_effector,
                                int *end_state,                         
                                int mut_step,
                                Mutation *mut_record,                                
                                float benefit1,
                                float benefit2,
                                float max_t_bias)
{
    int i,j,protein_id,counter;
    int N_output;
    int gene_copy_id[MAX_OUTPUT_GENES],N_gene_copies[MAX_OUTPUT_PROTEINS];
    float total_translation_rate;    
    float dt_prime, dt_rest, t_on, t_off, t_on_copy,t_off_copy;   
    float penalty=bmax/Ne_saturate;   
    float cost_of_expression;     
    float Ne, Ne_next;   
    int concurrent;
    float endtime;    
                      
    /* compute the total cost of translation across all genes  */
    total_translation_rate=0.0;
#if NO_REGULATION_COST
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {    
        if(genotype->which_protein[i]==genotype->nproteins-1)           
            total_translation_rate += genotype->translation[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation[i]*(float)state->mRNA_under_transl_delay_num[i];
    } 
#else
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {     
        total_translation_rate += genotype->translation_rate[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation_rate[i]*(float)state->mRNA_under_transl_delay_num[i];
    } 
#endif  
    cost_of_expression=total_translation_rate*c_transl;
    state->cumulative_cost+=cost_of_expression*dt;
    /*list the genes that encode effectors*/
    counter=0;
    for(i=0;i<genotype->n_output_proteins;i++)
    {        
        protein_id=genotype->output_protein_id[i];            
        for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)         
        {
            gene_copy_id[counter]=genotype->protein_pool[protein_id][1][j];
            counter++;
        }
    }    
#if POOL_EFFECTORS 
    N_output=1;
    N_gene_copies[0]=genotype->n_output_genes;    
#else
    N_output=genotype->n_output_proteins;
    for(i=0;i<genotype->n_output_proteins;i++)     
        N_gene_copies[i]=genotype->protein_pool[genotype->output_protein_id[i]][0][0];    
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
#elif !PEAK_SEARCH  
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
}


/* 
 * update both the protein concentration and current cell size *
 * 
 */
void update_protein_number_and_fitness( Genotype *genotype,
                                        CellState *state,                                   
                                        GillespieRates *rates,
                                        float dt, 
//                                        float tdevelopment,
//                                        char effect_of_effector,
                                        int *end_state,                                   
                                        int mut_step,
                                        Mutation *mut_record,
                                        float benefit1,
                                        float benefit2,
                                        float max_t_in_bias)
{
    int i,j,protein_id,counter;
    float ct, ect, one_minus_ect;
    float N_output_molecules_bf_dt[MAX_OUTPUT_GENES];  
    float N_output_molecules_aft_dt[MAX_OUTPUT_GENES];  
    FILE *fperror;     
    
    /* store the numbers of the effector proteins encoded by each copy of gene before updating*/  
    counter=0;
    for(i=0;i<genotype->n_output_proteins;i++)
    {        
        protein_id=genotype->output_protein_id[i];            
        for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)         
        {
            N_output_molecules_bf_dt[counter]=state->gene_specific_protein_number[genotype->protein_pool[protein_id][1][j]];
            counter++;
        }
    }
    /* update protein numbers*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {     
        ct = genotype->protein_decay_rate[i]*dt;
        ect = exp(-ct);
        if (fabs(ct)<EPSILON) one_minus_ect=ct;
        else one_minus_ect = 1.0-ect;      
        /* get the new protein concentration for this gene */
        state->gene_specific_protein_number[i]=ect*state->gene_specific_protein_number[i]+state->protein_synthesis_index[i]*one_minus_ect ;        
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
    for(i=0;i<genotype->n_output_proteins;i++)
    {        
        protein_id=genotype->output_protein_id[i];            
        for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)         
        {
            N_output_molecules_aft_dt[counter]=state->gene_specific_protein_number[genotype->protein_pool[protein_id][1][j]];
            counter++;
        }
    }
    /* now find out the protein numbers at end of dt interval and compute growth rate */   
    calc_instantaneous_fitness(genotype, 
                                state, 
                                N_output_molecules_bf_dt, 
                                N_output_molecules_aft_dt,
                                dt, 
//                                tdevelopment,
//                                effect_of_effector,
                                end_state,                                                   
                                mut_step,
                                mut_record,
                                benefit1,
                                benefit2,
                                max_t_in_bias);  
    if(*end_state==0)/*use 0 to indicate abnormal behavior of the program.
                   *I am expecting rounding error to raise this flag.
                   *In case this flag is raised, quit the current replicate
                   *of growth and rerun a replicate.*/
    {
        fperror=fopen(error_file,"a+");
        LOG("at mut step %d",mut_step);
        fclose(fperror);
        return;
    }      
}

/*
 * 
 * Functions that handle each possible Gillespie event 
 *
 */
int Gillespie_event_mRNA_decay(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
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

void Gillespie_event_repressed_to_intermediate(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
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
    state->gene_state[gene_id]=INTERMEDIATE;
}

void Gillespie_event_intermediate_to_repressed(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
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
    state->gene_state[gene_id]=REPRESSED;
}

void Gillespie_event_intermediate_to_active(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
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
    state->gene_state[gene_id] =ACTIVE;  
}

void Gillespie_event_active_to_intermediate(Genotype *genotype, CellState *state,GillespieRates *rates, RngStream RS)
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
    state->gene_state[gene_id]=INTERMEDIATE;   
}

void Gillespie_event_transcription_init(GillespieRates *rates, CellState *state, Genotype *genotype, float dt, RngStream RS)
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
    candidate_t=state->t+dt+TRANSCRIPTION_TIME;
    concurrent=check_concurrence(   candidate_t,
                                    state->mRNA_transl_init_time_end_head,
                                    state->mRNA_transcr_time_end_head,
                                    state->signal_on_head,
                                    state->signal_off_head,
                                    state->burn_in_growth_rate_head,
                                    state->t_to_update_probability_of_binding,
                                    state->change_signal_strength_head);
    while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
    {
        candidate_t+=TIME_OFFSET;
        concurrent=check_concurrence(   candidate_t,
                                        state->mRNA_transl_init_time_end_head,
                                        state->mRNA_transcr_time_end_head,
                                        state->signal_on_head,
                                        state->signal_off_head,
                                        state->burn_in_growth_rate_head,
                                        state->t_to_update_probability_of_binding,                                        
                                        state->change_signal_strength_head);        
    }    
    add_fixed_event(gene_id, candidate_t,&(state->mRNA_transcr_time_end_head), &(state->mRNA_transcr_time_end_tail));
    /* increase the number mRNAs being transcribed */
    (state->mRNA_under_transc_num[gene_id])++;                      
}
/*
 * END
 * Functions that handle each possible Gillespie event 
 */

/*copy genotype from the acestor to offsprings*/
void clone_genotype(Genotype *genotype_templet, Genotype *genotype_clone)
{
    int i, j;           
    /*reset which_cluster for the clone*/
    for(i=0;i<NGENES;i++)
        genotype_clone->which_cluster[i]=-1;
    /*copy which_cluster and cis-reg sequence*/
    for(i=0; i< genotype_templet->ngenes;i++)
    {
        genotype_clone->which_cluster[i]=genotype_templet->which_cluster[i];            
        memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_templet->cisreg_seq[i][0],CISREG_LEN*sizeof(char));                    
        genotype_clone->recalc_TFBS[i]=1;                
    }    
    /*reset clone's cisreg_cluster*/
    i=0;
    while(genotype_clone->cisreg_cluster[i][0]!=-1)
    {
        j=0;
        while(genotype_clone->cisreg_cluster[i][j]!=-1)
        {
            genotype_clone->cisreg_cluster[i][j]=-1;
            j++;
        }
        i++;
    }        
    /*then copy from templet*/
    i=0;
    while(genotype_templet->cisreg_cluster[i][0]!=-1)
    {
        j=0;
        while(genotype_templet->cisreg_cluster[i][j]!=-1)
        {
            genotype_clone->cisreg_cluster[i][j]=genotype_templet->cisreg_cluster[i][j];
            j++;
        }
        i++;
    }
    /*reset clone's information*/
    for(i=0;i<NGENES;i++)
    {
        genotype_clone->which_protein[i]=-1;
        genotype_clone->min_N_activator_to_transc[i]=MAX_BINDING+1;
    }    
    /*reset clone's protein_pool*/
    for(i=0;i<NPROTEINS;i++)
    {            
        for(j=0;j<NGENES;j++)
            genotype_clone->protein_pool[i][1][j]=-1;
        genotype_clone->protein_pool[i][0][0]=0;            
    }
    /*copy from templet's protein_pool*/
    for(i=0;i<genotype_templet->nproteins;i++)
    {            
        genotype_clone->protein_pool[i][0][0]=genotype_templet->protein_pool[i][0][0];            
        for(j=0;j<genotype_templet->protein_pool[i][0][0];j++)
            genotype_clone->protein_pool[i][1][j]=genotype_templet->protein_pool[i][1][j];                     
    }
    /*reset clone's output protein id*/
    for(i=0;i<MAX_OUTPUT_PROTEINS;i++) 
        genotype_clone->output_protein_id[i]=-1;    
    /*then copy*/
    for(i=0;i<genotype_templet->n_output_proteins;i++)
        genotype_clone->output_protein_id[i]=genotype_templet->output_protein_id[i];           
    /* copy binding sites' sequences*/  
    for(i=0; i < genotype_templet->ngenes; i++) 
    {          
        for(j=0;j<TF_ELEMENT_LEN;j++)
        {    
            genotype_clone->tf_binding_seq[i][j]=genotype_templet->tf_binding_seq[i][j];
            genotype_clone->tf_binding_seq_rc[i][j]=genotype_templet->tf_binding_seq_rc[i][j];
        }
    }
    /*copy kinetic constants*/
    for(i=0; i < genotype_templet->ngenes; i++) 
    {            
        genotype_clone->mRNA_decay_rate[i]=genotype_templet->mRNA_decay_rate[i];
        genotype_clone->protein_decay_rate[i]=genotype_templet->protein_decay_rate[i];
        genotype_clone->translation_rate[i]=genotype_templet->translation_rate[i];            
        genotype_clone->active_to_intermediate_rate[i]=genotype_templet->active_to_intermediate_rate[i];
        genotype_clone->which_protein[i]=genotype_templet->which_protein[i];
        genotype_clone->min_N_activator_to_transc[i]=genotype_templet->min_N_activator_to_transc[i];  
    } 
    /* copy TF information*/
    for(i=0;i<NPROTEINS;i++)
    {
        genotype_clone->protein_identity[i][0]=genotype_templet->protein_identity[i][0];
        genotype_clone->protein_identity[i][1]=genotype_templet->protein_identity[i][1];
        genotype_clone->Kd[i]=genotype_templet->Kd[i];
    }    
    /* copy gene and protein numbers*/
    genotype_clone->ngenes=genotype_templet->ngenes; 
    genotype_clone->nproteins=genotype_templet->nproteins;
    genotype_clone->n_output_genes=genotype_templet->n_output_genes;
    genotype_clone->n_output_proteins=genotype_templet->n_output_proteins;
    genotype_clone->N_act=genotype_templet->N_act;
    genotype_clone->N_rep=genotype_templet->N_rep;
    
#if SET_BS_MANUALLY
    for(i=2;i<genotype_templet->ngenes;i++)
    {
        genotype_clone->binding_sites_num[i]=genotype_templet->binding_sites_num[i];
        genotype_clone->max_hindered_sites[i]=genotype_templet->max_hindered_sites[i];
        for(j=0;j<genotype_templet->binding_sites_num[i];j++)
        {
            genotype_clone->all_binding_sites[i][j].tf_id=genotype_templet->all_binding_sites[i][j].tf_id;
            genotype_clone->all_binding_sites[i][j].N_hindered=genotype_templet->all_binding_sites[i][j].N_hindered;
            genotype_clone->all_binding_sites[i][j].Kd=genotype_templet->all_binding_sites[i][j].Kd; 
            genotype_clone->all_binding_sites[i][j].mis_match=genotype_templet->all_binding_sites[i][j].mis_match;
        }
    }
#endif
}

/*
 * run the model for a specified cell for a single timestep:
 */
void do_single_timestep(Genotype *genotype, 
                        CellState *state,                         
                        GillespieRates *rates,                        
                        float signal_strength,
                        float response_amplification,
                        float minimal_peak_response,
                        float benefit1,
                        float benefit2,
                        float max_t_in_bias,    
                        float duration_signal_on,
                        float duration_signal_off,                    
                        RngStream RS,                       
                        int mut_step,
                        Mutation *mut_record,
                        int *end_state,                
                        int thread_id,
                        float tdevelopment,
                        float *signal_profile) 
{    
    int event, UPDATE_WHAT;     
    float fixed_time; 
    float dt;
    float x;
    FILE *fperror;
//    char *effect_of_effector;

    x = expdev(RS);        /* draw random number */
    dt = x/rates->total_Gillespie_rate;
    if (dt < 0.0) 
    {	
        fperror=fopen(error_file,"a+");
        LOG("negative dt at mut_step %d\n",mut_step);
        fclose(fperror);
        *end_state=0; /*use 0 to indicate abnormal behavior of the program.
                       *I am expecting rounding error to raise this flag.
                       *In case this flag is raised, quit the current replicate
                       *of growth and rerun a replicate.*/
        return;    
    }
    /* first check to see if a fixed event occurs in current t->dt window, or in tdevelopment if running for a fixed development time */
    fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment;
    event = does_fixed_event_end(state, fixed_time);
    while(event!=0)
    {           
        /*after doing fixed event, return a flag to indicate whether mandatorily update Pact or Prep*/
        UPDATE_WHAT=do_fixed_event( genotype, 
                                    state, 
                                    rates, 
                                    &dt, 
//									tdevelopment,
                                    event, 
                                    signal_strength,
                                    duration_signal_on, 
                                    duration_signal_off, 
                                    response_amplification,
                                    minimal_peak_response,
                                    benefit1,
                                    benefit2,
                                    max_t_in_bias,
//                                    effect_of_effector,
                                    end_state, 
                                    mut_step, 
                                    mut_record,
                                    signal_profile);        
        if(*end_state==0) 
            return; 
        state->t += dt;  /* advance time by the dt */        
#if CAUTIOUS //watch for overlapping fixed events -- one of which will be skipped and cause error
        float dt_copy,x_copy,rate_subtotal_copy,t_copy;  
        int event2;      
        fixed_time=state->t+dt;
        event2=does_fixed_event_end(state, fixed_time);
        if(event2!=0)
        {
            fperror=fopen(error_file,"a+");
            LOG("overlapping fixed events at mut_step %d", mut_step);
            fclose(fp);
        }                       
        t_copy=state->t; 
        dt_copy=dt;
        x_copy=x; 
        rate_subtotal_copy=rates->total_Gillespie_rate; 
#endif        
        x -= dt*rates->total_Gillespie_rate;  /* we've been running with rates->total_Gillespie_rate for dt, so substract it from x*/        
        calc_all_rates(genotype, state, rates,tdevelopment,UPDATE_WHAT,thread_id);  /* update rates->total_Gillespie_rate and re-compute a new dt */      
        dt = x/rates->total_Gillespie_rate;
        /*deal with rounding error*/
        if(dt<0.0)
        {  	
#if CAUTIOUS // this rounding error can happen very often, therefore the error_log can be huge
            fperror=fopen(error_file,"a+");
            LOG("rounding error in dt at mut_step %d\n",mut_step);            	
            fclose(fp);
#endif 
            dt=TIME_OFFSET; 
        }
        fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment; 
        /* check to see there aren't more fixed events to do */
        event = does_fixed_event_end(state, fixed_time);                                    
    } 
  /* no remaining fixed events to do in dt, now do stochastic events */  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (state->t+dt < tdevelopment)
    {
        update_protein_number_and_fitness(genotype, state, rates, dt, end_state, mut_step, mut_record,benefit1, benefit2, max_t_in_bias);
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, RS, end_state, mut_step, mut_record);
        if(*end_state==0)
            return;  
        /* Gillespie step: advance time to next event at dt */
        state->t += dt;
        calc_all_rates(genotype,state,rates,tdevelopment,UPDATE_WHAT,thread_id);        
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - state->t;
        /* final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, dt, end_state, mut_step, mut_record, benefit1, benefit2, max_t_in_bias);
        if(*end_state==0)
            return; 
        state->t+=dt;
#if SELECT_SENSITIVITY_AND_PRECISION
        state->t = tdevelopment;
		state->sensitivity[2] = state->sensitivity[0];
		state->precision[2] = state->precision[0];
#endif
        /* advance to end of development (this exits the outer while loop) */   
    }
}

/* while there are fixed events
     occuring in current t->dt window */
int do_fixed_event(Genotype *genotype, 
                    CellState *state, 
                    GillespieRates *rates, 
                    float *dt,  
//					float tdevelopment,
                    int event,
                    float signal_on_strength,
                    float duration_signal_on,
                    float duration_signal_off,     
                    float response_amplification,
                    float minimal_peak_response,
                    float benefit1,
                    float benefit2,
                    float max_t_in_bias,
//                    char *effect_of_effector,
                    int *end_state,                  
                    int mut_step,
                    Mutation *mut_record,
                    float *signal_profile)
{     
    int i,j, protein_id, return_value;
    return_value=DO_NOTHING;
    switch (event) 
    {
        case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, state, rates, genotype,end_state,mut_step,mut_record,benefit1,benefit2, max_t_in_bias);
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, dt, end_state, mut_step,mut_record,benefit1,benefit2, max_t_in_bias);
            state->cell_activated=1;
            break;
        case 3:     /* turn signal off*/ 
            *dt = state->signal_off_head->time - state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, mut_step, mut_record, benefit1,benefit2, max_t_in_bias);
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));
//            if(fixed_effector_effect)
//                *effect_of_effector=init_effector_effect;
//            else
//                *effect_of_effector='d';
#if SELECT_SENSITIVITY_AND_PRECISION
			state->precision[0] = state->protein_number[genotype->nproteins - 1]; //reset running precision
#endif
            break;
        case 4:     /*turn signal on*/
            *dt = state->signal_on_head->time - state->t;   
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, mut_step, mut_record, benefit1,benefit2, max_t_in_bias);
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
            state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;
            state->recording_basal_response=0;
            state->found_gradient=1;
#if POOL_EFFECTORS
            state->threshold_response_to_bias[0]=(minimal_peak_response>state->basal_response[0]*response_amplification)?minimal_peak_response:state->basal_response[0]*response_amplification;
#endif           
//            if(fixed_effector_effect)                               
//                *effect_of_effector=init_effector_effect;            
//            else                 
//                *effect_of_effector='b'; 
#if SELECT_SENSITIVITY_AND_PRECISION
			state->precision[1] = state->precision[0]; //record the max precision to signal change 1
			state->sensitivity[1] = state->sensitivity[0]; //record the max sensitivity to signal change 1
			state->sensitivity[0] = 0.0;	// reset running sensitivity
#endif
            break;	
        case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, mut_step, mut_record,benefit1,benefit2, max_t_in_bias);
            for(i=0;i<MAX_OUTPUT_PROTEINS;i++)
                state->cumulative_fitness_after_burn_in[i]=state->cumulative_fitness[i];           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6: /* mandatorily updating Pact and Prep*/
            *dt=state->t_to_update_probability_of_binding-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, mut_step, mut_record, benefit1,benefit2, max_t_in_bias);
            break;
        case 7: /* update signal strength */
            *dt=state->change_signal_strength_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, mut_step, mut_record, benefit1,benefit2, max_t_in_bias);
            state->protein_number[N_SIGNAL_TF-1]=signal_profile[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            break;  
        case 8: /* sample concentration*/           
            *dt=state->sampling_point_end_head->time-state->t;
#if CHEMOTAXIS            
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL, benefit1,benefit2, max_t_in_bias);
            state->recording_basal_response=1;
            delete_fixed_event_from_head(&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
#else
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL, benefit1,benefit2, max_t_in_bias);
            delete_fixed_event_from_head(&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
            state->sampled_response[state->N_samples]=0.0;
            for(i=0;i<genotype->n_output_proteins;i++)
            {                
                protein_id=genotype->output_protein_id[i];            
                for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
                    state->sampled_response[state->N_samples]+=state->gene_specific_protein_number[genotype->protein_pool[protein_id][1][j]];
            } 
            state->N_samples++;
#endif           
            break;          
    }  
    return return_value;
}

/*
 *Gillespie events
 */
int do_Gillespie_event(Genotype *genotype,
                        CellState *state,
                        GillespieRates *rates,
                        float dt,                    
                        RngStream RS,
                        int *end_state,                        
                        int mut_step,
                        Mutation *mut_record)
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

/*Free linked tables*/
void free_fixedevent(CellState *state)
{
    FixedEvent *temp1, *temp2;
    /*signal_on_starts_end*/
    temp1=state->signal_on_head;
    while(temp1){		
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signal_on_head=NULL;
    state->signal_on_tail=NULL;
    /*signal_off_starts_end*/
    temp1=state->signal_off_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signal_off_head=NULL;
    state->signal_off_tail=NULL;
    /*mRNA_transcr_time_end*/
    temp1=state->mRNA_transcr_time_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;           
            free(temp2);	
    }
    state->mRNA_transcr_time_end_head=NULL;
    state->mRNA_transcr_time_end_tail=NULL;
    /*mRNA_transl_init_time_end*/
    temp1=state->mRNA_transl_init_time_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->mRNA_transl_init_time_end_head=NULL;
    state->mRNA_transl_init_time_end_tail=NULL;
    /*burn_in_growth_rate_end*/
    temp1=state->burn_in_growth_rate_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->burn_in_growth_rate_head=NULL;
    state->burn_in_growth_rate_tail=NULL;
    /*change_singal_strength_head*/
    temp1=state->change_signal_strength_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->change_signal_strength_head=NULL;
    state->change_signal_strength_tail=NULL;
    /*sampling_point_end_head*/
    temp1=state->sampling_point_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->sampling_point_end_head=NULL;
    state->sampling_point_end_tail=NULL;
}

/**
 *Calculate the fintess of a given genotype.
 *Essentially calling do_single_timestep until tdevelopment and calculate 
 *average growth rate over tdevelopment.
 */
void calc_cellular_fitness(  Genotype *genotype,
                            int init_mRNA[NGENES],
                            float init_protein_number[NPROTEINS],
                            RngStream RS_parallel[N_THREADS],                           
                            int mut_step,
                            float GR1[N_REPLICATES],
                            float GR2[N_REPLICATES],
                            Mutation *mut_record)       
{ 
    /*Making clones of a genotype, and have the clones run in parallel*/
    #pragma omp parallel num_threads(N_THREADS) 
    {
        int thread_ID=omp_get_thread_num();
//        int thread_ID=0;
        int i,j,k;
        int N_replicates_per_thread=N_REPLICATES/N_THREADS;
        int end_state;
        char effect_of_effector;
        Genotype genotype_clone;
        CellState state_clone;
        GillespieRates rate_clone;
        int init_mRNA_clone[NGENES]; 
        float init_protein_number_clone[NGENES];
        float gr1[N_replicates_per_thread],gr2[N_replicates_per_thread];   
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];     
        float *signal_profile;
#if CAUTIOUS
        FILE *fperror;
#endif
        /*alloc space for linked tables and set default values for parameters, in genotype*/
        initialize_cache(&genotype_clone);  
        /*initialize growth rate to 0*/
        for(i=0;i<N_replicates_per_thread;i++)
        {
            gr1[i]=0.0;
            gr2[i]=0.0;
        }
        /*clone genotype and initial mRNA and protein numbers*/
       #pragma omp critical
        { 
            genotype_clone.ngenes=genotype->ngenes;          
            genotype_clone.nproteins=genotype->nproteins;
            clone_genotype(genotype, &genotype_clone);
#if FORCE_OR_GATE           
            for(i=0;i<11;i++)
                genotype_clone.N_motifs[i]=genotype->N_motifs[i];
            for(i=0;i<NGENES;i++)
            {
                genotype_clone.gene_in_core_C1ffl[i]=genotype->gene_in_core_C1ffl[i];
                for(j=0;j<NPROTEINS;j++)
                    genotype_clone.TF_in_core_C1ffl[i][j]=genotype->TF_in_core_C1ffl[i][j];
            }
#endif  
            for(j=0; j < NGENES; j++) 
            {  
                init_mRNA_clone[j] = init_mRNA[j];
                init_protein_number_clone[j] = init_protein_number[j];
            }
        #if SET_BS_MANUALLY
            for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            {
                genotype_clone.binding_sites_num[i]=genotype->binding_sites_num[i];
                genotype_clone.max_hindered_sites[i]=genotype->max_hindered_sites[i];
                for(j=0;j<genotype->binding_sites_num[i];j++)
                {
                    genotype_clone.all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i][j].tf_id;
                    genotype_clone.all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i][j].N_hindered;
                    genotype_clone.all_binding_sites[i][j].Kd=genotype->all_binding_sites[i][j].Kd;                    
                }
            }   
        #endif
        } 
        /*Set initial mRNA and protein number using given values*/
        for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)  //expression of the sensor TF is not considered in the model           
            mRNA[j] = init_mRNA_clone[j];                       
        for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
        {
            for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0]; //split the initial protein number equally to different copies
                                                                                                                                        //this is to make sure all proteins have equal initial numbers
        }        
    #if !SET_BS_MANUALLY
        calc_all_binding_sites(&genotype_clone);
    #endif        
        /* now calc growth rate under two environments*/
        for(i=0;i<N_replicates_per_thread;i++) /* env 1, usually a constant signal that matches env*/
        {	 
            effect_of_effector=env1_initial_effect_of_effector; // whether the effector is beneficial or deleterious
            end_state=1; //flag to show whether do_single_time_step encounters an error. 1 means no error, 0 means error           
#if EXTERNAL_SIGNAL
            j=RngStream_RandInt(RS_parallel[thread_ID],0,99);
            signal_profile=&(signal_profile_matrix[thread_ID][j][0]);
#else
            signal_profile=NULL;            
#endif            
            /*initialize mRNA and protein numbers, and gene states, in cell*/
            initialize_cell(&genotype_clone,
                            &state_clone,                             
                            genotype_clone.ngenes, 
                            genotype_clone.nproteins,
                            genotype_clone.protein_pool,
                            genotype_clone.mRNA_decay_rate, 
                            mRNA, 
                            protein, 
                            RS_parallel[thread_ID],
                            env1_t_development,
                            env1_t_signal_off);
            /*set how the environment signal should change during simulation*/
            set_signal(&state_clone,
                        env1_t_signal_on,
                        env1_t_signal_off,
                        signal_profile,
                        env1_t_development,
                        env1_signal1_strength);          
            /*growth starts at 0 min*/
            state_clone.t = 0.0;
            /*calcualte the rates of cellular activity based on the initial cellular state*/
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, env1_t_development,INITIALIZATION,thread_ID); 
            /*run growth simulation until tdevelopment or encounter an error*/
            while(state_clone.t<env1_t_development && end_state==1) 
            {
                do_single_timestep(&genotype_clone, 
                                    &state_clone, 
                                    &rate_clone,                               
                                    env1_signal2_strength,
                                    env1_response_amplification,
                                    env1_minimal_peak_response,
                                    env1_benefit1,
                                    env1_benefit2,
                                    env1_max_duration_bias,
                                    env1_t_signal_on, 
                                    env1_t_signal_off,
                                    RS_parallel[thread_ID],                                   
                                    mut_step,  
                                    mut_record,
                                    &end_state,                           
                                    thread_ID,
                                    env1_t_development,
                                    signal_profile);
            } 
            if(end_state==1) // no error
            {               
                gr1[i]=calc_replicate_fitness(&state_clone,1,env1_t_development,env1_t_signal_off,env1_response_amplification,env1_minimal_peak_response,genotype_clone.n_output_proteins);
            #if CAUTIOUS                
                if(gr1[i]<0.0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("negative growth rate at replicate %d in test 1 thread %d at mut step %d\n", i, thread_ID, mut_step);
                    fclose(fp);
                }              
            #endif
            }
            else
            {                
            	/* if do_single_timestep throws out an error. Re-run replicates*/
                i--;
#if CAUTIOUS
                fperror=fopen(error_file,"a+");
                LOG("Rerun replicates at replicate %d in test 1 thread %d at mut step %d\n", i, thread_ID, mut_step);
                fclose(fp);
#endif
            }
            /*free linked tables*/
#if PEAK_SEARCH
            free(state_clone.sampled_response);
#endif
            free_fixedevent(&state_clone);           
        }            
        /*******************env2*********************/
        for(i=0;i<N_replicates_per_thread;i++) 
        {	 
            effect_of_effector=env2_initial_effect_of_effector;  
            end_state=1;           
            signal_profile=NULL;            
            initialize_cell(&genotype_clone,
                            &state_clone, 
                            genotype_clone.ngenes, 
                            genotype_clone.nproteins,
                            genotype_clone.protein_pool,
                            genotype_clone.mRNA_decay_rate, 
                            mRNA, 
                            protein, 
                            RS_parallel[thread_ID],
                            env2_t_development,
                            env2_t_signal_off);
            set_signal( &state_clone,
                        env2_t_signal_on,   
                        env2_t_signal_off,
                        signal_profile,
                        env2_t_development,
                        env2_signal1_strength);  
            state_clone.t = 0.0;
            calc_all_rates( &genotype_clone, 
                            &state_clone, 
                            &rate_clone,
                            env2_t_development,
                            INITIALIZATION,thread_ID);	
            while(state_clone.t<env2_t_development && end_state==1)
            {
                do_single_timestep(&genotype_clone, 
                                    &state_clone, 
                                    &rate_clone,               
                                    env2_signal2_strength,
                                    env2_response_amplification,
                                    env2_minimal_peak_response,
                                    env2_benefit1,
                                    env2_benefit2,
                                    env2_max_duration_bias,
                                    env2_t_signal_on,
                                    env2_t_signal_off,
                                    RS_parallel[thread_ID],                                  
                                    mut_step,
                                    mut_record,
                                    &end_state,                           
                                    thread_ID,
                                    env2_t_development,
                                    signal_profile);
            } 
            if(end_state==1)
            {
                gr2[i]=calc_replicate_fitness(&state_clone,2,env2_t_development,env2_t_signal_off,env2_response_amplification,env2_minimal_peak_response,genotype_clone.n_output_proteins);
#if CAUTIOUS             
                if(gr2[i]<0.0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("negative growth rate at replicate %d in test 2 thread %d at mut step %d\n", i, thread_ID, mut_step);
                    fclose(fp);
                }               
#endif
            }
            else
            {                
                i--;
#if CAUTIOUS
                fperror=fopen(error_file,"a+");
                LOG("Rerun replicates at replicate %d in test 2 thread %d at mut step %d\n", i, thread_ID, mut_step);
                fclose(fp);
#endif
            } 
#if PEAK_SEARCH
            free(state_clone.sampled_response);
#endif
            free_fixedevent(&state_clone);            
        }    
        /*free linked tables*/
        for(j=0;j<NGENES;j++)
            free(genotype_clone.all_binding_sites[j]);
        /*pool growth rates from each thread*/
       #pragma omp critical
        {
            j=0;
            for(i=thread_ID*N_replicates_per_thread;i<(thread_ID+1)*N_replicates_per_thread;i++)
            {
                GR1[i]=gr1[j];
                GR2[i]=gr2[j];
                j++;
            }
        }        
    } 
}

float calc_replicate_fitness(CellState *state, int which_env, float t_development, float t_signal_on, float response_amplitude, float minimal_peak_response, int N_output)
{
    float fitness, peak_response;
    int i;
    float max_response,fifty_percent_response,ninty_percent_recovery;
    int pos_max_response;
    float pos_fifty_percent_response,pos_ninty_percent_recovery;
    FILE *fp;
#if SELECT_SENSITIVITY_AND_PRECISION
    fitness = (state->sensitivity[which_env] > Ne_saturate) ? 1.0 : (state->sensitivity[which_env] / Ne_saturate);
    fitness -= (state->precision[which_env] > Ne_saturate) ? 0.5 : 0.5*(state->precision[1which_env] / Ne_saturate);
    if(t_development>100.0)
    {
        fitness += (state->sensitivity[which_env] > Ne_saturate) ? 1.0 : (state->sensitivity[which_env] / Ne_saturate);
        fitness -= (state->precision[which_env] > Ne_saturate) ? 0.5 : 0.5*(state->precision[which_env] / Ne_saturate);
    }
    fitness -= state->cumulative_cost / env1_t_development;
#elif SELECT_ON_DURATION
    if(state->Duration_pulse<=opt_pulse_duration)    // fitness is proportional to pulse duration            
        fitness=state->Duration_pulse/opt_pulse_duration;                
    else                
        fitness=1.0-0.5*(state->Duration_pulse-opt_pulse_duration)/opt_pulse_duration; //but is penalized for overly long duration  
    fitness-=0.5*state->T_pulse_on/opt_pulse_duration;
//                gr2[i]+=(state->T_pulse_on<tolerable_delay_bf_pulse)?0.5*(tolerable_delay_bf_pulse-state->T_pulse_on)/tolerable_delay_bf_pulse:0.0; // and is also penalized by late pulse
    fitness -= state->cumulative_cost / t_development;
#elif REALLY_COMPLETECATE
    fitness=0.1*bmax*state->cumulative_basal_benefit/saturate_cumulative_response_from_pulse-state->cumulative_cost/t_development;
    fitness+=(state->cumulative_advanced_benefit<saturate_cumulative_response_from_pulse)?state->cumulative_advanced_benefit/saturate_cumulative_response_from_pulse*bmax:bmax;
    fitness-=0.5*bmax*state->cumulative_damage/saturate_cumulative_response_from_pulse;
#elif PEAK_SEARCH
    find_max(&(state->sampled_response[0]),0,state->N_samples,&max_response,&pos_max_response);
    if(max_response==0.0 || fabs(max_response-state->sampled_response[0])/state->sampled_response[0]<=EPSILON)//"flat" or monotonous decrease
    {
        fitness=0.0-state->cumulative_cost/t_development;
    }
    else 
    {
        if(pos_max_response==state->N_samples-1)//monotonous increase
        {           
            pos_ninty_percent_recovery=TIME_INFINITY;
        }    
        else //there is a pulse. Flat tail is taken care of by find_x, making pos_mid2=TIME_INFINITY
        {  
            ninty_percent_recovery=max_response*0.1+state->sampled_response[0]*0.9;  
            find_x(&(state->sampled_response[0]),pos_max_response,state->N_samples-1,ninty_percent_recovery,&pos_ninty_percent_recovery,1); // look for mid point aft the peak
        }
        fifty_percent_response=(max_response+state->sampled_response[0])*0.5;
        find_x(&(state->sampled_response[0]),0,pos_max_response,fifty_percent_response,&pos_fifty_percent_response,0); // look for mid point bf the peak
        peak_response=(minimal_peak_response>(response_amplitude*state->sampled_response[0]))?minimal_peak_response:(response_amplitude*state->sampled_response[0]);
//        fitness=(max>saturate_pulse_amplitude)?1.0:max/saturate_pulse_amplitude;
        fitness=max_response/(max_response+peak_response/9.0);
        fitness*=(t_development-t_signal_on-pos_fifty_percent_response)/(t_development-t_signal_on);
//        fitness*=exp(-pow((sampling_interval*(pos_mid2-pos_mid1)-opt_pulse_duration)/sd_opt_pulse_duration,2.0));
        fitness*=opt_pulse_duration*9.0/(opt_pulse_duration*9.0+sampling_interval*(pos_ninty_percent_recovery-pos_fifty_percent_response));
        fitness=fitness-state->cumulative_cost/t_development;          
    }             
#elif CHEMOTAXIS
    if(POOL_EFFECTORS)
        fitness=(state->cumulative_benefit-state->cumulative_cost)/t_development;
#else    
    if(POOL_EFFECTORS)
        fitness=(state->cumulative_fitness[0]-state->cumulative_fitness_after_burn_in[0])/(t_development-duration_of_burn_in_growth_rate); 
    else
    {
        fitness=0.0;
        for(i=0;i<N_output;i++)
            fitness+=(state->cumulative_fitness[i]-state->cumulative_fitness_after_burn_in[i])/(t_development-duration_of_burn_in_growth_rate); 
        fitness/=N_output;
    }
    fitness-=state->cumulative_cost/t_development;
#endif
    return fitness;
}


/*To get protein number and fitness at given timepoint, we add a special fixed event (
 *sampling_point_end)to take a snapshot of cell state at given timepoint. The additional fixed event 
 *requires some revision to do_singe_timestep etc. The functions that require revision
 *are compiled only if necessary.
 */
#if JUST_PLOTTING 
int does_fixed_event_end_plotting(  CellState *state,
                                    float t) 
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
    t6 = state->sampling_point_end_head ? state->sampling_point_end_head->time : TIME_INFINITY;
    t7 = state->t_to_update_probability_of_binding;
    t8 = state->change_signal_strength_head ? state->change_signal_strength_head->time : TIME_INFINITY;
    if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5) && (t1<=t6) && (t1<=t7)&&(t1<=t8))
    {
        retval = 1;	
    }
    else if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5) && (t2<=t6) && (t2<=t7)&&(t2<=t8)) 
    {
        retval = 2;
    }            
    else if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5) && (t3<=t6) && (t3<=t7)&&(t3<=t8)) 
    {
        retval = 3;
    }
    else if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5) && (t4<=t6) && (t4<=t7)&&(t4<=t8)) 
    {
        retval = 4;
    }                    
    else if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4) && (t5<=t6) && (t5<=t7)&&(t5<=t8))
    {
        retval = 5;
    }
    else if((t6 <= t1) && (t6 <= t) && (t6 <= t2) && (t6 <= t3) && (t6<=t4) && (t6<=t5) && (t6<=t7)&&(t6<=t8))
    {
        retval=6;
    }
    else if((t7 <= t1) && (t7 <= t) && (t7 <= t2) && (t7 <= t3) && (t7<=t4) && (t7<=t5) && (t7<=t6)&&(t7<=t8))
    {
        retval=7;
    }
    else if((t8 <= t1) && (t8 <= t) && (t8 <= t2) && (t8 <= t3) && (t8<=t4) && (t8<=t5) && (t8<=t6)&&(t8<=t7))            
    {
        retval=8;
    }
    else    
    {
        retval = 0;
    }                       
    return retval;
}

int do_fixed_event_plotting(    Genotype *genotype, 
                                CellState *state, 
                                GillespieRates *rates, 
                                float *dt,                                 
                                int event, 
                                float signal_on_strength,
                                float duration_signal_on,
                                float duration_signal_off,     
                                float response_amplification,
                                float minimal_peak_response,
                                float benefit1,
                                float benefit2,
                                float max_t_in_bias,                             
                                int *end_state,
                                float *signal_profile)
{      
    int i, j, return_value,protein_id;
    return_value=0;
    switch (event) 
    {        
        case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, state, rates, genotype,end_state,0,NULL,benefit1,benefit2, max_t_in_bias);
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, dt, end_state, 0,NULL,benefit1,benefit2, max_t_in_bias);
            break;
        case 3:     /* turn signal off*/ 
            *dt = state->signal_off_head->time - state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL, benefit1,benefit2, max_t_in_bias);
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));
//            if(fixed_effector_effect)
//                *effect_of_effector=init_effector_effect;
//            else
//                *effect_of_effector='d';
#if SELECT_SENSITIVITY_AND_PRECISION
			state->precision[0] = state->protein_number[genotype->nproteins - 1]; //reset running precision
#endif
            break;
        case 4:     /*turn signal on*/
            *dt = state->signal_on_head->time - state->t;   
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL, benefit1,benefit2, max_t_in_bias);
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
            if(state->protein_number[N_SIGNAL_TF-1]<signal_on_strength)
                state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;
            else
                state->protein_number[N_SIGNAL_TF-1]*=2.0;
            state->recording_basal_response=0;
            state->found_gradient=1;
#if POOL_EFFECTORS
            state->threshold_response_to_bias[0]=(minimal_peak_response>state->basal_response[0]*response_amplification)?minimal_peak_response:state->basal_response[0]*response_amplification;
#endif           
//            if(fixed_effector_effect)                               
//                *effect_of_effector=init_effector_effect;            
//            else                 
//                *effect_of_effector='b'; 
#if SELECT_SENSITIVITY_AND_PRECISION
			state->precision[1] = state->precision[0]; //record the max precision to signal change 1
			state->sensitivity[1] = state->sensitivity[0]; //record the max sensitivity to signal change 1
			state->sensitivity[0] = 0.0;	// reset running sensitivity
#endif
            break;	
        case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-state->t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL,benefit1,benefit2, max_t_in_bias);
            for(i=0;i<MAX_OUTPUT_PROTEINS;i++)
                state->cumulative_fitness_after_burn_in[i]=state->cumulative_fitness[i];           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6:
            *dt=state->sampling_point_end_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL, benefit1,benefit2, max_t_in_bias);
            delete_fixed_event_from_head(&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
            if(state->t>60.0)
            {
                for(i=0;i<genotype->n_output_proteins;i++)
                {                
                    protein_id=genotype->output_protein_id[i];            
                    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++) 
                        state->sampled_response[state->N_samples]+=state->gene_specific_protein_number[genotype->protein_pool[protein_id][1][j]];
                } 
                state->N_samples++;  
            }
            
            break;
        case 7: /* update force to update Pact and Prep*/
            *dt=state->t_to_update_probability_of_binding-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL, benefit1,benefit2, max_t_in_bias);
            break;  
        case 8: /* update signal strength */
            *dt=state->change_signal_strength_head->time-state->t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, end_state, 0, NULL, benefit1,benefit2, max_t_in_bias);
            state->protein_number[N_SIGNAL_TF-1]=signal_profile[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            break;
    }    
    return return_value;
}

void do_single_timestep_plotting(   Genotype *genotype, 
                                    CellState *state,                         
                                    GillespieRates *rates,                                                           
                                    float signal_strength,
                                    float response_amplification,
                                    float minimal_peak_response,
                                    float benefit1,
                                    float benefit2,
                                    float max_t_in_bias,    
                                    float duration_signal_on,
                                    float duration_signal_off,       
                                    float (*phenotype)[N_TIMEPOINTS],                                    
                                    float fitness[N_TIMEPOINTS],                                    
                                    RngStream RS,                                    
                                    int *timepoint,                                   
                                    int *end_state,
                                    int thread_id,
                                    float tdevelopment,
                                    float *signal_profile)
{    
    int event,i,UPDATE_WHAT;      
    float fixed_time; 
    float dt;
    float x;    
    x = expdev(RS);        
    dt = x/rates->total_Gillespie_rate;    
    fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment;
    event = does_fixed_event_end_plotting(state, fixed_time);
    while(event!=0)
    {           
        UPDATE_WHAT=do_fixed_event_plotting(genotype, 
                                            state, 
                                            rates, 
                                            &dt, 
                                            event, 
                                            signal_strength,
                                            duration_signal_on, 
                                            duration_signal_off, 
                                            response_amplification,
                                            minimal_peak_response,
                                            benefit1,
                                            benefit2,
                                            max_t_in_bias,
                                            end_state,
                                            signal_profile);
        if(*end_state==0)
            return;        
        if(event==6)/*time to take a snapshot of fitness and protein number*/
        {    
            for(i=0;i<genotype->nproteins;i++)
                phenotype[i][*timepoint]=state->protein_number[i];            
//            fitness[*timepoint]=state->instantaneous_fitness;          
            (*timepoint)++;    
        }              
        fixed_time=state->t+dt; 
        state->t += dt;                     
        x -= dt*rates->total_Gillespie_rate;       
        calc_all_rates(genotype, state, rates, tdevelopment,UPDATE_WHAT,thread_id);  
        dt = x/rates->total_Gillespie_rate;	
        if(dt<0.0)	
            dt=TIME_OFFSET;         
        /* check to see there aren't more fixed events to do */
        fixed_time = (state->t+dt<tdevelopment)?(state->t+dt):tdevelopment;         
        event = does_fixed_event_end_plotting(state, fixed_time);                                    
    } 
    /* no remaining fixed events to do in dt, now do stochastic events */  
    /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (state->t+dt < tdevelopment)
    { 
        update_protein_number_and_fitness(genotype, state, rates, dt, end_state, 0, NULL, benefit1, benefit2, max_t_in_bias);  
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, RS, end_state, 0, NULL);
        if(*end_state==0)
            return;  
        /* Gillespie step: advance time to next event at dt */
        state->t += dt;
        calc_all_rates(genotype,state,rates,tdevelopment,UPDATE_WHAT,thread_id);
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - state->t;
        /* final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, dt, end_state, 0, NULL, benefit1, benefit2, max_t_in_bias);
       
        if(*end_state==0)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        state->t = tdevelopment;
    }
}


void calc_cellular_fitness_plotting( Genotype *genotype,
                                        int init_mRNA[NGENES],
                                        float init_protein_number[NPROTEINS],
                                        RngStream RS[N_THREADS])
{     
    float phenotypeA[N_REPLICATES][genotype->nproteins][N_TIMEPOINTS];
    float phenotypeB[N_REPLICATES][genotype->nproteins][N_TIMEPOINTS];
    float fitnessA[N_REPLICATES][N_TIMEPOINTS];
    float fitnessB[N_REPLICATES][N_TIMEPOINTS];
    float fitnessA2[N_REPLICATES];
    float fitnessB2[N_REPLICATES];
    int l,m,n;
    char filename1[32],filename2[32];
    FILE *fp1,*fp2;
    
    for(l=0;l<N_REPLICATES;l++)
    {
        for(m=0;m<genotype->nproteins;m++)
        {
            for(n=0;n<N_TIMEPOINTS;n++)
            {
                phenotypeA[l][m][n]=0.0;
                phenotypeB[l][m][n]=0.0;               
            }
        }
    }    
    for(l=0;l<N_REPLICATES;l++)
    {   
        fitnessA2[l]=0.0;
        fitnessB2[l]=0.0;
        for(n=0;n<N_TIMEPOINTS;n++)
            {
                fitnessA[l][n]=0.0;
                fitnessB[l][n]=0.0;               
            }        
    }
   
    #pragma omp parallel num_threads(N_THREADS)
    {
        int thread_ID=omp_get_thread_num();
//        int thread_ID=0;
        int i,j,timepoint,k;    
        int end_state;  
        int N_replicates_per_thread=N_REPLICATES/N_THREADS;
//        int N_replicates_per_thread=40;
        char effect_of_effector;
        Genotype genotype_clone;
        CellState state_clone;
        GillespieRates rate_clone;
        float init_mRNA_clone[NGENES]; 
        float init_protein_number_clone[NGENES];
        float t;
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];
        float *signal_profile;
        initialize_cache(&genotype_clone);
        #pragma omp critical
        {  
            clone_genotype(genotype, &genotype_clone);
            for(j=0; j < NGENES; j++) 
            {  
                init_mRNA_clone[j] = init_mRNA[j];
                init_protein_number_clone[j] = init_protein_number[j];
            }
        #if SET_BS_MANUALLY
            for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            {
                genotype_clone.binding_sites_num[i]=genotype->binding_sites_num[i];
                genotype_clone.max_hindered_sites[i]=genotype->max_hindered_sites[i];
                for(j=0;j<genotype->binding_sites_num[i];j++)
                {
                    genotype_clone.all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i][j].tf_id;
                    genotype_clone.all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i][j].N_hindered;
                    genotype_clone.all_binding_sites[i][j].Kd=genotype->all_binding_sites[i][j].Kd;                    
                }
            }   
        #endif
        }        
    #if !SET_BS_MANUALLY        
        calc_all_binding_sites(&genotype_clone); 
    #endif   
        for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)             
                mRNA[j] = init_mRNA_clone[j];                       
        for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
        {
            for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0];                
        }  
        
        for(i=0;i<N_replicates_per_thread;i++)        
        {             
            effect_of_effector=env1_initial_effect_of_effector;
            end_state=1;
            signal_profile=NULL;
                     
            initialize_cell(&genotype_clone,
                            &state_clone,                             
                            genotype_clone.ngenes, 
                            genotype_clone.nproteins,
                            genotype_clone.protein_pool,
                            genotype_clone.mRNA_decay_rate, 
                            mRNA, 
                            protein, 
                            RS[thread_ID],
                            env1_t_development,
                            env1_t_signal_off);            
            set_signal(&state_clone,
                        env1_t_signal_on,
                        env1_t_signal_off,
                        signal_profile,
                        env1_t_development,
                        env1_signal1_strength);         
            state_clone.t = 0.0;
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, env1_t_development,INITIALIZATION,thread_ID);
            timepoint=0;
            while(state_clone.t<env1_t_development && end_state==1)
            { 
                do_single_timestep_plotting(&genotype_clone, 
                                            &state_clone, 
                                            &rate_clone,                                                                                                                                
                                            env1_signal2_strength,
                                            env1_response_amplification,
                                            env1_minimal_peak_response,
                                            env1_benefit1,
                                            env1_benefit2,
                                            env1_max_duration_bias,
                                            env1_t_signal_on, 
                                            env1_t_signal_off,
                                            phenotypeA[thread_ID*N_replicates_per_thread+i],                                           
                                            fitnessA[thread_ID*N_replicates_per_thread+i],  
                                            RS[thread_ID],                                           
                                            &timepoint,                                           
                                            &end_state,
                                            thread_ID,
                                            env1_t_development,
                                            signal_profile);        
            }
            if(end_state==0)
                i--;        
            else
                fitnessA2[thread_ID*N_replicates_per_thread+i]=calc_replicate_fitness(&state_clone,1,env1_t_development,env1_t_signal_off,env1_response_amplification,env1_minimal_peak_response,genotype_clone.n_output_proteins);
            free_fixedevent(&state_clone);
#if PEAK_SEARCH
            free(state_clone.sampled_response);
#endif
        }
        for(i=0;i<N_replicates_per_thread;i++)        
        {              
            effect_of_effector=env2_initial_effect_of_effector; 
            end_state=1;
            signal_profile=NULL;         
            initialize_cell(&genotype_clone,
                            &state_clone, 
                            genotype_clone.ngenes, 
                            genotype_clone.nproteins,
                            genotype_clone.protein_pool,
                            genotype_clone.mRNA_decay_rate, 
                            mRNA, 
                            protein, 
                            RS[thread_ID],
                            env2_t_development,
                            env2_t_signal_off);
            set_signal( &state_clone,
                        env2_t_signal_on,   
                        env2_t_signal_off,
                        signal_profile,
                        env2_t_development,
                        env2_signal1_strength);    
            state_clone.t = 0.0;
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, env2_t_development,INITIALIZATION,thread_ID);        
            timepoint=0;
            while(state_clone.t<env2_t_development && end_state==1)
            {
                do_single_timestep_plotting(&genotype_clone, 
                                            &state_clone, 
                                            &rate_clone,                                                                                                                                       
                                            env2_signal2_strength,
                                            env2_response_amplification,
                                            env2_minimal_peak_response,
                                            env2_benefit1,
                                            env2_benefit2,
                                            env2_max_duration_bias,
                                            env2_t_signal_on,
                                            env2_t_signal_off,
                                            phenotypeB[thread_ID*N_replicates_per_thread+i],
                                            fitnessB[thread_ID*N_replicates_per_thread+i],                                            
                                            RS[thread_ID],                                          
                                            &timepoint,                                            
                                            &end_state,
                                            thread_ID,
                                            env2_t_development,
                                            signal_profile);
            }
            if(end_state==0)
                i--;
            else
                fitnessB2[thread_ID*N_replicates_per_thread+i]=calc_replicate_fitness(&state_clone,2,env2_t_development,env2_t_signal_off,env2_response_amplification,env2_minimal_peak_response,genotype_clone.n_output_proteins);            
            free_fixedevent(&state_clone);
#if PEAK_SEARCH
            free(state_clone.sampled_response);
#endif
        }        
    }//end of parallel section   
    /*Output protein numbers at each timepoint*/
    /*Each protein has its own file: A is env1 and B env2*/
    /*In each file, columns are replicates, rows are timepoints*/
    for(l=0;l<genotype->nproteins;l++)
    {
        snprintf(filename1,sizeof(char)*32,"phenotypeA_%i",l);
        snprintf(filename2,sizeof(char)*32,"phenotypeB_%i",l);
        fp1=fopen(filename1,"w");
        fp2=fopen(filename2,"w");
        for(m=0;m<N_TIMEPOINTS;m++)
        {
            for(n=0;n<N_REPLICATES;n++)
            {
                fprintf(fp1,"%f ",phenotypeA[n][l][m]);
                fprintf(fp2,"%f ",phenotypeB[n][l][m]);
            }
            fprintf(fp1,"\n");
            fprintf(fp2,"\n");
        }
        fclose(fp1);
        fclose(fp2);
    }
    /*Output growth rate*/
    /*Each environment has its file: A is env1 and B env2*/
    /*In each file, columns are replicates, rows are timepoints*/
    fp1=fopen("fitnessA","w");
    fp2=fopen("fitnessB","w");    
//    for(m=0;m<N_TIMEPOINTS;m++)
//    {
//        for(n=0;n<N_REPLICATES;n++)
//        {
//            fprintf(fp1,"%f ",fitnessA[n][m]);
//            fprintf(fp2,"%f ",fitnessB[n][m]);
//        }
//        fprintf(fp1,"\n");
//        fprintf(fp2,"\n");
//    }  
    for(l=0;l<N_REPLICATES;l++)
    {
        fprintf(fp1,"%f\n",fitnessA2[l]);
        fprintf(fp2,"%f\n",fitnessB2[l]);
    }
    fclose(fp1);
    fclose(fp2);
    
    fp1=fopen("output_protein_id.txt","w");
    for(l=0;l<genotype->n_output_proteins;l++)
        fprintf(fp1,"%d\n",genotype->output_protein_id[l]);
    fclose(fp1);
}
#endif /*End of plotting functions*/


/*
 ************ begin of mutation functions **************
 */
/*mutate an effector gene to normal tf gene*/
void mut_effector2TF(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene,protein_id;
    while(1)
    {
        which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
        if(genotype->protein_identity[genotype->which_protein[which_gene]][1]!=-1)
            break;
    }
    mut_record->which_gene=which_gene;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'t');
    genotype->n_output_genes--;
}

void reproduce_effector2TF(Genotype *genotype, Mutation *mut_record)
{
    int which_gene, protein_id;
    which_gene=mut_record->which_gene;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'t');
    genotype->n_output_genes--;
}

/*single nucleic acid substitution in cis-reg*/
void mut_substitution(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_nucleotide, which_gene,i,N_BS_bf_mutation,flag_difference;
    char *Genome,nucleotide; 
    AllTFBindingSites *container; // used to store binding sites before substitution
                                  // we compare whether substituion changes any binding sites,
                                  // in order to determine whether mutation creates a unique cis-reg     
                                  // sequence whose binding configuration always needs computation    
    /*points to the current cis-reg*/
    Genome= &genotype->cisreg_seq[0][0]; 
    /*whether to simulate bias in the frequency of different substitutions*/
    #if !SIMPLE_SUBSTITUTION 
        float random;
        while(1)
        {
            /*this is the nucleic acid that's going to be mutated*/ 
            which_nucleotide=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
            random=RngStream_RandU01(RS);
            if (Genome[which_nucleotide]=='a'||Genome[which_nucleotide]=='t')
            {
                 if(random<=0.2)break;  // 20% of mutations hit A or T
            }	            
            else
            { 
                 if(random>0.2)break; //80% hit C or G
            }
        }  
        /*this is gene hit by the mutation*/
        which_gene=which_nucleotide/CISREG_LEN;
        /*calculate and store the distribution of binding sites before mutation*/
        calc_all_binding_sites_copy(genotype,which_gene);
        container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
            container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
            container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
        }
        /*generate new nucleic acid*/
        random=RngStream_RandU01(RS);
        switch (Genome[which_nucleotide])
        {
            case 'a':
                if(random<=0.2)
                    Genome[which_nucleotide]='c'; // 20% of mutation to A changes it to C
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='t'; // 30% changes it to T
                    else
                        Genome[which_nucleotide]='g'; // 50% changes it to G
                }
                break;
            case 't':
                if(random<=0.2)
                    Genome[which_nucleotide]='g';
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='a';
                    else
                        Genome[which_nucleotide]='c';
                }
                break;        
            case 'c':
                if(random<=0.2)
                    Genome[which_nucleotide]='g';
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='t';
                    else
                        Genome[which_nucleotide]='a';
                }
                break;           
            case 'g':  
                if(random<=0.2)
                    Genome[which_nucleotide]='c';
                else
                {
                    if((random-0.2)<=0.3)
                        Genome[which_nucleotide]='a';
                    else
                        Genome[which_nucleotide]='t';
                }
                break;
        }   
    #else // if do not simulate bias in substitution
        /*this is the nucleic acid that's going to be mutated*/        
        which_nucleotide=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
        /*this is gene hit by the mutation*/
        which_gene=which_nucleotide/CISREG_LEN;
        /*calculate and store the distribution of binding sites before mutation*/
        calc_all_binding_sites_copy(genotype,which_gene);
        container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
            container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
            container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
        }
        /*generate new nucleic acid*/
        nucleotide=set_base_pair(RngStream_RandU01(RS));
        while(nucleotide==Genome[which_nucleotide]) // make sure we get a different nucleic acid
            nucleotide=set_base_pair(RngStream_RandU01(RS));
        Genome[which_nucleotide]=nucleotide;
    #endif    
    /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster*/
    calc_all_binding_sites_copy(genotype,which_gene);    
    if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])
    {
        update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1);  
    }
    else
    {
        /*assuming no change in binding sites*/
        flag_difference=0;
        /*comparing binding sites pair by pair*/
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            if(container[i].BS_pos!=genotype->all_binding_sites[which_gene][i].BS_pos ||
                container[i].tf_id!=genotype->all_binding_sites[which_gene][i].tf_id ||
                container[i].mis_match!=genotype->all_binding_sites[which_gene][i].mis_match)
            {    
                flag_difference=1;
                break;
            }
        }
        if(flag_difference==1)
           update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
    free(container);
    /*record mutation info*/
    mut_record->which_nucleotide=which_nucleotide;
    mut_record->which_gene=which_gene;
    mut_record->nuc_diff[0]=Genome[which_nucleotide];
    genotype->recalc_TFBS[which_gene]=1;  /*the binding sites on this cis-reg needs to be updated*/    
}

/*
 * Replay substitution. Used when replay the evolution for analysis
 */
void reproduce_substitution(Genotype *genotype, Mutation *mut_record)
{    
    char *Genome;
    int i, N_BS_bf_mutation,which_gene,flag_difference;    
    AllTFBindingSites *container;    
    /*get the mutated gene from record*/
    which_gene=mut_record->which_gene;
    /*calculate and store the distribution of binding sites before mutation*/
    calc_all_binding_sites_copy(genotype,which_gene); 
    container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
    N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
    for(i=0;i<N_BS_bf_mutation;i++)
    {
        container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
        container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
        container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
    }
    /*apply the mutation from record*/
    Genome= &genotype->cisreg_seq[0][0];    
    Genome[mut_record->which_nucleotide]=mut_record->nuc_diff[1]; 
    /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster*/
    calc_all_binding_sites_copy(genotype,which_gene);    
    if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])
    {
        update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1);  
    }
    else
    {
        flag_difference=0;
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            if(container[i].BS_pos!=genotype->all_binding_sites[which_gene][i].BS_pos ||
                container[i].tf_id!=genotype->all_binding_sites[which_gene][i].tf_id ||
                container[i].mis_match!=genotype->all_binding_sites[which_gene][i].mis_match)
            {    
                flag_difference=1;
                break;
            }
        }
        if(flag_difference==1)
           update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
    free(container);    
    genotype->recalc_TFBS[which_gene]=1;  /*recalc TFBS*/
}

/**
 *Deletion whole cis-reg sequence
 */
void mut_whole_gene_deletion(Genotype *genotype, Mutation *mut_record, RngStream RS) // any gene can be deleted
{
    int i,j;
    int which_gene, protein_id;
    char *temp;        
    /* check which genes can be deleted*/      
    if(genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF>1)//if there are more than one copies of genes of non-output genes  
    {
        if(genotype->n_output_genes==1)//if the effector have only one copies
        {            
            while(1) //the last copy of the effector cannot be deleted. choose a different gene
            {   
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                protein_id=genotype->which_protein[which_gene];
                if(genotype->protein_identity[protein_id][1]==NON_OUTPUT_PROTEIN)
                    break;
            }
        }
        else //otherwise any gene can be deleted
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    }
    else //if there's only one copy of non-output gene
    {        
        while(1)
        {
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
            protein_id=genotype->which_protein[which_gene];
            if(genotype->protein_identity[protein_id][1]!=NON_OUTPUT_PROTEIN)//obviously we can only delete a copy of the effector gene
                break;                           
        }
    }
    /*record mutation info*/
    mut_record->which_gene=which_gene;     
    /*points the first nucleotide of the cis-reg sequence to be deleted*/
    temp = &genotype->cisreg_seq[which_gene][0];   
    /* shift the cisreg array to overwrite the cis_seq to be deleted */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-1-which_gene);i++) 
    {				
        *temp=*(temp+CISREG_LEN);         
        temp++;				
    }    
    /* if the to-be-deleted gene is a tf gene*/
    protein_id=genotype->which_protein[which_gene];      
    /*if this tf has only one copy of gene, then we'll delete the binding seq of this tf and remove the tf from locus_specific_TF_behavior*/ 
    if(genotype->protein_pool[protein_id][0][0]==1) 
    {            
        /* shift the tf_reg array to overwrite the binding sequence to be deleted */
        temp=&genotype->tf_binding_seq[protein_id][0];     
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-1-protein_id);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        }
        /* shift the tf_reg_rc array to overwrite the binding sequence to be deleted */    
        temp=&genotype->tf_binding_seq_rc[protein_id][0];
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-1-protein_id);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        } 
        /*UPDATE locus_specific_TF_behavior*/
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
        {
            for(j=N_SIGNAL_TF;j<genotype->nproteins;j++)
                genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i][j+1];
        }    
    }     
    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
    for(i=which_gene;i<genotype->ngenes-1;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];            
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->translation_rate[i]=genotype->translation_rate[i+1];  
        for(j=0;j<NPROTEINS;j++)
            genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i+1][j];
        genotype->recalc_TFBS[i]=1;        
    }
    /* if the to-be-deleted gene is an effector gene, change n_output_genes as well*/
    if(genotype->protein_identity[protein_id][1]!=NON_OUTPUT_PROTEIN)
        genotype->n_output_genes--;    
    /* now change protein_pool and cisreg_cluster*/   
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);  
    genotype->ngenes--;
}

void reproduce_whole_gene_deletion(Genotype *genotype, Mutation *mut_record) // any gene can be deleted
{    
    int which_gene, i,j;
    char *temp;    
    int protein_id;       
    which_gene=mut_record->which_gene;     
    temp = &genotype->cisreg_seq[which_gene][0];	
    for(i=0;i<CISREG_LEN*(genotype->ngenes-which_gene-1);i++) 
    {				
        *temp=*(temp+CISREG_LEN);         
        temp++;				
    }
    protein_id=genotype->which_protein[which_gene];
    if(genotype->protein_pool[protein_id][0][0]==1) 
    {  
        temp=&genotype->tf_binding_seq[protein_id][0];       
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        }            
        temp=&genotype->tf_binding_seq_rc[protein_id][0];
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        }   
    }
    for(i=which_gene;i<genotype->ngenes-1;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];            
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->translation_rate[i]=genotype->translation_rate[i+1];   
        for(j=0;j<NPROTEINS;j++)
            genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i+1][j];
        genotype->recalc_TFBS[i]=1;
    }
    if(genotype->protein_identity[protein_id][1]!=NON_OUTPUT_PROTEIN)
        genotype->n_output_genes--; 
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);
    genotype->ngenes--;   
}

/*
 * Duplicate a whole cis-reg sequence
 */
void mut_duplication(Genotype *genotype, Mutation *mut_record, RngStream RS) 
{
    int which_gene, i, protein_id;
    char *temp1, *temp2;   
    
    if(genotype->n_output_genes>=MAX_OUTPUT_GENES) // too many effector gene copies
    {     
        //note that it's not possible to have too many effector gene copies and too many tf gene copies at the same time
        //because that will make duplication rate 0.
        while(1)
        {
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
            protein_id=genotype->which_protein[which_gene];
            if(genotype->protein_identity[protein_id][1]==-1)//if which_gene is non-output
                break;
        }        
    }
    else
    {
        if(genotype->ngenes-genotype->n_output_genes>=MAX_NON_OUTPUT_GENES) // too many non-output gene copies
        {
            while(1)
            {
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                protein_id=genotype->which_protein[which_gene];
                if(genotype->protein_identity[protein_id][1]!=-1)
                    break;
            }    
        }            
        else //any gene can be duplicated
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1); 
    }   
    /*record mutation info*/
    mut_record->which_gene=which_gene;  
    /* copy the promoter*/
    temp1=&genotype->cisreg_seq[which_gene][0]; /* points to the gene to be duplicated*/
    temp2=&genotype->cisreg_seq[genotype->ngenes][0]; /* points to the end of the effector gene */
    /* shift the sequences of the effector gene CISREG_LEN bp */
    for(i=0;i<CISREG_LEN;i++)
        *temp2++=*temp1++;           
    /* copy and paste info to the slot*/
    genotype->min_N_activator_to_transc[genotype->ngenes]=genotype->min_N_activator_to_transc[which_gene];
    genotype->active_to_intermediate_rate[genotype->ngenes]=genotype->active_to_intermediate_rate[which_gene];
    genotype->mRNA_decay_rate[genotype->ngenes]=genotype->mRNA_decay_rate[which_gene];
    genotype->protein_decay_rate[genotype->ngenes]=genotype->protein_decay_rate[which_gene];
    genotype->translation_rate[genotype->ngenes]=genotype->translation_rate[which_gene]; 
    for(i=0;i<NPROTEINS;i++)
        genotype->locus_specific_TF_behavior[genotype->ngenes][i]=genotype->locus_specific_TF_behavior[which_gene][i];
    genotype->recalc_TFBS[genotype->ngenes]=1;
    /* update protein_pool*/
    protein_id=genotype->which_protein[which_gene];    
    update_protein_pool(genotype,protein_id,which_gene,'d'); 
    /* update cisreg_cluster*/    
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    /* update gene numbers*/  
    if(genotype->protein_identity[protein_id][1]!=-1)//note duplication do not change nproteins           
        genotype->n_output_genes++;     
    genotype->ngenes++;
}

void reproduce_gene_duplication(Genotype *genotype, Mutation *mut_record) //any gene can be duplicated
{   
    int which_gene, i, protein_id;
    char *temp1, *temp2; 
    /*get the gene to be duplicated from record*/
    which_gene=mut_record->which_gene;      
    temp1=&genotype->cisreg_seq[which_gene][0]; 
    temp2=&genotype->cisreg_seq[genotype->ngenes][0];  
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;    
    genotype->min_N_activator_to_transc[genotype->ngenes]=genotype->min_N_activator_to_transc[which_gene];
    genotype->active_to_intermediate_rate[genotype->ngenes]=genotype->active_to_intermediate_rate[which_gene];
    genotype->mRNA_decay_rate[genotype->ngenes]=genotype->mRNA_decay_rate[which_gene];
    genotype->protein_decay_rate[genotype->ngenes]=genotype->protein_decay_rate[which_gene];
    genotype->translation_rate[genotype->ngenes]=genotype->translation_rate[which_gene]; 
    for(i=0;i<NPROTEINS;i++)
        genotype->locus_specific_TF_behavior[genotype->ngenes][i]=genotype->locus_specific_TF_behavior[which_gene][i];
    genotype->recalc_TFBS[genotype->ngenes]=1;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'d');     
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    if(genotype->protein_identity[protein_id][1]!=-1)
        genotype->n_output_genes++;       
    genotype->ngenes++;     
}

/*
 * Mutation to the binding sequence of a TF gene
 */
void mut_binding_sequence(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene, which_nucleotide, protein_id, i;
    char nucleotide;      
    char *tf_binding_seq, *tf_binding_seq_rc,*temp1,*temp2;
    /*get which gene to mutate*/
    which_gene=RngStream_RandInt(RS,0,genotype->ngenes-1);  
    /*if this TF has more than one copies of gene, then the mutation adds a new tf
     *which requires a new slot in tf_binding_seq and tf_binding_seq_rc to store the new binding seq*/
    protein_id=genotype->which_protein[which_gene];    
    if(genotype->protein_pool[protein_id][0][0]>1)
    {
        /*points to the first nucleotide of the binding sequence*/
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        /*points to an empty slot*/
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        /*copy the binding sequences to empty slots*/
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }
        /*point tf_binding_seq and tf_binding_seq_rc to the new slot so that we can apply mutation later*/
        tf_binding_seq=&genotype->tf_binding_seq[genotype->nproteins][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        /*update protein pool*/
        update_protein_pool(genotype,protein_id,which_gene,'c');          
        /* Update locus_specific_TF_behavior: assuming mutation to binding sequence changes nothing to the locus-specific behavior */   
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else /*if this tf has only one copy of gene, no new slot is required*/
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
    }
    /*mutate the binding sequence*/
    /*mutation only changes one nucleotide in the binding sequence*/
    which_nucleotide=RngStream_RandInt(RS,0,TF_ELEMENT_LEN-1);
    nucleotide=set_base_pair(RngStream_RandU01(RS));        
    while (nucleotide == tf_binding_seq[which_nucleotide])
        nucleotide=set_base_pair(RngStream_RandU01(RS));
    tf_binding_seq[which_nucleotide]=nucleotide;    
    /*record mutation info*/
    mut_record->which_gene=which_gene;
    mut_record->which_nucleotide=which_nucleotide;
    mut_record->nuc_diff[0]=nucleotide;    
    /* update the reverse complement sequence*/
    switch (nucleotide)
    {
        case 'g':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
    }  
    /* The binding sites on every promoter needs recalculation */
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=1;
    /*decide whether to update cisreg clusters. Mutation to binding seq may differentiate bs distributions among genes in a cluster*/       
    int new_clusters[NGENES][NGENES];  
    int genes_in_cluster[NGENES];
    int N_genes_in_cluster,no_difference,reference_gene,gene_to_be_sorted;
    int N_new_clusters,N_genes_in_new_cluster,j,k;
    calc_all_binding_sites(genotype); 
    i=N_SIGNAL_TF;
    while(genotype->cisreg_cluster[i][0]!=-1)/*check each cisreg cluster*/
    {       
        for(j=0;j<NGENES;j++)
        {
            for(k=0;k<NGENES;k++)
                new_clusters[j][k]=-1;
        }  
        /*copy an original cluster and count genes in the cluster*/
        N_genes_in_cluster=0;
        while(genotype->cisreg_cluster[i][N_genes_in_cluster]!=-1)
        {
            genes_in_cluster[N_genes_in_cluster]=genotype->cisreg_cluster[i][N_genes_in_cluster];
            N_genes_in_cluster++;
        }
        N_new_clusters=0; 
        /*the while loop below sort genes in a cluster into groups based on whether they have the same BS distributions*/ 
        /*We use one gene in the original cluster as a reference, and sort all genes 
         *that have the same binding sites as the reference gene into a new cluster. 
         *Genes that are different from the reference gene are sorted again similarly 
         *through iterations. After all genes are sorted, we check if the original 
         *cluster turns into multiple new cluster.*/
        while(N_genes_in_cluster>0)//run until every gene in the original cluster has been sorted into new clusters
        {
            reference_gene=genes_in_cluster[0];//compare each gene to the first gene in cluster            
            N_genes_in_new_cluster=0;
            which_gene=0;//start comparison from the first gene in the cluster
            while(which_gene<N_genes_in_cluster)
            {
                no_difference=1;
                gene_to_be_sorted=genes_in_cluster[which_gene];                
                if(genotype->binding_sites_num[gene_to_be_sorted]==genotype->binding_sites_num[reference_gene])
                {
                    for(j=0;j<genotype->binding_sites_num[reference_gene];j++)
                    {
                        if(genotype->all_binding_sites[reference_gene][j].BS_pos!=genotype->all_binding_sites[gene_to_be_sorted][j].BS_pos ||
                            genotype->all_binding_sites[reference_gene][j].tf_id!=genotype->all_binding_sites[gene_to_be_sorted][j].tf_id ||
                            genotype->all_binding_sites[reference_gene][j].mis_match!=genotype->all_binding_sites[gene_to_be_sorted][j].mis_match)
                        {
                            no_difference=0;
                            break;
                        }
                    }
                }
                else
                    no_difference=0;                
                if(no_difference)//if the gene has the same binding sites as the reference gene
                {
                    /*put the gene into the new cluster*/
                    new_clusters[N_new_clusters][N_genes_in_new_cluster]=genes_in_cluster[which_gene];
                    N_genes_in_new_cluster++;
                    /*shift to remove the gene from the copy of the original cluster*/
                    for(j=which_gene;j<N_genes_in_cluster-1;j++)
                        genes_in_cluster[j]=genes_in_cluster[j+1];
                    N_genes_in_cluster--;                
                }
                else //otherwise, compare next gene with the reference gene                
                    which_gene++;                
            }
            N_new_clusters++;
        }
        if(N_new_clusters!=1)//if the original cluster turns into multiple clusters
            update_cisreg_cluster(genotype,-1,'c',new_clusters,N_new_clusters,i);
        i++;
    }
    /* Calling calc_all_binding_sites reset recalc_TFBS to 0, so we need to turn them on again */
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=1;  
}

void reproduce_mut_binding_sequence(Genotype *genotype, Mutation *mut_record)
{    
    int which_gene, which_nucleotide, protein_id, i;    
    char *tf_binding_seq, *tf_binding_seq_rc, *temp1,*temp2;     
    which_gene=mut_record->which_gene;
    protein_id=genotype->which_protein[which_gene];  
    if(genotype->protein_pool[protein_id][0][0]>1)
    {        
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }      
        tf_binding_seq=&genotype->tf_binding_seq[genotype->nproteins][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[genotype->nproteins][0];        
        update_protein_pool(genotype,protein_id,which_gene,'c');  
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else 
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
    }        			
    which_nucleotide=mut_record->which_nucleotide;        
    tf_binding_seq[which_nucleotide]=mut_record->nuc_diff[1];
    switch (tf_binding_seq[which_nucleotide])
    {
        case 'g':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
    }     
    for(i=0;i<genotype->ngenes;i++)
        genotype->recalc_TFBS[i]=1;
    calc_all_binding_sites(genotype);
    int new_clusters[NGENES][NGENES],genes_in_cluster[NGENES];
    int N_genes_in_cluster,no_difference,reference_gene,gene_to_be_sorted;
    int N_new_clusters,N_genes_in_new_cluster,j,k;
    i=N_SIGNAL_TF;
    while(genotype->cisreg_cluster[i][0]!=-1)
    {        
        N_new_clusters=0;
        N_genes_in_cluster=0;
        for(j=0;j<NGENES;j++)
        {
            for(k=0;k<NGENES;k++)
                new_clusters[j][k]=-1;
        }  
        while(genotype->cisreg_cluster[i][N_genes_in_cluster]!=-1)
        {            
            genes_in_cluster[N_genes_in_cluster]=genotype->cisreg_cluster[i][N_genes_in_cluster];
            N_genes_in_cluster++;
        }
        while(N_genes_in_cluster>0)
        {
            reference_gene=genes_in_cluster[0];
            N_genes_in_new_cluster=0;
            which_gene=0;
            while(which_gene<N_genes_in_cluster)
            {
                no_difference=1;
                gene_to_be_sorted=genes_in_cluster[which_gene];                
                if(genotype->binding_sites_num[gene_to_be_sorted]==genotype->binding_sites_num[reference_gene])
                {
                    for(j=0;j<genotype->binding_sites_num[reference_gene];j++)
                    {
                        if(genotype->all_binding_sites[reference_gene][j].BS_pos!=genotype->all_binding_sites[gene_to_be_sorted][j].BS_pos ||
                            genotype->all_binding_sites[reference_gene][j].tf_id!=genotype->all_binding_sites[gene_to_be_sorted][j].tf_id ||
                            genotype->all_binding_sites[reference_gene][j].mis_match!=genotype->all_binding_sites[gene_to_be_sorted][j].mis_match)
                        {
                            no_difference=0;
                            break;
                        }
                    }
                }
                else
                    no_difference=0;                
                if(no_difference)
                {
                    new_clusters[N_new_clusters][N_genes_in_new_cluster]=genes_in_cluster[which_gene];
                    N_genes_in_new_cluster++;
                    for(j=which_gene;j<N_genes_in_cluster-1;j++)
                        genes_in_cluster[j]=genes_in_cluster[j+1];
                    N_genes_in_cluster--;                
                }
                else 
                    which_gene++;
                
            }
            N_new_clusters++;
        }
        if(N_new_clusters!=1)
            update_cisreg_cluster(genotype,-1,'c',new_clusters,N_new_clusters,i);
        i++;
    }    
    for(i=0;i<genotype->ngenes;i++)
        genotype->recalc_TFBS[i]=1;   
}

/* Mutations to the rate of mRNA_decay, translation, protein_decay, and pic_disassembly 
 * will be mutated. 
 */
void mut_kinetic_constant(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    float random;    
    int which_gene;
    float total_mut_flux=   proportion_mut_kdis+
                            proportion_mut_mRNA_decay+
                            proportion_mut_protein_decay+
                            proportion_mut_translation_rate+
                            proportion_mut_cooperation;
    
    random=RngStream_RandU01(RS);
    which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);    
    /*record mutation info*/
    mut_record->which_gene=which_gene;
    /********************************** mut kdis *******************************/
    if(random<=proportion_mut_kdis/total_mut_flux) 
    {  
        genotype->active_to_intermediate_rate[which_gene]=mut_make_new_value(genotype->active_to_intermediate_rate[which_gene],
                                                                            miu_ACT_TO_INT_RATE,
                                                                            sigma_ACT_TO_INT_RATE,
                                                                            MAX_ACT_TO_INT_RATE,
                                                                            MIN_ACT_TO_INT_RATE,
                                                                            RS,
                                                                            mut_record);
        /*record mutation info*/
        mut_record->kinetic_type=0;
        mut_record->kinetic_diff=genotype->active_to_intermediate_rate[which_gene];
    }
    else /******************************** mut mRNAdecay **********************/
    {
        random-=proportion_mut_kdis/total_mut_flux;
        if(random<=proportion_mut_mRNA_decay/total_mut_flux) 
        {       
            genotype->mRNA_decay_rate[which_gene]=mut_make_new_value(genotype->mRNA_decay_rate[which_gene],miu_mRNA_decay,sigma_mRNA_decay, MAX_MRNA_DECAY, MIN_MRNA_DECAY,RS,mut_record);
            /*record mutation info*/
            mut_record->kinetic_type=1;
            mut_record->kinetic_diff=genotype->mRNA_decay_rate[which_gene];
        }
        else /*************************** mut translation *********************/
        {
            random-=proportion_mut_mRNA_decay/total_mut_flux;
            if(random<=proportion_mut_translation_rate/total_mut_flux) 
            {       
                genotype->translation_rate[which_gene]=mut_make_new_value(genotype->translation_rate[which_gene],miu_translation_init,sigma_translation_init, MAX_TRANSLATION_RATE, MIN_TRANSLATION_RATE,RS,mut_record);        
                mut_record->kinetic_type=2;
                mut_record->kinetic_diff=genotype->translation_rate[which_gene];
            }
            else /********************* mut protein decay **********************/
            {  
                random-=proportion_mut_translation_rate/total_mut_flux;
                if(random<=proportion_mut_protein_decay/total_mut_flux)
                {                   
                    genotype->protein_decay_rate[which_gene]=mut_make_new_value(genotype->protein_decay_rate[which_gene],miu_protein_decay,sigma_protein_decay,MAX_PROTEIN_DECAY,MIN_PROTEIN_DECAY,RS,mut_record);                  
                    mut_record->kinetic_type=3;
                    mut_record->kinetic_diff=genotype->protein_decay_rate[which_gene];
                }
                else /****************** mut cooperation ***********************/
                {
                    while(genotype->which_protein[which_gene]!=genotype->nproteins-1)
                        which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                    mut_record->which_gene=which_gene;
                    if(genotype->min_N_activator_to_transc[which_gene]==1)
                        genotype->min_N_activator_to_transc[which_gene]=2;
                    else
                        genotype->min_N_activator_to_transc[which_gene]=1;
                    mut_record->kinetic_type=4;
                    mut_record->kinetic_diff=(float)genotype->min_N_activator_to_transc[which_gene];
                    update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
                }
            } 
        }
    }
}

float mut_make_new_value(float old_val, float miu, float sigma, float upper_bound, float lower_bound, RngStream RS, Mutation *mut_record)
{
    float new_val;
    new_val=old_val;
    while(new_val==old_val) 
    {
        new_val=old_val*pow(10.0,sigma*gasdev(RS)+mutational_regression_rate*(miu-log10(old_val)));
        while(new_val>=upper_bound || new_val<=lower_bound)
        {
            new_val=old_val*pow(10.0,sigma*gasdev(RS)+mutational_regression_rate*(miu-log10(old_val)));
            mut_record->N_hit_bound++;
        }            
    }
    return new_val;
}


void reproduce_mut_kinetic_constant(Genotype *genotype, Mutation *mut_record)
{    
    int which_gene;            
    which_gene=mut_record->which_gene; /* which gene */ 
    switch (mut_record->kinetic_type)
    {
        case 0: /* mut kdis */
            genotype->active_to_intermediate_rate[which_gene]=mut_record->kinetic_diff;
            break;
        case 1: /* mut mRNAdecay */ 
            genotype->mRNA_decay_rate[which_gene]=mut_record->kinetic_diff;
            break;
        case 2: /* mut translation */                      
            genotype->translation_rate[which_gene]=mut_record->kinetic_diff;           
            break;        
        case 3: /* mut protein decay */
            genotype->protein_decay_rate[which_gene]=mut_record->kinetic_diff;            
            break; 
        case 4:
            genotype->min_N_activator_to_transc[which_gene]=(int)mut_record->kinetic_diff;
            update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
            break;
    }
}

void mut_locus_specific_tf_behavior(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int i, gene_id, protein_id,which_cluster;
    /*which locus to mutate*/
    gene_id=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    /*which TF is flipped*/
    protein_id=RngStream_RandInt(RS,0,genotype->nproteins-1);
    /*flip the chosen TF*/
    genotype->locus_specific_TF_behavior[gene_id][protein_id]=(genotype->locus_specific_TF_behavior[gene_id][protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR;
    /*if gene_id does not have a unique cis-reg, this mutation separate gene_id from its original cisreg cluster*/
    which_cluster=genotype->which_cluster[gene_id];
    if(genotype->cisreg_cluster[i][1]!=NA)
        update_cisreg_cluster(genotype,gene_id,'l',NULL,-1,-1);
    
    mut_record->which_gene=gene_id;
    mut_record->which_protein=protein_id;   
    genotype->recalc_TFBS[gene_id]=1;    
}

void mut_identity(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int tf_id,protein_id,i;   
    char *tf_binding_seq,*tf_binding_seq_rc,*temp1,*temp2; 
    /*which tf gene to mutate*/
    tf_id = RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1); // the first sensor tf must be an activator, therefore is not subject to mutation
    protein_id=genotype->which_protein[tf_id];  
    /*save record*/
    mut_record->which_gene=tf_id;       
    /* if this tf gene has more than one copies, the mutation adds a new protein*/
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        /*give the new protein its own binding sequence*/
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0]; 
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        } 
        update_protein_pool(genotype,protein_id,tf_id,'e');
        /* update_protein_pool put the new protein at nproteins and then increases nproteins by 1, 
        * so the new protein is at nproteins-1 now. Note that N_act and N_rep is updated in update_protein_pool*/ 
        genotype->protein_identity[genotype->nproteins-1][0]=(genotype->protein_identity[genotype->nproteins-1][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR;    
        /*the locus specific behavior is not flipped by this mutation*/
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->protein_identity[genotype->nproteins-1][0];
    }
    else //otherwise we just flip the property of an exisiting TF
    {
        genotype->protein_identity[protein_id][0]=(genotype->protein_identity[protein_id][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id][0]==ACTIVATOR)//mutate to activator
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else //otherwise to repressor
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][protein_id]=genotype->protein_identity[protein_id][0];
    }    
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=1; /* recalculate binding sites on every promoter */        
}

void reproduce_mut_identity(Genotype *genotype, Mutation *mut_record)
{
    int tf_id, protein_id,i;
    char *tf_binding_seq,*tf_binding_seq_rc,*temp1,*temp2;    
    tf_id = mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }     
        update_protein_pool(genotype,protein_id,tf_id,'e');  
        genotype->protein_identity[genotype->nproteins-1][0]=(genotype->protein_identity[genotype->nproteins-1][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->protein_identity[genotype->nproteins-1][0];
    }
    else
    {
        genotype->protein_identity[protein_id][0]=(genotype->protein_identity[protein_id][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id][0]==ACTIVATOR)
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][protein_id]=genotype->protein_identity[protein_id][0];
    }    
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=1;   
}

/*
 * mutate affinity of TF
 */
void mut_Kd(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int i;   
    float new_Kd; 
    int tf_id,protein_id;
    char *tf_binding_seq,*tf_binding_seq_rc,*temp1,*temp2;
    /*which TF to mutate*/
    tf_id=RngStream_RandInt(RS,0,genotype->ngenes-1);
    protein_id=genotype->which_protein[tf_id];   
    /*generate a new koff */    
    new_Kd=mut_make_new_value(genotype->Kd[protein_id],miu_Kd,sigma_Kd,MAX_KD,MIN_KD,RS,mut_record);  
    /* if this tf gene has more than one copies, the mutation adds a new protein*/    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }
        update_protein_pool(genotype,protein_id,tf_id,'f');
        /* update_protein_pool put the new protein at nproteins and then increases nproteins by 1, 
         * so the new protein is at nproteins-1 now.*/
        genotype->Kd[genotype->nproteins-1]=new_Kd;  
        /* Update locus_specific_TF_behavior: assuming mutation to Kd changes nothing to the locus-specific behavior */   
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }                                                                                                           
    else
        genotype->Kd[protein_id]=new_Kd;     
    /*record mutation*/
    mut_record->which_gene=tf_id;
    mut_record->kinetic_diff=new_Kd; 
    /* recalculate binding sites on every promoter */
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;
}

void reproduce_mut_Kd(Genotype *genotype, Mutation *mut_record)
{   
    int tf_id,protein_id,i;
    char *tf_binding_seq, *tf_binding_seq_rc, *temp1, *temp2;    
    tf_id=mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }    
        update_protein_pool(genotype,protein_id,tf_id,'f');  
        genotype->Kd[genotype->nproteins-1]=mut_record->kinetic_diff;  
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else
        genotype->Kd[protein_id]=mut_record->kinetic_diff;
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;   
}

/*
 *Wrapper of all mutation functions
 */
void mutate(Genotype *genotype, RngStream RS, Mutation *mut_record)
{
    int i;
#if !SET_BS_MANUALLY    
    draw_mutation(genotype, &(mut_record->mut_type),RS);
#else   
    float random;
    random=RngStream_RandU01(RS);
    if(random<proportion_mut_koff)
        mut_record->mut_type='f';
    else
        mut_record->mut_type='k';    
#endif
    /*reset records of nucleotides*/
    for(i=0;i<3;i++)
        mut_record->nuc_diff[i]='\0';
    switch (mut_record->mut_type)
    {
        case 's': //substitution in cis-reg       		
            mut_substitution(genotype,mut_record,RS);
            break;        
        case 'w': // whole gene deletion.          
            mut_whole_gene_deletion(genotype,mut_record,RS);            
            break;
        case 'd': // Whole gene duplication                  
            mut_duplication(genotype,mut_record,RS);            
            break;
        case 'c': //binding sequence 
            mut_binding_sequence(genotype,mut_record,RS);            
            break;
        case 'k': //mutations to kinetic constants        
            mut_kinetic_constant(genotype, mut_record,RS);           
            break;
        case 'e': //activator to repressor or the reverse
            mut_identity(genotype, mut_record, RS);            
            break;
        case 'f': //mutations to the koff of a tf
            mut_Kd(genotype,mut_record,RS);            
            break;
        case 't': //effector to regular TF
            mut_effector2TF(genotype,mut_record,RS);
            break;
    }
}

/* this function perform mutation indicated by input. Used to reproduce the genotype following a serial of mutations*/
void reproduce_mutate(Genotype *genotype, Mutation *mut_record,RngStream RS)
{ 
    switch (mut_record->mut_type)
    {
        case 's': //substitution        		
            reproduce_substitution(genotype,mut_record);
            break;       
        case 'w': // whole gene deletion.         
            reproduce_whole_gene_deletion(genotype,mut_record);
            break;
        case 'd': // Whole gene duplication                 
            reproduce_gene_duplication(genotype,mut_record);
            break;        
        case 'c': //binding sequence        
            reproduce_mut_binding_sequence(genotype,mut_record);
            break;
        case 'k': //mutations in kinetic constants        
            reproduce_mut_kinetic_constant(genotype, mut_record);            
            break;
        case 'e': //changing the identity of a TF
            reproduce_mut_identity(genotype, mut_record);
            break;
        case 'f':
            reproduce_mut_Kd(genotype,mut_record);
            break; 
        case 't': //effector to regular TF
            reproduce_effector2TF(genotype,mut_record);
            break;
    }   
}

/* this function calculates the probability of different mutations based on
 * the current genotype, and draws a mutation accordingly.
 */
void draw_mutation(Genotype *genotype, char *mut_type, RngStream RS)
{
    float random;
    float tot_mut_rate=0.0;
    float tot_subs_rate, tot_dup_rate, tot_sil_rate, tot_mut_kin_rate, tot_mut_identity_rate, tot_mut_binding_seq_rate, tot_mut_koff_rate, tot_effector2TF; 
    int N_target_genes; 
    
    /* duplication rate*/    
    tot_dup_rate=0.0;
    N_target_genes=genotype->ngenes-N_SIGNAL_TF;//NA to sensor TF 
    if(genotype->ngenes-genotype->n_output_genes>=MAX_NON_OUTPUT_GENES)//too many non-output genes
        N_target_genes-=genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF; //do not duplicate non-output gene anymore
    if(genotype->n_output_genes>=MAX_OUTPUT_GENES)//too many effector gene
        N_target_genes-=genotype->n_output_genes;//do not duplicate effector gene anymore
    tot_dup_rate=N_target_genes*DUPLICATION;
    tot_mut_rate+=tot_dup_rate; 
    
    /* silencing rate*/ 
    N_target_genes=0;
    if(genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF>1)//if there's more than one copy of non-output genes
        N_target_genes+=genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF;
    if(genotype->n_output_genes>1)
        N_target_genes+=genotype->n_output_genes; 
    tot_sil_rate=N_target_genes*SILENCING; 
    tot_mut_rate+=tot_sil_rate;    
    
    /* calc total susbtitution rate*/
    tot_subs_rate=(genotype->ngenes-N_SIGNAL_TF)*CISREG_LEN*SUBSTITUTION; //NA to the sensor TF gene
    tot_mut_rate+=tot_subs_rate;   

    /* mut in kinetic constants */    
    tot_mut_kin_rate=(genotype->ngenes-N_SIGNAL_TF)*MUTKINETIC*(proportion_mut_kdis+proportion_mut_mRNA_decay+proportion_mut_protein_decay+proportion_mut_translation_rate); // NA to the sensor TF gene
    tot_mut_kin_rate+=(genotype->n_output_genes)*MUTKINETIC*proportion_mut_cooperation;
    tot_mut_rate+=tot_mut_kin_rate;
    
    /* mut in binding seq*/
    tot_mut_binding_seq_rate=genotype->ngenes*MUTKINETIC*proportion_mut_binding_seq; // NA to the effector genes
    tot_mut_rate+=tot_mut_binding_seq_rate;    
    
    /* mut in identity*/
    tot_mut_identity_rate=(genotype->ngenes-N_SIGNAL_TF)*MUTKINETIC*proportion_mut_identity; // NA to the sensor TF gene
    tot_mut_rate+=tot_mut_identity_rate;
    
    /* mut in koff*/
    tot_mut_koff_rate=genotype->ngenes*MUTKINETIC*proportion_mut_koff;  
    tot_mut_rate+=tot_mut_koff_rate;  
    
    /*effector to regular TF*/
    if(genotype->n_output_genes>1)
        tot_effector2TF=genotype->n_output_genes*MUTKINETIC*proportion_effector2TF;
    else
        tot_effector2TF=0.0;
    tot_mut_rate+=tot_effector2TF;
    
    /*Draw a mutation based on the above rates*/
    random=RngStream_RandU01(RS);    
    if(random<=tot_subs_rate/tot_mut_rate)    
        *mut_type='s';                      /* substitution*/
    else
    {   
        random-=tot_subs_rate/tot_mut_rate;
        if(random<= tot_dup_rate/tot_mut_rate)
            *mut_type='d';                          /* gene duplication */
        else
        {
            random-=tot_dup_rate/tot_mut_rate;
            if(random<=tot_sil_rate/tot_mut_rate)
                *mut_type='w';                              /* gene deletion*/  
            else
            {
                random-=tot_sil_rate/tot_mut_rate;                
                if(random<=tot_mut_kin_rate/tot_mut_rate)
                    *mut_type='k';                  /* mut kinetic const*/
                else
                {
                    random-=tot_mut_kin_rate/tot_mut_rate;
                    if(random<=tot_mut_binding_seq_rate/tot_mut_rate)
                        *mut_type='c';              /* mut binding seq*/
                    else
                    {
                        random-=tot_mut_binding_seq_rate/tot_mut_rate;
                        if(random<=tot_mut_koff_rate/tot_mut_rate)
                            *mut_type='f';          /* mut koff*/
                        else
                        {
                            random-=tot_mut_koff_rate/tot_mut_rate;
                            if(random<=tot_effector2TF/tot_mut_rate)
                                *mut_type='t';       /* effector to regular TF*/
                            else
                                *mut_type='e';           /* mut identity of a TF */
                        }    
                    }
                }               
            }
        }
    }  
}

/*
 *We use which_protein to look for the protein encoded by a given gene copy,
 *and protein_pool to look for genes encoding a given protein. These two tables
 *are updated upon mutations.
 */
void update_protein_pool(Genotype *genotype, int which_protein, int which_gene, char mut_type)
{
    int i, j, gene_id, which_output; 
    /*Protein_pool stores the numbers of gene copies that encode a given protein, and the ids of these gene copies.
     *One important thing is that the genes encoding a given protein are not stored by the order of their ids in protein_pool.
     *To delete a gene, which might be the only gene encoding a given protein, we shift protein_pool to overwrite the to-be-deleted gene
     *We need to update the ids of the remaining genes and proteins
     *For gene duplication, the new gene is always add to the end of the list of genes encoding a given protein.
     *A new protein is also add to the end of protein_pool
     *which_protein can be updated easily; changing which protein a gene encodes is always easy. For deletion, 
     *we just shift the array to overwrite the to-be-deleted gene and update the ids of the remaining genes.*/ 
    switch (mut_type)
    {
        case 'w':/*a whole gene deletion*/        
            if(genotype->protein_pool[which_protein][0][0]==1) /* if this is the only gene copy,we also need to delete a protein*/
            {   
                /*
                 * UPDATE protein_pool for protein>=which_protein
                 */
                /*shift protein>which_protein to overwrite the to-be-deleted protein*/
                for(i=which_protein;i<genotype->nproteins-1;i++)  
                {   
                    /*reset the portion of protein_pool to be overwritten*/
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                        genotype->protein_pool[i][1][j]=-1;
                    /*overwrite*/
                    genotype->protein_pool[i][0][0]=genotype->protein_pool[i+1][0][0]; //this is number of gene copies encoding a protein
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                    {
                        gene_id=genotype->protein_pool[i+1][1][j];//these are the gene copies encoding a protein
                        /*note that deletion changes the ids of the remaining genes!!! Any gene that is greater than which_gene is reduced by one*/
                        genotype->protein_pool[i][1][j]=(gene_id>which_gene)?gene_id-1:gene_id;
                    }   
                } 
                /*
                 * UPDATE which_protein
                 */
                /* update which_protein for gene<which_gene in which_protein*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i]<which_protein)?genotype->which_protein[i]:genotype->which_protein[i]-1;//the deletion also changes the ids of proteins
                /* shift and update which_protein for gene>=which_gene in which_protein*/                
                for(i=which_gene;i<genotype->ngenes-1;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i+1]>which_protein)?genotype->which_protein[i+1]-1:genotype->which_protein[i+1];  
                /*
                 * UPDATE the number of activators or that of repressors
                 */
                if(genotype->protein_identity[which_protein][0]==ACTIVATOR) 
                    genotype->N_act--;
                else
                    genotype->N_rep--;
                /* remove which_protein from protein_identity and Kd */              
                for(i=which_protein;i<genotype->nproteins-1;i++)
                {
                    genotype->protein_identity[i][0]=genotype->protein_identity[i+1][0];                    
                    genotype->Kd[i]=genotype->Kd[i+1];                    
                }                           
                /*
                 * UPDATE output_protein_id
                 */
                which_output=genotype->protein_identity[which_protein][1];                
                if(which_output!=-1)/*if which_protein is an effector protein*/
                {         
                    /*shift output_protein_id and update the id of proteins in output_protein_id*/
                    for(i=which_output;i<genotype->n_output_proteins-1;i++)
                        genotype->output_protein_id[i]=genotype->output_protein_id[i+1]-1; //output_protein_id is ordered ascendingly
                    /*update protein_identity similarly*/
                    for(i=which_protein;i<genotype->nproteins-1;i++)
                        genotype->protein_identity[i][1]=(genotype->protein_identity[i+1][1]==-1)?-1:genotype->protein_identity[i+1][1]-1;                    
                    genotype->n_output_proteins--;
                }
                else 
                {
                    /*just update the id of protein in output_protein_id*/
                    for(i=0;i<genotype->n_output_proteins;i++)
                        genotype->output_protein_id[i]=(genotype->output_protein_id[i]<which_protein)?genotype->output_protein_id[i]:genotype->output_protein_id[i]-1;
                    /*remove which_protein from protein_identity*/
                    for(i=which_protein;i<genotype->nproteins-1;i++)
                        genotype->protein_identity[i][1]=genotype->protein_identity[i+1][1];
                }  
                /*one less protein*/
                genotype->nproteins--;
                /* in the case, all genes need to recalc binding sites*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)                    
                    genotype->recalc_TFBS[i]=1; /* recalc BS */ 
            }  
            else /*if the protein has more than one genes*/
            {
                /*
                 * UPDATE protein_pool for protein>=which_protein
                 */
                /* find where is this which_gene*/
                i=0;
                while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++; 
                /*shift protein_pool to overwrite which_gene, and update ids of genes*/
                j=i;
                for(;i<genotype->protein_pool[which_protein][0][0]-1;i++)
                    genotype->protein_pool[which_protein][1][i]=(genotype->protein_pool[which_protein][1][i+1]>which_gene)?genotype->protein_pool[which_protein][1][i+1]-1:genotype->protein_pool[which_protein][1][i+1];// deletion changes the ids of genes!!!
                /*also update the ids for genes before j*/
                for(i=0;i<j;i++)
                    genotype->protein_pool[which_protein][1][i]=(genotype->protein_pool[which_protein][1][i]>which_gene)?genotype->protein_pool[which_protein][1][i]-1:genotype->protein_pool[which_protein][1][i];                            
                /*one less copy encoding which_protein*/
                genotype->protein_pool[which_protein][0][0]--;                
                /*update the ids of genes in protein>which_protein*/                
                for(i=which_protein+1;i<genotype->nproteins;i++)  
                {   
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                        genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]>which_gene)?genotype->protein_pool[i][1][j]-1:genotype->protein_pool[i][1][j];//note that deletion changes the ids of genes!!!
                }
                /*
                 * UPDATE which_protein
                 */
                /*shift which_protein to delete which_gene*/
                for(i=which_gene;i<genotype->ngenes-1;i++)
                    genotype->which_protein[i]=genotype->which_protein[i+1];                
            }
            /*
             * UPDATE protein_pool for protein<which_protein
             */           
            for(i=N_SIGNAL_TF;i<which_protein;i++)
            {
                for(j=0;j<genotype->protein_pool[i][0][0];j++)
                    genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]<which_gene)?genotype->protein_pool[i][1][j]:genotype->protein_pool[i][1][j]-1;
            }            
            break;
        case 'd': /*a gene duplication*/
            /* add it to protein_pool, but do not change nproteins*/    
            genotype->protein_pool[which_protein][1][genotype->protein_pool[which_protein][0][0]]=genotype->ngenes; //append newly duplicated gene to the end
            genotype->protein_pool[which_protein][0][0]++;             
            /*update which_protein*/          
            genotype->which_protein[genotype->ngenes]=which_protein;                         
            break;
        case 'c': /*mutation in tf binding seq, creating a new tf*/
            /* remove this copy of gene from the original protein_pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            /*shift to delete which_gene*/
            for(;i<genotype->protein_pool[which_protein][0][0]-1;i++) 
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            /*one less gene copy to encoding which_protein*/
            genotype->protein_pool[which_protein][0][0]--;                               
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins; //put the new protein to the end
            genotype->protein_pool[genotype->nproteins][0][0]=1;
            genotype->protein_pool[genotype->nproteins][1][0]=which_gene;
            /* make Kd for the new protein*/
            genotype->Kd[genotype->nproteins]=genotype->Kd[which_protein];  
            /* update activator or repressor numbers, and protein_identity*/
            if(genotype->protein_identity[which_protein][0]==ACTIVATOR) //mutation to binding seq does not change the identity of a tf
                genotype->N_act++;
            else
                genotype->N_rep++;
            genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];  
            /* Does this gene encodes effector? */                      
            if(genotype->protein_identity[which_protein][1]!=-1) /*Yes!.-1 for non-output proteins*/
            {
                /*the new protein still encodes effector*/
                genotype->protein_identity[genotype->nproteins][1]=genotype->n_output_proteins;
                /*effector is now encoded by an extra protein*/
                genotype->output_protein_id[genotype->n_output_proteins]=genotype->nproteins;                   
                genotype->n_output_proteins++;
            }
            else
            {
                genotype->protein_identity[genotype->nproteins][1]=-1;
            }              
            /* finally, update protein numbers*/
            genotype->nproteins++;            
            break;            
        case 'e': /*mutation in the identity of a TF, creating a new tf*/
            /* remove this copy of gene from the original protein*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0]-1;i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            }
            genotype->protein_pool[which_protein][0][0]--;                       
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins; 
            genotype->protein_pool[genotype->nproteins][0][0]=1;
            genotype->protein_pool[genotype->nproteins][1][0]=which_gene;
            /* update Kd*/
            genotype->Kd[genotype->nproteins]=genotype->Kd[which_protein];
            /* update protein_identity*/
            if(genotype->protein_identity[which_protein][0]==ACTIVATOR) 
                genotype->N_rep++;  /* an activator turns into a repressor */
            else
                genotype->N_act++;
            genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];
            /* Does this gene encodes effector? */                      
            if(genotype->protein_identity[which_protein][1]!=-1) //-1 for non-output proteins
            {
                /*the new protein still encodes effector*/
                genotype->protein_identity[genotype->nproteins][1]=genotype->n_output_proteins;
                /*effector is now encoded by an extra protein*/
                genotype->output_protein_id[genotype->n_output_proteins]=genotype->nproteins;                   
                genotype->n_output_proteins++;
            }
            else
            {
                genotype->protein_identity[genotype->nproteins][1]=-1;
            }  
            /* finally, update protein numbers*/
            genotype->nproteins++;           
            break;            
        case 'f': /*mutation in tf Kd. We only call update_protein_pool if >1 gene 
                   *copies encode the same TF, in which case a new tf is created*/
            /* remove this copy of gene from the original protein pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0]-1;i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; /* rearrange data array */
            }
            genotype->protein_pool[which_protein][0][0]--; 
            /* create a new protein and link it to this gene*/
            genotype->which_protein[which_gene]=genotype->nproteins; 
            genotype->protein_pool[genotype->nproteins][0][0]=1;
            genotype->protein_pool[genotype->nproteins][1][0]=which_gene;            
            /* update protein_identity*/
            if(genotype->protein_identity[which_protein][0]==ACTIVATOR) 
                genotype->N_act++;  
            else
                genotype->N_rep++;
            genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];  
            /* Does this gene encodes effector? */                      
            if(genotype->protein_identity[which_protein][1]!=-1) //-1 for non-output proteins
            {
                /*the new protein still encodes effector*/
                genotype->protein_identity[genotype->nproteins][1]=genotype->n_output_proteins;
                /*effector is now encoded by an extra protein*/
                genotype->output_protein_id[genotype->n_output_proteins]=genotype->nproteins;                   
                genotype->n_output_proteins++;
            }
            else
            {
                genotype->protein_identity[genotype->nproteins][1]=-1;
            }
            /* finally, update protein numbers*/
            genotype->nproteins++;
            /* NOTE: this mutation does not change the number of genes*/
            break;
        case 't': /*an effector gene mutate into regular TF*/
            if(genotype->protein_pool[which_protein][0][0]>1) /*if there are more than 1 gene encoding the TF*/
            {
                /*the mutation creates a new non-output protein*/
                gene_id=0;
                while(genotype->protein_pool[which_protein][1][gene_id]!=which_gene) gene_id++;
                /*remove this gene from the original protein pool*/
                for(i=gene_id;i<genotype->protein_pool[which_protein][0][0];i++)
                    genotype->protein_pool[which_protein][1][i]=genotype->protein_pool[which_protein][1][i+1];
                genotype->protein_pool[which_protein][0][0]--;
                /*make a new protein*/
                genotype->protein_pool[genotype->nproteins][0][0]=1;
                genotype->protein_pool[genotype->nproteins][1][0]=which_gene;
                genotype->which_protein[which_gene]=genotype->nproteins;
                genotype->Kd[genotype->nproteins]=genotype->Kd[which_protein];
                genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];
                genotype->protein_identity[genotype->nproteins][1]=-1;    
                if(genotype->protein_identity[genotype->nproteins][0]==ACTIVATOR)
                    genotype->N_act++;
                else
                    genotype->N_rep++;
                genotype->nproteins++;
            }
            else /*just remove the gene from output_protein_id*/
            {
                which_output=genotype->protein_identity[which_protein][1];
                for(i=0;i<genotype->n_output_proteins;i++)
                    genotype->protein_identity[genotype->output_protein_id[i]][1]=(genotype->protein_identity[genotype->output_protein_id[i]][1]<which_output)?
                                                                                    genotype->protein_identity[genotype->output_protein_id[i]][1]:genotype->protein_identity[genotype->output_protein_id[i]][1]-1;
                for(i=which_output;i<genotype->n_output_proteins-1;i++)
                    genotype->output_protein_id[i]=genotype->output_protein_id[i+1];
                
                genotype->protein_identity[which_protein][1]=-1;
                genotype->n_output_proteins--;
            }            
            break;
    }
}


/*To reduce the amount of calculation on the probability of binding distributions, we group gene copies
 *that are created by whole gene duplication. We call such a group a cis-reg cluster because gene copies 
 *in the group should have the same cis-reg sequence. For each cis-reg cluster we only need to calculate
 *the probability of binding distributions once. However, substitutions in cis-reg sequence can create/remove
 *binding sites, therefore we need to check whether a gene copy is still in the original cis-reg cluster 
 *after mutation.We use cisreg_cluster and which_cluster to track the bi-way relation between a gene and 
 *a cis-reg cluster.*/
void update_cisreg_cluster(Genotype *genotype, int which_gene, char mut_type, int new_clusters[NGENES][NGENES], int N_new_clusters, int original_cluster_id)
{
    /*In a cis-reg cluster, gene copies are ordered ascendingly by their ids. There are no empty slots in the list 
     *of gene copies. Empty slots after the list are marked by -1. We do not track the number of gene copies in a
     *cluster. In order to count the number of gene copies correctly, marking the empty slots accurately become important.
     *Therefore, we always set all the slots in a cluster to -1, when deleting or overwriting the cluster. 
     */
    int cisreg_seq_cluster_id, cisreg_seq_cluster_id_copy, i, j, last_cluster; 
    
    if(mut_type=='c')//a mutation to the binding sequence of a TF happened. this can create new clusters.
    {    
        /*reset the original cisreg_cluster*/
        i=0;
        while(genotype->cisreg_cluster[original_cluster_id][i]!=-1)
        {
            genotype->cisreg_cluster[original_cluster_id][i]=-1;
            i++;
        }
        /*assign the first new cluster to the original cluster*/
        i=0;
        while(new_clusters[0][i]!=-1)
        {
            genotype->cisreg_cluster[original_cluster_id][i]=new_clusters[0][i];
            genotype->which_cluster[new_clusters[0][i]]=original_cluster_id;
            i++;
        }
        /*find a empty slot in cisreg_cluster for the rest of new clusters*/    
        last_cluster=N_SIGNAL_TF;
        while(genotype->cisreg_cluster[last_cluster][0]!=-1)last_cluster++;//an empty cluster has all of its slots marked by -1
        /*assign the rest of new clusters into empty slots*/       
        for(i=1;i<N_new_clusters;i++)
        {
            /*reset the empty cluster, just in case*/
            j=0;
            while(genotype->cisreg_cluster[last_cluster][j]!=-1)j++;
            /*assign new clusters*/
            j=0;
            while(new_clusters[i][j]!=-1)
            {
                genotype->cisreg_cluster[last_cluster][j]=new_clusters[i][j];
                genotype->which_cluster[new_clusters[i][j]]=last_cluster;
                j++;
            }
            last_cluster++;
        }        
        return;
    } 
    /*find the cis-reg cluster of which gene*/
    cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
    if(mut_type!='d')/*not a duplication*/
    {   
        if(mut_type!='w') /*not a gene deletion. Substitution or locus_specific_TF_behavior just kicked one gene copy out of a cluster*/
        {
            if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=-1) //if a cluster contains more than 1 gene
            {   
                /*then remove which_gene from the original cluster*/        
                i=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=which_gene) i++;
                /*shift the slots to overwrite which_gene*/
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=-1)
                {
                    genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1];
                    i++;
                }                     
                /* and create a new cluster*/
                last_cluster=N_SIGNAL_TF; // start from the first non-signal-tf gene
                while(genotype->cisreg_cluster[last_cluster][0]!=-1) last_cluster++; 
                genotype->cisreg_cluster[last_cluster][0]=which_gene;
                genotype->which_cluster[which_gene]=last_cluster;
            }
        }
        else /*is gene deletion*/
        {               
            if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=-1) //if a cluster contains more than 1 gene
            {   
                /*then remove which_gene from the original cluster*/   
                /*find which_gene in cluster*/
                i=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=which_gene) i++;
                /*shift to overwrite which_gene, and update ids of the remaining genes*/
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=-1)
                {      
                    /*note that gene copies are ordered ascendingly*/
                    genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1]<0)?-1:(genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1]-1); 
                    i++;
                }               
                /*take care of clusters>cisreg_seq_cluster_id*/
                cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id+1;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=-1)
                {
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=-1)
                    {
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]<which_gene)?genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]:genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]-1;
                        i++;
                    }
                    cisreg_seq_cluster_id_copy++;
                }              
                /*update which_cluster*/
                for(i=which_gene;i<genotype->ngenes-1;i++)                              
                    genotype->which_cluster[i]=genotype->which_cluster[i+1];               
            }
            else//if a cluster contains only 1 gene
            {   
                /*need to shift cisreg_cluster*/ 
                cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=-1)
                {
                    /*reset cluster=cisreg_seq_cluster_id_copy*/
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=-1)
                    {
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=-1;
                        i++;
                    }
                    /*then copy from cluster=cisreg_seq_cluster_id_copy+1*/                    
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]!=-1)
                    {                        
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]<which_gene)?genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]:genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]-1; /*note deletion changes gene thread_ID*/
                        i++;
                    }
                    cisreg_seq_cluster_id_copy++;
                }
                /*shift which_cluster and update cluster ids for gene>=which_gene*/                
                for(i=which_gene;i<genotype->ngenes-1;i++)                                
                    genotype->which_cluster[i]=(genotype->which_cluster[i+1]<cisreg_seq_cluster_id)?genotype->which_cluster[i+1]:genotype->which_cluster[i+1]-1;
                /*update cluster ids for gene<which_gene*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_cluster[i]=(genotype->which_cluster[i]>cisreg_seq_cluster_id)?genotype->which_cluster[i]-1:genotype->which_cluster[i];
            }
            /*update gene ids in clusters<cisreg_seq_cluster_id*/
            for(i=N_SIGNAL_TF;i<cisreg_seq_cluster_id;i++)
            {
                j=0;
                while(genotype->cisreg_cluster[i][j]!=-1)
                {
                    genotype->cisreg_cluster[i][j]=(genotype->cisreg_cluster[i][j]<which_gene)?genotype->cisreg_cluster[i][j]:genotype->cisreg_cluster[i][j]-1;
                    j++;
                }
            }                                            
        }
    }
    else /*gene duplication*/
    { 
        cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
        /*find an empty slot in cisreg_seq_cluster_id*/
        i=0;
        while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=-1) i++;
        /*the duplicate is always appended to the end of gene list. Note that 
         *ngenes has not been updated when update_cisreg_cluster is called.*/
        genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->ngenes;
        /*update which_cluster*/
        genotype->which_cluster[genotype->ngenes]=cisreg_seq_cluster_id;
    }   
}

/****************** end of mutation functions *********************************/


/*
 *Set default values and allocate space for variables in Genotype
 */
void initialize_cache(Genotype *genotype)
{
    int j,k;
    FILE *fperror;  
    /*Initialize variables that applies to loci*/
    for(j=0;j<NGENES;j++)
    {
        genotype->which_protein[j]=-1;                
        genotype->recalc_TFBS[j]=1;
        genotype->which_cluster[j]=-1; 
        genotype->min_N_activator_to_transc[j]=MAX_BINDING+1; /*by default a gene cannot be turned on. 
                                                       *MAX_BINDING is the maximum number of tf that 
                                                       *can bind to a cis-reg sequence.*/
        genotype->Kd[j]=-1.0;
        for(k=0;k<NGENES;k++)        
            genotype->cisreg_cluster[j][k]=-1;
        for(k=0;k<NPROTEINS;k++)
            genotype->locus_specific_TF_behavior[j][k]=NON_TF;
    }    
    for(j=0;j<NGENES;j++)
        genotype->cisreg_cluster[NGENES][j]=-1;
    /* initialize variables that applies to protein */
    for(j=0;j<NPROTEINS;j++)
    {
        genotype->protein_pool[j][0][0]=0;
        for(k=0;k<NGENES;k++)        
            genotype->protein_pool[j][1][k]=-1;
    }
    for(j=0;j<NPROTEINS;j++)
    {
        genotype->protein_identity[j][0]=-1;
        genotype->protein_identity[j][1]=-1;
    }   
    for(j=0;j<MAX_OUTPUT_PROTEINS;j++)        
        genotype->output_protein_id[j]=-1;
    /* alloc space for binding sites*/
    genotype->N_allocated_elements=MAXELEMENTS;
    for(j=0;j<NGENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        if (!(genotype->all_binding_sites[j])) 
        {
            fperror=fopen(error_file,"a+");
            LOG("Failed to allocate space\n");
            fclose(fperror);            
            exit(-1);
        }
    }
    /*Initialize binding sites summary*/
    for(j=N_SIGNAL_TF;j<NGENES;j++)
    {
        genotype->N_act_BS[j]=0;
        genotype->N_rep_BS[j]=0;
        genotype->binding_sites_num[j]=0;
    }    
}


/**
 * Given the fitness of the resident and a mutant, decide whether the mutant can replace the resident
 */
float try_fixation(Genotype *resident, Genotype *mutant, int N_measurement_resident, int N_measurement_mutant, int *fixation, RngStream RS)
{  
    float s, ref;
	s = (mutant->fitness - resident->fitness) / fabs(resident->fitness);
    if(s>=minimal_selection_coefficient)
        *fixation=1;
    else          
        *fixation=0;  

    return s;
}

void replay_mutations(  Genotype *genotype_ori,
                        Genotype *genotype_ori_copy,
                        FILE *file_mutation,
                        Mutation *mut_record,
                        int replay_N_steps,
                        RngStream RS)
{
    int i;
    FILE *fp;  

    
    for(i=0;i<replay_N_steps;i++)
    {
        calc_all_binding_sites(genotype_ori);
//        if(i%OUTPUT_INTERVAL==0)
//            summarize_binding_sites(genotype_ori,i); 
//        find_i1ffl(genotype_ori);
//        print_core_i1ffls(genotype_ori);
        clone_genotype(genotype_ori,genotype_ori_copy);
        fscanf(file_mutation,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                    &(mut_record->which_gene),
                                                    &(mut_record->which_nucleotide), 
                                                    mut_record->nuc_diff,               
                                                    &(mut_record->kinetic_type),
                                                    &(mut_record->kinetic_diff));
        reproduce_mutate(genotype_ori_copy,mut_record,RS);        
        clone_genotype(genotype_ori_copy,genotype_ori); 
//        printf("%d, %c, %f\n",i,mut_record->mut_type,mut_record->kinetic_diff);
    }
    calc_all_binding_sites(genotype_ori);
//    summarize_binding_sites(genotype_ori,i); 
//    find_i1ffl(genotype_ori);
//    print_core_i1ffls(genotype_ori);
}

#if NEUTRAL
void evolve_neutrally(  Genotype *genotype_ori,
                        Genotype *genotype_ori_copy,                      
                        Mutation *mut_record,
                        int max_mut_steps,                       
                        RngStream RS_main)
{
    int i=0;    
    FILE *fp;
    DUPLICATION=1.5e-7;                 
    SILENCING = 1.3e-7+0.2e-7;
    while(i<max_mut_steps)
    {        
        clone_genotype(genotype_ori,genotype_ori_copy);
        mutate(genotype_ori_copy,RS_main,mut_record);  
        fp=fopen(mutation_file,"a+");
        fprintf(fp,"%c %d %d '%s' %d %a\n",
                mut_record->mut_type,    
                mut_record->which_gene,
                mut_record->which_nucleotide,
                mut_record->nuc_diff,
                mut_record->kinetic_type,
                mut_record->kinetic_diff);
        fclose(fp);       
        clone_genotype(genotype_ori_copy,genotype_ori);   
        i++;
        calc_all_binding_sites(genotype_ori);
//        find_ffl(genotype_ori); 
//        print_core_c1ffls(genotype_ori);        
        /*output network topology every OUTPUT_INTERVAL steps*/ 
        if(i%OUTPUT_INTERVAL==0 && i!=0)
            summarize_binding_sites(genotype_ori,i);        
        /*output a summary of simulation every step*/
        output_genotype(genotype_ori, i);
    } 
}

#endif

#if JUST_PLOTTING
void run_plotting(  Genotype *genotype_ori,
                    Genotype *genotype_ori_copy,
                    int init_mRNA[NGENES],
                    float init_protein_number[NGENES],
                    RngStream RS_parallel[N_THREADS],
                    Mutation *mut_record,
                    FILE *fp,
                    int replay_N_steps)
{      
    replay_mutations(genotype_ori, genotype_ori_copy, fp, mut_record, replay_N_steps,RS_parallel[0]); 
    #if !SET_BS_MANUALLY
//         calc_all_binding_sites(genotype_ori);  
    #endif
//    summarize_binding_sites(genotype_ori,1);   
//    exit(0);
    /* conditions under which the phenotype and fitness is measured */    
    env1_t_development=899.9;
    env2_t_development=899.9;
    opt_pulse_duration=30.0;
    sampling_interval=1.0; 
//    saturate_cumulative_response_from_pulse=4.0*Ne_saturate*opt_pulse_duration;
//    tolerable_delay_bf_pulse=40.0;
    duration_of_burn_in_growth_rate=0.0;
    env1_minimal_peak_response=100000.0;
    env2_minimal_peak_response=100000.0;
    env1_response_amplification=10000.0;
    env2_response_amplification=10000.0;
    env1_benefit1=0.1;
    env1_benefit2=1.0;
    env1_max_duration_bias=20.0;            
    env2_benefit1=0.1;
    env2_benefit2=1.0;
    env2_max_duration_bias=20.0;
    env1_signal1_strength=10.0;
    env1_signal2_strength=500.0;
    env2_signal1_strength=10.0;
    env2_signal2_strength=500.0;
    env1_t_signal_on=120.0;    
    env1_t_signal_off=60.0;     
    env2_t_signal_on=120.0;
    env2_t_signal_off=60.0;
    env1_initial_effect_of_effector='d';
    env2_initial_effect_of_effector='d';            
    env1_fixed_effector_effect=0;    
    env2_fixed_effector_effect=0; 
    recalc_new_fitness=1;
    env1_occurence=0.5;
    env2_occurence=0.5;            
    calc_cellular_fitness_plotting(  genotype_ori, 
                                    init_mRNA, 
                                    init_protein_number,
                                    RS_parallel);    
}
#endif

#if PLOT_ALTERNATIVE_FITNESS
void plot_alternative_fitness(  Genotype *genotype_ori,
                                Genotype *genotype_ori_copy,
                                int init_mRNA[NGENES],
                                float init_protein_number[NGENES],
                                RngStream RS_parallel[N_THREADS],
                                Mutation *mut_record,
                                FILE *file_mutation,
                                int replay_N_steps)
{
    int i,j;
    int steps_to_be_recalc[MAX_MUT_STEP],number_of_steps;
    FILE *alternative_fitness,*file_steps_to_be_recalc;   
    env1_t_development=149.9;                     // global variable
    env2_t_development=149.9; 
    env1_signal_strength=10000.0;
    env2_signal_strength=10000.0;
    duration_of_burn_in_growth_rate=0.0;    // global variable    
    env1_t_signal_on=200.0;     
    env1_t_signal_off=0.0;
    env2_t_signal_on=20.0;
    env2_t_signal_off=130.0;
    env1_initial_effect_of_effector='b';    
    env2_initial_effect_of_effector='d'; 
    env1_fixed_effector_effect=0;
    env2_fixed_effector_effect=1;
    recalc_new_fitness=5;                   // global variable, make sure its value is smaller than MAX_RECALC_FITNESS        
    env1_occurence=0.5;                     // global variable
    env2_occurence=0.5;      
    float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];  
    for(i=0;i<replay_N_steps;i++)
    {
        calc_all_binding_sites(genotype_ori);        
        clone_genotype(genotype_ori,genotype_ori_copy);
        fscanf(file_mutation,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                    &(mut_record->which_gene),
                                                    &(mut_record->which_nucleotide), 
                                                    mut_record->nuc_diff,               
                                                    &(mut_record->kinetic_type),
                                                    &(mut_record->kinetic_diff));
        reproduce_mutate(genotype_ori_copy,mut_record,NULL);        
        clone_genotype(genotype_ori_copy,genotype_ori); 
        find_ffl(genotype_ori);
//        print_core_c1ffls(genotype_ori);
//        summarize_binding_sites(genotype_ori,i);
        if(i==20670)
        {
            for(j=0;j<recalc_new_fitness;j++)                    
            {    
                calc_cellular_fitness(   genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,                                        
                                        0,
                                        GR1[j],
                                        GR2[j],
                                        mut_record); 
            }

            calc_genotype_fitness( genotype_ori,
                                &(GR1[0]),
                                &(GR2[0]),
                                recalc_new_fitness);  

            alternative_fitness=fopen("alternative_fitness.txt","a+");
            fprintf(alternative_fitness,"%d %.10f %.10f %.10f %.10f %.10f %.10f\n",i,
                    genotype_ori->fitness,
                    sqrt(genotype_ori->sq_SE_fitness),
                    genotype_ori->avg_GR1,
                    genotype_ori->avg_GR2,
                    sqrt(genotype_ori->sq_SE_GR1),
                    sqrt(genotype_ori->sq_SE_GR2));
            fclose(alternative_fitness);            
        }
    }    
}
#endif

void run_simulation(    Genotype *genotype_ori, 
                        Genotype *genotype_ori_copy, 
                        float init_protein_number[NPROTEINS],
                        int init_mRNA[NGENES],                    
                        Mutation *mut_record, 
                        int init_N_tot_mutations,
                        int init_step,
                        RngStream RS_main,
                        RngStream RS_parallel[N_THREADS])
{
    FILE *fp;
    int max_mut_steps,run_burn_in,N_tot_trials,first_step,end_state=0;    
    first_step=init_step;
    N_tot_trials=init_N_tot_mutations; 
    int i;
    /* first, run burn-in */
    if(BURN_IN_I)
    {
        run_burn_in=1;      
        max_mut_steps=BURN_IN_I;    
        env1_t_development=119.9;
        env2_t_development=119.9;                 // global variable
        duration_of_burn_in_growth_rate=0.0;// global variable        
        env1_response_amplification=5.0;
        env2_response_amplification=5.0;
        env1_benefit1=1.0;
        env1_benefit2=2.0;
        env1_max_duration_bias=20.0;            
        env2_benefit1=1.0;
        env2_benefit2=2.0;
        env2_max_duration_bias=20.0;
        env1_signal1_strength=1000.0;
        env1_signal2_strength=1000.0;
        env2_signal1_strength=1000.0;
        env2_signal2_strength=1000.0;    
        env1_t_signal_on=70.0;    
        env1_t_signal_off=20.0;     
        env2_t_signal_on=70.0;
        env2_t_signal_off=20.0;
        env1_initial_effect_of_effector='b';
        env2_initial_effect_of_effector='b';
        env1_fixed_effector_effect=0;    
        env2_fixed_effector_effect=0;            // global variable
        recalc_new_fitness=5;               // global variable, make sure its value is smaller than MAX_RECALC_FITNESS
        env1_occurence=0.5;                 // global variable
        env2_occurence=0.5;                 // global variable
        DUPLICATION=1.5e-7;                 
        SILENCING = 3.0e-7;
        miu_ACT_TO_INT_RATE=1.27;
        miu_mRNA_decay=-1.49;       
        miu_translation_init=0.408;
      
        fp=fopen(RuntimeSumm,"a+");
        fprintf(fp,"**********Burn-in_I conditions**********\n");
        fprintf(fp,"BURN_IN_I=%d\n",BURN_IN_I);                
        fprintf(fp,"N_REPLICATES=%d\n",N_REPLICATES);        
        fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
        fprintf(fp,"env1_t_development=%f, env2_t_development=%f\n",env1_t_development,env2_t_development);
        fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);        
        fprintf(fp,"env1: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env1_t_signal_on, env1_t_signal_off, env1_initial_effect_of_effector, env1_fixed_effector_effect, env1_occurence);
        fprintf(fp,"env2: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env2_t_signal_on, env2_t_signal_off, env2_initial_effect_of_effector, env2_fixed_effector_effect, env2_occurence);
        fclose(fp); 
        
        end_state=evolve_N_steps( genotype_ori, 
                        genotype_ori_copy,
                        &first_step, 
                        max_mut_steps, 
                        &N_tot_trials,                      
                        init_protein_number,
                        init_mRNA,                         
                        mut_record, 
                        RS_main,
                        RS_parallel,
                        run_burn_in);         
       
        fp=fopen(RuntimeSumm,"a+");
        fprintf(fp,"Burn_in_I completes after the %dth step.\n",first_step);
        fclose(fp);
        
        if(end_state==-1)
            return;  
        if(init_step<BURN_IN_I)
        {
            env1_occurence=0.33;                     
            env2_occurence=0.67;   
            float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];
            for(i=1;i<recalc_new_fitness;i++)                    
                calc_cellular_fitness(  genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,                                    
                                        BURN_IN_I,
                                        GR1[i],
                                        GR2[i],
                                        mut_record);  
            calc_genotype_fitness( genotype_ori,
                                    &(GR1[0]),
                                    &(GR2[0]),
                                    recalc_new_fitness); 
            calc_all_binding_sites(genotype_ori); 

            /*calculate the number of c1-ffls every step*/
            find_i1ffl(genotype_ori); 
            print_core_i1ffls(genotype_ori); 
            /*output network topology every OUTPUT_INTERVAL steps*/
            if(i%OUTPUT_INTERVAL==0 && i!=0) 
                summarize_binding_sites(genotype_ori,i);   

            /*output a summary of simulation every step*/
            output_genotype(genotype_ori, BURN_IN_I);        
            /* output rng seeds*/
    #if OUTPUT_RNG_SEEDS
            unsigned long seeds[6];
            RngStream_GetState(RS_main,seeds);
            fp=fopen("RngSeeds.txt","a+");
            fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]);            
            for(i=0;i<N_THREADS;i++)
            {
                RngStream_GetState(RS_parallel[i],seeds);
                fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]); 
            }
            fprintf(fp,"\n");
            fclose(fp);
    #endif
            /*output precise fitness*/
            fp=fopen("precise_fitness.txt","a+");
            fprintf(fp,"%d %a %a %a %a %a %a\n",N_tot_trials,                                                                       
                                                genotype_ori->fitness,
                                                genotype_ori->sq_SE_fitness,
                                                genotype_ori->avg_GR1,
                                                genotype_ori->avg_GR2,
                                                genotype_ori->sq_SE_GR1,
                                                genotype_ori->sq_SE_GR2);
            fclose(fp);        
            /* marks the last step at which all state of the program has been output*/
            fp=fopen("saving_point.txt","w");
            fprintf(fp,"%d %d\n",i,N_tot_trials);
            fclose(fp);  
        }
    }    
    
    /* post-burn-in simulations*/
    run_burn_in=0;
    max_mut_steps=MAX_MUT_STEP;    
    env1_t_development=119.9;
    env2_t_development=119.9;
    opt_pulse_duration=10.0;
    sampling_interval=1.0;          
//    saturate_pulse_amplitude=Ne_saturate;
//    sd_opt_pulse_duration=10.0;
//    saturate_cumulative_response_from_pulse=50.0*Ne_saturate*opt_pulse_duration;
//    tolerable_delay_bf_pulse=40.0;
    duration_of_burn_in_growth_rate=0.0;
    env1_minimal_peak_response=10000.0;
    env2_minimal_peak_response=10000.0;
    env1_response_amplification=1000.0;
    env2_response_amplification=1000.0;
    env1_benefit1=0.1;
    env1_benefit2=1.0;
    env1_max_duration_bias=20.0;            
    env2_benefit1=0.1;
    env2_benefit2=1.0;
    env2_max_duration_bias=20.0;
    env1_signal1_strength=10.0;
    env1_signal2_strength=500.0;
    env2_signal1_strength=500.0;
    env2_signal2_strength=25000.0;
    env1_t_signal_on=300.0;    
    env1_t_signal_off=60.0;     
    env2_t_signal_on=300.0;
    env2_t_signal_off=60.0;
    env1_initial_effect_of_effector='d';
    env2_initial_effect_of_effector='d';
    env1_fixed_effector_effect=0;    
    env2_fixed_effector_effect=0; 
    recalc_new_fitness=5; // make sure its value is smaller than MAX_RECALC_FITNESS
    env1_occurence=0.5;
    env2_occurence=0.5;                       // global variable
    DUPLICATION=1.5e-7;                 
    SILENCING = 1.5e-7;
    miu_ACT_TO_INT_RATE=1.27;
    miu_mRNA_decay=-1.49;       
    miu_translation_init=0.408;
    
//     for(i=N_SIGNAL_TF;i<genotype_ori->ngenes;i++)
//     {
//         if(genotype_ori->which_protein[i]==genotype_ori->nproteins-1)
//             genotype_ori->min_act_to_transc[i]=2;
//     }

    fp=fopen(RuntimeSumm,"a+");
    fprintf(fp,"**********Post-burn-in conditions**********\n");
    fprintf(fp,"second phase steps=%d\n",max_mut_steps);                
    fprintf(fp,"N_replicates=%d\n",N_REPLICATES);        
    fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
    fprintf(fp,"env1_t_development=%f,env2_t_development=%f\n",env1_t_development,env2_t_development);        
    fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);         
//    fprintf(fp,"env1: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env1_t_signal_on, env1_t_signal_off, env1_initial_effect_of_effector, env1_fixed_effector_effect, env1_occurence);
//    fprintf(fp,"env2: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env2_t_signal_on, env2_t_signal_off, env2_initial_effect_of_effector, env2_fixed_effector_effect, env2_occurence);       
    fprintf(fp,"Background signal strength=%f\n",background_signal_strength);
//    fprintf(fp,"Signal off strength=%f, env1 signal on strength=%f, env2 signal on strength=%f \n",signal_off_strength,env1_signal_strength,env2_signal_strength);
//    fprintf(fp,"env1 init effecto effect %c, fixed? %d. env2 init effecto effect %c, fixed? %d.\n",env1_initial_effect_of_effector,env1_fixed_effector_effect,env2_initial_effect_of_effector,env2_fixed_effector_effect);
    fclose(fp);
    
    end_state=evolve_N_steps( genotype_ori, 
                    genotype_ori_copy,
                    &first_step, 
                    max_mut_steps, 
                    &N_tot_trials,                   
                    init_protein_number,
                    init_mRNA,                 
                    mut_record, 
                    RS_main,
                    RS_parallel,
                    run_burn_in); 
#if !SET_BS_MANUALLY
    calc_all_binding_sites(genotype_ori);
    summarize_binding_sites(genotype_ori,max_mut_steps); /*snapshot of the final distribution binding sites */
#endif
}

void continue_simulation(   Genotype *genotype_ori, 
                            Genotype *genotype_ori_copy,                            
                            int replay_N_steps,                       
                            float init_protein_number[NPROTEINS],
                            int init_mRNA[NGENES],                
                            Mutation *mut_record, 
                            RngStream RS_main,
                            RngStream RS_parallel[N_THREADS])
{
    int i,j,N_tot_mutations;    
    unsigned long rng_seeds[N_THREADS+1][6];
    char buffer[200]; 
    FILE *fp,*fperror;
    
    
    /*delete the incomplete lines in the output files*/
    tidy_output_files(output_file,mutation_file);
    
    /* set genotype based on previous steps*/
    fp=fopen(mutation_file,"r");
    if(fp!=NULL)
        replay_mutations(genotype_ori, genotype_ori_copy, fp, mut_record, replay_N_steps,RS_main);
    else
    {
        fperror=fopen(error_file,"a+");
        LOG("cannot open mutation_file\n");
        fclose(fperror);
        exit(-2);
    }
    fclose(fp);

    /* set random number seeds*/
    fp=fopen("RngSeeds.txt","r");
    if(fp!=NULL)
    {
        for(i=0;i<replay_N_steps;i++)
        {
            for(j=0;j<N_THREADS;j++)        
            {
                fscanf(fp,"%lu %lu %lu %lu %lu %lu ",
                        &(rng_seeds[j][0]),
                        &(rng_seeds[j][1]),
                        &(rng_seeds[j][2]),
                        &(rng_seeds[j][3]),
                        &(rng_seeds[j][4]),
                        &(rng_seeds[j][5]));
            }
            fscanf(fp,"%lu %lu %lu %lu %lu %lu \n",
                    &(rng_seeds[N_THREADS][0]),
                    &(rng_seeds[N_THREADS][1]),
                    &(rng_seeds[N_THREADS][2]),
                    &(rng_seeds[N_THREADS][3]),
                    &(rng_seeds[N_THREADS][4]),
                    &(rng_seeds[N_THREADS][5]));
        }
    }
    else
    {
        fperror=fopen(error_file,"a+");
        LOG("cannot open RngSeeds.txt\n");
        fclose(fperror);
        exit(-2);
    }
    fclose(fp);
    RngStream_SetSeed(RS_main,rng_seeds[0]);
    for(i=0;i<N_THREADS;i++)
        RngStream_SetSeed(RS_parallel[i],rng_seeds[i+1]);
    
    /* load fitness,N_tot_mutations*/
    fp=fopen("precise_fitness.txt","r");
    if(fp!=NULL)
    {  
        for(i=0;i<replay_N_steps-1;i++)
            fgets(buffer,200,fp);
            
        fscanf(fp,"%d %a %a %a %a %a %a\n",&N_tot_mutations,                
                                            &(genotype_ori->fitness),
                                            &(genotype_ori->sq_SE_fitness),
                                            &(genotype_ori->avg_GR1),
                                            &(genotype_ori->avg_GR2),
                                            &(genotype_ori->sq_SE_GR1),
                                            &(genotype_ori->sq_SE_GR2));
    }
    else
    {
        fperror=fopen(error_file,"a+");
        LOG("cannot open precise_fitness.txt\n");
        fclose(fperror);
        exit(-2);
    }        
    fclose(fp);
    /*continue running simulation*/
    run_simulation( genotype_ori, 
                    genotype_ori_copy,  
                    init_protein_number,
                    init_mRNA,                           
                    mut_record, 
                    N_tot_mutations,
                    replay_N_steps+1,
                    RS_main,
                    RS_parallel);    
}

void calc_genotype_fitness(Genotype *genotype,
                            float (*GR1)[N_REPLICATES],
                            float (*GR2)[N_REPLICATES],
                            int N_recalc_fitness)
{
    float avg_GR1=0.0;
    float avg_GR2=0.0;   
    float sum_sq_diff_GR1=0.0;
    float sum_sq_diff_GR2=0.0;   
    float sum_sq_diff_mean_GR=0.0;
    float diff_GR1,diff_GR2,sq_SE_GR1,sq_SE_GR2;    
    int counter=0;
    int i,j;

    for(i=0;i<N_recalc_fitness;i++)
    {
        for(j=0;j<N_REPLICATES;j++)
        {

            avg_GR1+=GR1[i][j];
            avg_GR2+=GR2[i][j];
            genotype->fitness_measurement[counter]=env1_occurence*GR1[i][j]+env2_occurence*GR2[i][j];
            counter++;
        }
    }
    avg_GR1=avg_GR1/(N_recalc_fitness*N_REPLICATES);
    avg_GR2=avg_GR2/(N_recalc_fitness*N_REPLICATES);    
    for(i=0;i<N_recalc_fitness;i++)
    {
        for(j=0;j<N_REPLICATES;j++)
        {
            diff_GR1=GR1[i][j]-avg_GR1;
            diff_GR2=GR2[i][j]-avg_GR2;
            sum_sq_diff_GR1+=pow(diff_GR1,2.0);
            sum_sq_diff_GR2+=pow(diff_GR2,2.0);
            sum_sq_diff_mean_GR+=pow(diff_GR1*env1_occurence+diff_GR2*env2_occurence,2.0);

        }
    }
    sq_SE_GR1=sum_sq_diff_GR1/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));
    sq_SE_GR2=sum_sq_diff_GR2/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));     
    genotype->avg_GR1=avg_GR1;
    genotype->avg_GR2=avg_GR2;
    genotype->sq_SE_GR1=sq_SE_GR1*N_recalc_fitness*N_REPLICATES;
    genotype->sq_SE_GR2=sq_SE_GR2*N_recalc_fitness*N_REPLICATES;
    genotype->fitness=env1_occurence*avg_GR1+env2_occurence*avg_GR2;
    genotype->sq_SE_fitness=sum_sq_diff_mean_GR/(N_recalc_fitness*N_REPLICATES-1)/(N_recalc_fitness*N_REPLICATES);
}

int evolve_N_steps(Genotype *genotype_ori, 
                    Genotype *genotype_ori_copy,
                    int *init_step, 
                    int max_steps, 
                    int *N_tot_trials,                 
                    float init_protein_number[NPROTEINS],
                    int init_mRNA[NGENES],                 
                    Mutation *mut_record, 
                    RngStream RS_main,
                    RngStream RS_parallel[N_THREADS],
                    int run_burn_in)
{
    int i,j;
    int fixation;
    int N_trials;
    float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];
    float score;  
    FILE *fp;
    /*If burn_in evolution is set, allow burn_in to end early when certain criteria is met*/
#if QUICK_BURN_IN  
    int ready_to_evolve;
    if(run_burn_in)
        ready_to_evolve=0;
#endif
    for(i=(*init_step);i<=max_steps;i++)
    {             
        fixation=0;      
        N_trials=0;          
        
        if(*N_tot_trials>MAX_MUTATIONS) /*This is an alternative termination condition, although we usually use MAX_MUT_STEP as termination condition*/
        {
            *init_step=i;
            return 0;
        }
        while(1) /*try mutations until one replaces the current genotype*/
        {	
            N_trials++;
            (*N_tot_trials)++;
            if(N_trials>MAX_TRIALS) /*Tried too many mutation in one step.*/
            {
                fp=fopen(output_file,"a+");                              
                fprintf(fp,"Tried %d mutations, yet none could fix\n",MAX_TRIALS);
                fclose(fp); 
                summarize_binding_sites(genotype_ori,i);
                return -1;
            }
            /*do mutation on a copy of the current genotype*/
            clone_genotype(genotype_ori,genotype_ori_copy); 
            mutate(genotype_ori_copy,RS_main,mut_record);
#if OUTPUT_MUTANT_DETAILS
            /*record every mutation*/
            fp=fopen("MUT_Detail.txt","a+");
            fprintf(fp,"%d %d %c %d %d '%s' %d %a\n",
                    i,
                    *N_tot_trials,
                    mut_record->mut_type,
                    mut_record->which_gene,
                    mut_record->which_nucleotide,
                    mut_record->nuc_diff,
                    mut_record->kinetic_type,
                    mut_record->kinetic_diff);
            fclose(fp); 
#endif
            #if !SET_BS_MANUALLY		
                calc_all_binding_sites(genotype_ori_copy);
            #endif
            MAXELEMENTS=genotype_ori_copy->N_allocated_elements;
            /*calculate the fitness of the mutant at low resolution*/
            calc_cellular_fitness(  genotype_ori_copy,
                                    init_mRNA,
                                    init_protein_number,
                                    RS_parallel,                                   
                                    i,
                                    GR1[0],
                                    GR2[0],
                                    mut_record);
            calc_genotype_fitness( genotype_ori_copy,
                                    &(GR1[0]),
                                    &(GR2[0]),
                                    1); // N_recalc_fitness=1 implicitly when calculating the fitness of a mutant.
#if OUTPUT_MUTANT_DETAILS
            /*record low-resolution fitness*/
            fp=fopen("Mut_detail_fitness.txt","a+");
            fprintf(fp,"%.10f %.10f\n",
                    genotype_ori_copy->fitness,
                    genotype_ori_copy->sq_SE_fitness);
            fclose(fp);        
#endif
            /*Can the mutant replace the current genotype?*/
            score=try_fixation(genotype_ori, genotype_ori_copy, N_REPLICATES*recalc_new_fitness, N_REPLICATES, &fixation, RS_main); 
            /*If yes, record relevant info*/
            if(fixation==1)
            {                    
                fp=fopen(output_file,"a+");
                fprintf(fp,"%d %d %d %c %f ",i, *N_tot_trials, N_trials,mut_record->mut_type,score);
                fclose(fp);
                fp=fopen(mutation_file,"a+");
                fprintf(fp,"%c %d %d '%s' %d %a\n",
                        mut_record->mut_type,    
                        mut_record->which_gene,
                        mut_record->which_nucleotide,
                        mut_record->nuc_diff,
                        mut_record->kinetic_type,
                        mut_record->kinetic_diff);
                fclose(fp);
                break;
            }
        }
        /*replace the current genotype by overwriting it*/
        clone_genotype(genotype_ori_copy,genotype_ori);        
        #if !SET_BS_MANUALLY
            calc_all_binding_sites(genotype_ori); 
        #endif
        /*increase the accuracy of the fitness of the new genotype*/ 
        if(i!=BURN_IN_I)
        {
            for(j=1;j<recalc_new_fitness;j++)                    
                calc_cellular_fitness(   genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,                                    
                                        0,
                                        GR1[j],
                                        GR2[j],
                                        mut_record);  
            calc_genotype_fitness( genotype_ori,
                                    &(GR1[0]),
                                    &(GR2[0]),
                                    recalc_new_fitness); 
           
            #if !SET_BS_MANUALLY
                calc_all_binding_sites(genotype_ori); 
            #endif            
            /*calculate the number of c1-ffls every step*/
            find_i1ffl(genotype_ori); 
            print_core_i1ffls(genotype_ori); 
            /*output network topology every OUTPUT_INTERVAL steps*/
            if(i%OUTPUT_INTERVAL==0 && i!=0) 
                summarize_binding_sites(genotype_ori,i);   
       
            /*output a summary of simulation every step*/
            output_genotype(genotype_ori, i);        
            /* output rng seeds*/
            #if OUTPUT_RNG_SEEDS
                unsigned long seeds[6];
                RngStream_GetState(RS_main,seeds);
                fp=fopen("RngSeeds.txt","a+");
                fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]);            
                for(j=0;j<N_THREADS;j++)
                {
                    RngStream_GetState(RS_parallel[j],seeds);
                    fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]); 
                }
                fprintf(fp,"\n");
                fclose(fp);
            #endif
            /*output precise fitness*/
            fp=fopen("precise_fitness.txt","a+");
            fprintf(fp,"%d %a %a %a %a %a %a\n",*N_tot_trials,                                                                       
                                                genotype_ori->fitness,
                                                genotype_ori->sq_SE_fitness,
                                                genotype_ori->avg_GR1,
                                                genotype_ori->avg_GR2,
                                                genotype_ori->sq_SE_GR1,
                                                genotype_ori->sq_SE_GR2);
            fclose(fp);        
            /* marks the last step at which all state of the program has been output*/
            fp=fopen("saving_point.txt","w");
            fprintf(fp,"%d %d\n",i,*N_tot_trials);
            fclose(fp);   
        }
        /* determine if QUICK_BURN_IN criteria is met*/
        #if QUICK_BURN_IN
            if(run_burn_in)
            {            
                if(genotype_ori->fitness>0.3*env1_occurence+0.7*env2_occurence)
                    ready_to_evolve++;
                else
                    ready_to_evolve=0;
                if(ready_to_evolve>=10)
                    break;
            }
        #endif
    } 
    *init_step=i;
    return 0;
}

#if SET_BS_MANUALLY
void draw_network(Genotype* genotype)
{
    /* This allows the user to draw a network with certain topology. Drawing is basically 
     * specifying the number and identities of binding sites on a cis-regulatory sequence,
     * but the actual sequence of the cis-regulatory region and the binding sequences of
     * TFs are still randomly generated (the actually topology is overridden by the user).
     * This means if calc_all_binding_sites is called, the network will revert to its
     * actual topology. The user also needs to specify other parameters that are associated
     * with network topology, such as protein_pool, which_protein, ngenes, nproteins, etc. 
     * Some of these parameters needs to altered in initialize_genotype. In addition, kinetic
     * parameters of each gene are not affected by the "drawing"; they are still generated in
     * initialize_genotype_fixed. 
     */
    
    /* drawing a c1-FFL with two copies of the selection gene*/   
//    int N_BS_per_gene[5]={0,0,1,2,2};
//    int tf_id_per_site[5][5]={{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{1,-1,-1,-1,-1},{2,2,-1,-1,-1},{2,3,-1,-1,-1}};
//    int hindrance_per_site[5][5]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
//    float Kd_per_site[5][5]={{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8}};
    int N_BS_per_gene[6]={0,0,1,2,2,2};
    int tf_id_per_site[6][5]={{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{1,-1,-1,-1,-1},{2,2,-1,-1,-1},{2,3,-1,-1,-1},{2,3,-1,-1,-1}};
    int hindrance_per_site[6][5]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    float Kd_per_site[6][5]={{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8},{1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8}};    
    genotype->ngenes++; 
    genotype->which_protein[5]=4;
    genotype->protein_pool[4][0][0]=2;
    genotype->protein_pool[4][1][1]=5;
    genotype->mRNAdecay[5]=genotype->mRNAdecay[4];
    genotype->proteindecay[5]=genotype->proteindecay[4];
    genotype->translation[5]=genotype->translation[4];
    genotype->pic_disassembly[5]=genotype->pic_disassembly[4];
    genotype->cisreg_cluster[4][1]=5;
    genotype->which_cluster[5]=4;
    genotype->min_act_to_transc[5]=genotype->min_act_to_transc[4];
    genotype->min_act_to_transc[3]=2;
    set_binding_sites(genotype,N_BS_per_gene,tf_id_per_site,hindrance_per_site,Kd_per_site);
}
#endif

int init_run_pop(unsigned long int seeds[6], int CONTINUE)
{  
    int i,j,k,buffer_int;
    Genotype genotype_ori;
    Genotype genotype_ori_copy;    
    int init_mRNA[NGENES]; 
    float init_protein_number[NGENES];   
    int init_step=0;
    Mutation mut_record;
    FILE *fp;    
    RngStream RS_main,RS_parallel[N_THREADS];    
    /*create threads*/
    omp_set_num_threads(N_THREADS);    
    /* initialize random number seeds*/
    RngStream_SetPackageSeed(seeds);    
    RS_main=RngStream_CreateStream("Main");
    for(i=0; i < N_THREADS; i++)
        RS_parallel[i]=RngStream_CreateStream("");
    /* set initial numbers of mRNA and protein*/
    for(i=N_SIGNAL_TF; i < NGENES; i++) /* loop through tf genes*/
    {
        init_mRNA[i]=0;
        init_protein_number[i]=0.0;
    }
    /* initialize genotype */
    initialize_cache(&genotype_ori);
    initialize_cache(&genotype_ori_copy);    
    initialize_genotype(&genotype_ori, RS_main);
    #if SET_BS_MANUALLY  /*manually input binding sites info below*/    
        draw_network(&genotype_ori);
    #endif 
    genotype_ori_copy.ngenes=genotype_ori.ngenes;   
    genotype_ori_copy.nproteins=genotype_ori.nproteins; 
    
    /* initialize mut_record */
    mut_record.kinetic_diff=0.0;
    mut_record.kinetic_type=-1;
    mut_record.mut_type='\0';
    mut_record.nuc_diff[0]='\0';
    mut_record.nuc_diff[1]='\0';
    mut_record.nuc_diff[2]='\0';
    mut_record.which_gene=-1;
    mut_record.which_nucleotide=-1;
    /*different modes of simulation*/    
    if(CONTINUE) /* continue a simulation from a previously saved state*/
    {
        int replay_N_steps=0;
        fp=fopen("saving_point.txt","r");
        fscanf(fp,"%d %d",&replay_N_steps,&buffer_int);
        fclose(fp);
        /*Load external signal profile*/
        fp=fopen("noisy_signal.txt","r");
        if(fp!=NULL)
        {        
            for(j=0;j<200;j++)
            {
                for(k=0;k<15;k++)
                {           
                    fscanf(fp,"%f\n",&(signal_profile_matrix[0][j][k]));
                    for(i=1;i<N_THREADS;i++)
                        signal_profile_matrix[i][j][k]=signal_profile_matrix[0][j][k];                    
                }
            }
            fclose(fp);
        } 
        if(replay_N_steps!=0)
        {                      
            continue_simulation(&genotype_ori, 
                                &genotype_ori_copy,                                
                                replay_N_steps,                                
                                init_protein_number,
                                init_mRNA,                                
                                &mut_record, 
                                RS_main,
                                RS_parallel);                                                                    
            
        }  
    }
    else /* otherwise the simulation starts over from beginning*/
    {   
        /* record the initial network topology*/
        summarize_binding_sites(&genotype_ori,init_step); /*snapshot of the initial (0) distribution binding sites */   
        find_i1ffl(&genotype_ori); 
        print_core_i1ffls(&genotype_ori);

    #if JUST_PLOTTING 
        fp=fopen("MUT.txt","r");    
        if(fp!=NULL)
        {
            printf("LOAD MUTATION RECORD SUCCESSFUL!\n"); 
            run_plotting(   &genotype_ori,
                            &genotype_ori_copy,
                            init_mRNA,
                            init_protein_number,
                            RS_parallel,
                            &mut_record,                                
                            fp,
                            50000);
            fclose(fp);
        }
    #elif PLOT_ALTERNATIVE_FITNESS
        fp=fopen("MUT.txt","r");    
        if(fp!=NULL)
        {
            printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
            plot_alternative_fitness(&genotype_ori,&genotype_ori_copy,init_mRNA,init_protein_number,
                                        RS_parallel,&mut_record,fp,20671);
            fclose(fp);
        }
    #else     
        if(!SKIP_INITIAL_GENOTYPE)/* get the fitness of the initial genotype */ 
        {                     
            env1_t_development=119.9;
            env2_t_development=119.9;
            opt_pulse_duration=10.0;
            sampling_interval=1.0;          
//            saturate_pulse_amplitude=Ne_saturate;
//            sd_opt_pulse_duration=30.0;
//            saturate_cumulative_response_from_pulse=50.0*Ne_saturate*opt_pulse_duration;
//            tolerable_delay_bf_pulse=40.0;
            duration_of_burn_in_growth_rate=0.0;
            env1_minimal_peak_response=10000.0;
            env2_minimal_peak_response=10000.0;
            env1_response_amplification=1000.0;
            env2_response_amplification=1000.0;
            env1_benefit1=0.1;
            env1_benefit2=1.0;
            env1_max_duration_bias=20.0;            
            env2_benefit1=0.1;
            env2_benefit2=1.0;
            env2_max_duration_bias=20.0;
            env1_signal1_strength=10.0;
            env1_signal2_strength=500.0;
            env2_signal1_strength=500.0;
            env2_signal2_strength=25000.0;
            env1_t_signal_on=300.0;    
            env1_t_signal_off=60.0;     
            env2_t_signal_on=300.0;
            env2_t_signal_off=60.0;
            env1_initial_effect_of_effector='d';
            env2_initial_effect_of_effector='d';
            env1_fixed_effector_effect=0;    
            env2_fixed_effector_effect=0; 
            recalc_new_fitness=5; // make sure its value is smaller than MAX_RECALC_FITNESS
            env1_occurence=0.5;
            env2_occurence=0.5;    
            float GR1[recalc_new_fitness][N_REPLICATES],GR2[recalc_new_fitness][N_REPLICATES];
            /*Load external signal profile if there is one*/
            fp=fopen("noisy_signal.txt","r");
            if(fp!=NULL)
            {        
                for(j=0;j<200;j++)
                {
                    for(k=0;k<15;k++)
                    {           
                        fscanf(fp,"%f\n",&(signal_profile_matrix[0][j][k]));
                        for(i=1;i<N_THREADS;i++)
                            signal_profile_matrix[i][j][k]=signal_profile_matrix[0][j][k];                    
                    }
                }
                fclose(fp);
            } 
            for(i=0;i<recalc_new_fitness;i++)      
            {
                calc_cellular_fitness(  &genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,                                       
                                        0,
                                        GR1[i],
                                        GR2[i],
                                        &mut_record);                
            }
                                        
            calc_genotype_fitness(&genotype_ori,&(GR1[0]),&(GR2[0]),recalc_new_fitness);
            fp=fopen(output_file,"a+");
            fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step fixed_mutation score fitness se_fitness avg_GR1 avg_GR2 std_GR1 std_GR2 N_genes N_proteins N_activator N_repressor\n");
            fprintf(fp,"0 0 0 na na %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d \n",  
                    genotype_ori.fitness,
                    sqrt(genotype_ori.sq_SE_fitness),
                    genotype_ori.avg_GR1,
                    genotype_ori.avg_GR2,
                    sqrt(genotype_ori.sq_SE_GR1),
                    sqrt(genotype_ori.sq_SE_GR2),
                    genotype_ori.ngenes,
                    genotype_ori.nproteins,
                    genotype_ori.N_act,
                    genotype_ori.N_rep);
            fclose(fp); 
        }
        #if RUN_FULL_SIMULATION    
        run_simulation( &genotype_ori, 
                        &genotype_ori_copy,                      
                        init_protein_number,
                        init_mRNA,                        
                        &mut_record, 
                        0, // this is the number of total mutations that have been tried
                        1, // this tells the program from which step the simulation begins
                        RS_main,
                        RS_parallel);
    
        #else // no selection at all, just randomly shuffle network topology and kinetic parameters            
        int max_mut_steps=MAX_MUT_STEP;
        evolve_neutrally(   &genotype_ori,
                            &genotype_ori_copy,                          
                            &mut_record,
                            max_mut_steps,                           
                            RS_main);     
        #endif
    #endif
    }
    /*release memory*/
    release_memory(&genotype_ori,&genotype_ori_copy,&RS_main, RS_parallel);
    return 1;	
}

void release_memory(Genotype *genotype_ori,Genotype *genotype_ori_copy,RngStream *RS_main,RngStream RS_parallel[N_THREADS])
{
    int i;
    
    for(i=0;i<NGENES;i++)
    {
        free(genotype_ori->all_binding_sites[i]);
        free(genotype_ori_copy->all_binding_sites[i]);
    }
    RngStream_DeleteStream (RS_main);
    for(i=0;i<N_THREADS;i++)
        RngStream_DeleteStream (&RS_parallel[i]);    
}

void output_genotype(Genotype *genotype, int step_i)
{
    FILE *OUTPUT;
//    char filename[32];
//    int j,k;
    
    OUTPUT=fopen(output_file,"a+");
    fprintf(OUTPUT,"%.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d \n",  
            genotype->fitness,
            sqrt(genotype->sq_SE_fitness),
            genotype->avg_GR1,
            genotype->avg_GR2,
            sqrt(genotype->sq_SE_GR1),
            sqrt(genotype->sq_SE_GR2),
            genotype->ngenes,
            genotype->nproteins,
            genotype->N_act,
            genotype->N_rep);
    fclose(OUTPUT);    
//    OUTPUT=fopen(file_act_BS,"a+");	     
//    fprintf(OUTPUT,"%d %d %d %d %d %d %d %d %d\n",
//            step_i,
//            genotype->N_act_BS[2],
//            genotype->N_act_BS[3],
//            genotype->N_act_BS[4],
//            genotype->N_act_BS[5],
//            genotype->N_act_BS[6],
//            genotype->N_act_BS[7],
//            genotype->N_act_BS[8],
//            genotype->N_act_BS[9]);
//    fclose(OUTPUT);
//    OUTPUT=fopen(file_rep_BS,"a+");
//    fprintf(OUTPUT,"%d %d %d %d %d %d %d %d %d\n",
//            step_i,
//            genotype->N_rep_BS[2],
//            genotype->N_rep_BS[3],
//            genotype->N_rep_BS[4],
//            genotype->N_rep_BS[5],
//            genotype->N_rep_BS[6],
//            genotype->N_rep_BS[7],
//            genotype->N_rep_BS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 [8],
//            genotype->N_rep_BS[9]);
//    fclose(OUTPUT);
//    OUTPUT=fopen(file_all_BS,"a+");
//    fprintf(OUTPUT,"%d %d %d %d %d %d %d %d %d\n",
//            step_i,
//            genotype->binding_sites_num[2],
//            genotype->binding_sites_num[3],
//            genotype->binding_sites_num[4],
//            genotype->binding_sites_num[5],
//            genotype->binding_sites_num[6],
//            genotype->binding_sites_num[7],
//            genotype->binding_sites_num[8],
//            genotype->binding_sites_num[9]);
//    fclose(OUTPUT);
//    for(k=N_SIGNAL_TF;k<genotype->ngenes;k++)
//    {
//            snprintf(filename,sizeof(char)*32,"CIS_%i.txt",k);			
//            OUTPUT=fopen(filename,"a+");
//            fprintf(OUTPUT,"%d",step_i);
//            for(j=0;j<CISREG_LEN;j++)
//                    fprintf(OUTPUT,"%c",genotype->cisreg_seq[k][j]);
//            fprintf(OUTPUT,"\n");
//            fclose(OUTPUT);
//    }
//    for(k=0;k<genotype->ntfgenes;k++)
//    {
//            snprintf(filename,sizeof(char)*32,"TF_%i.txt",k);
//            OUTPUT=fopen(filename,"a+");
//            fprintf(OUTPUT,"%d",step_i);
//            for(j=0;j<TF_ELEMENT_LEN;j++)
//                    fprintf(OUTPUT,"%c",genotype->tf_binding_seq[k][j]);
//            fprintf(OUTPUT,"\n");
//            fclose(OUTPUT);
//            snprintf(filename,sizeof(char)*32,"TF_r_%i.txt",k);
//            OUTPUT=fopen(filename,"a+");
//            fprintf(OUTPUT,"%d",step_i);
//            for(j=0;j<TF_ELEMENT_LEN;j++)
//                    fprintf(OUTPUT,"%c",genotype->tf_binding_seq_rc[k][j]);
//            fprintf(OUTPUT,"\n");
//            fclose(OUTPUT);
//    }
}

void print_core_i1ffls(Genotype *genotype)
{
    FILE *fp; 
    fp=fopen("proportion_c1ffl.txt","a+");
    fprintf(fp,"%d %d %d %d %d %d %d\n",
            genotype->N_motifs[0],
            genotype->N_motifs[1],
            genotype->N_motifs[2],
            genotype->N_motifs[3],
            genotype->normalizer1,
            genotype->normalizer2,
            genotype->n_output_genes);                                                                                                
    fclose(fp); 
}

void print_all_protein_time_courses(int nprotein, TimeCourse **timecoursestart, TimeCourse **timecourselast, FILE *fp)
{
    int i,j,k;
    TimeCourse *temp;
    temp=timecoursestart[0];
    k=0;
    while(temp)
    {
        k++;
        temp=temp->next;
    }
    for(i=0;i<k;i++)
    {
        fprintf(fp,"%f  %f  ",timecoursestart[0]->time,timecoursestart[0]->concentration);
        timecoursestart[0]=timecoursestart[0]->next;
        for(j=1;j<nprotein;j++)
        {
            fprintf(fp,"%f  ",timecoursestart[j]->concentration);
            timecoursestart[j]=timecoursestart[j]->next;
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void summarize_binding_sites(Genotype *genotype,int step_i)
{
    FILE *OUTPUT1;
    int i,j,protein_id;
    int table[NGENES][NGENES];
//    int table_reduced_sites[NGENES][NGENES];
        
    for(i=0;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->ngenes;j++)
        {
            table[i][j]=0; 
//            table_reduced_sites[i][j]=0;
        }
    }
    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)        
    {
        //print_binding_sites_distribution(genotype,i,step_i);
//        resolve_overlapping_sites(genotype,i,table_reduced_sites[i]);
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {            
            table[i][genotype->all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
        }    
    }
    
    /*Output all binding sites*/ 
    OUTPUT1=fopen("summary_BS.txt","a+");
    fprintf(OUTPUT1,"step %d\n",step_i);
    fprintf(OUTPUT1,"Promoter ");    
    for(i=0;i<genotype->nproteins;i++)
    {
        if(genotype->protein_identity[i][0]==1)
        {
            if(genotype->protein_identity[i][1]!=-1)
                fprintf(OUTPUT1," A%d ",i);
            else
                fprintf(OUTPUT1," a%d ",i);
//            N_act++;
        }
        if(genotype->protein_identity[i][0]==0)
        {
            if(genotype->protein_identity[i][1]!=-1)
                fprintf(OUTPUT1," R%d ",i);
            else
                fprintf(OUTPUT1," r%d ",i);
//            N_rep++;
        }
    }
    fprintf(OUTPUT1,"MIN_ACT_TO_TURN_ON\n");
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        if(i<9)
            fprintf(OUTPUT1,"%d        ",i+1);
        else
            fprintf(OUTPUT1,"%d       ",i+1);
        
        for(j=0;j<genotype->nproteins;j++)
        {
            if(table[i][j]<10)
                fprintf(OUTPUT1," %d  ",table[i][j]);
            else
                fprintf(OUTPUT1," %d ",table[i][j]);
        }
        if(genotype->protein_identity[genotype->which_protein[i]][1]!=-1)
            fprintf(OUTPUT1,"E%d",genotype->which_protein[i]);
        else
        {
            if(genotype->protein_identity[genotype->which_protein[i]][0]==1)
                fprintf(OUTPUT1,"a%d",genotype->which_protein[i]); 
            if(genotype->protein_identity[genotype->which_protein[i]][0]==0)
                fprintf(OUTPUT1,"r%d",genotype->which_protein[i]);
        }
        fprintf(OUTPUT1," %d \n",genotype->min_N_activator_to_transc[i]);
    }
    fprintf(OUTPUT1,"\n"); 
    fprintf(OUTPUT1,"\n"); 
    fclose(OUTPUT1);
    
    /*Output non-overlapping sites*/
//    OUTPUT1=fopen("summary_non_overlapping_BS","a+");
//    fprintf(OUTPUT1,"step %d\n",step_i);
//    fprintf(OUTPUT1,"Promoter ");
//    for(i=0;i<genotype->nproteins-1;i++)
//    {
//        if(genotype->protein_identity[i]==1)
//        {
//            fprintf(OUTPUT1," A%d ",i);
////            N_act++;
//        }
//        if(genotype->protein_identity[i]==0)
//        {
//            fprintf(OUTPUT1," R%d ",i);
////            N_rep++;
//        }
//    }
//    fprintf(OUTPUT1,"\n");
//    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
//    {
//        if(i<9)
//            fprintf(OUTPUT1,"%d        ",i+1);
//        else
//            fprintf(OUTPUT1,"%d       ",i+1);
//        
//        for(j=0;j<genotype->nproteins-1;j++)
//        {
//            if(table[i][j]<10)
//                fprintf(OUTPUT1," %d  ",table[i][j]);
//            else
//                fprintf(OUTPUT1," %d ",table[i][j]);
//        }
//        if(genotype->which_protein[i]==genotype->nproteins-1)
//            fprintf(OUTPUT1,"S");
//        else
//        {
//            if(genotype->protein_identity[genotype->which_protein[i]]==1)
//                fprintf(OUTPUT1,"A%d",genotype->which_protein[i]); 
//            if(genotype->protein_identity[genotype->which_protein[i]]==0)
//                fprintf(OUTPUT1,"R%d",genotype->which_protein[i]);
//        }
//        fprintf(OUTPUT1,"\n");
//    }
//    fprintf(OUTPUT1,"\n"); 
//    fprintf(OUTPUT1,"\n"); 
//    fclose(OUTPUT1);
}

void print_binding_sites_distribution(Genotype *genotype, int which_gene, int init_or_end)
{
    FILE *fp;
    int i;
    char filename[32];
    snprintf(filename,sizeof(char)*32,"distribution_on_%i_%i.txt",which_gene,init_or_end);
    fp=fopen(filename,"w");
    fprintf(fp,"y tf_id pos AR\n");   
    for(i=0;i<genotype->binding_sites_num[which_gene];i++)
    {
        fprintf(fp,"0 %d %d %d\n",genotype->all_binding_sites[which_gene][i].tf_id,genotype->all_binding_sites[which_gene][i].BS_pos+8,genotype->protein_identity[genotype->all_binding_sites[which_gene][i].tf_id][0]);
    }
    fclose(fp);   
}

void resolve_overlapping_sites(Genotype *genotype, int which_gene, int N_non_overlapping_sites[NGENES])
{
    int i,head,tail;
    int temp[NGENES][100]; // assuming there are at most 100 BS of the same TF on a promoter.
    int temp2[NGENES]; // counter of the number of 
    
    for(i=0;i<genotype->ngenes-1;i++)
        temp2[i]=0;
    
    for(i=0;i<genotype->binding_sites_num[which_gene];i++)
    {
        temp[genotype->all_binding_sites[which_gene][i].tf_id][temp2[genotype->all_binding_sites[which_gene][i].tf_id]]=i;
        temp2[genotype->all_binding_sites[which_gene][i].tf_id]++;
    }  
    
    for(i=0;i<genotype->ngenes;i++)
    {
        if(temp2[i]>1)
        { 
            head=0;
            tail=1;
            N_non_overlapping_sites[i]=1;
            while(tail<temp2[i])
            {
                if(genotype->all_binding_sites[which_gene][temp[i][tail]].BS_pos-genotype->all_binding_sites[which_gene][temp[i][head]].BS_pos-1<TF_ELEMENT_LEN+2*HIND_LENGTH)
                    tail++;                
                else
                {
                    head=tail;
                    tail++;
                    N_non_overlapping_sites[i]++;
                }
            } 
        }
        else
            N_non_overlapping_sites[i]=temp2[i];
    }    
}

void find_i1ffl(Genotype *genotype)
{
    int i,j,k,l,cluster_size;
    int found_bs;
    int gene_id,gene_id_copy,site_id,protein_id,N_copies,N_activators, N_repressors;
    int repressors[NPROTEINS];
    int activators[NPROTEINS];
    int regulated_by_signal[NGENES];
//    int copies_reg_by_env[NGENES],copies_not_reg_by_env[NGENES],N_copies_reg_by_env,N_copies_not_reg_by_env;
    int N_act_effector_genes;
    int N_rep_non_effector_genes;
    int N_rep_effector_genes_per_protein[MAX_OUTPUT_PROTEINS],N_rep_effector_genes,N_rep_effector;
    int N_all_motifs;
    int N_motifs;    
  
    /*reset variables*/   
    genotype->normalizer1=0;
    genotype->normalizer2=0;   
    N_act_effector_genes=0;
    N_rep_non_effector_genes=0;
    N_rep_effector_genes=0;
    N_rep_effector=0;
    for(i=0;i<MAX_OUTPUT_PROTEINS;i++)
        N_rep_effector_genes_per_protein[i]=0;
    for(i=0;i<NGENES;i++)
    {
        genotype->gene_in_core_C1ffl[i]=0;
        for(j=0;j<NPROTEINS;j++)
            genotype->TF_in_core_C1ffl[i][j]=0;
    }     
    for(i=0;i<4;i++)    
        genotype->N_motifs[i]=0;  
    /*begin searching motifs*/
    if(genotype->N_rep>0)
    {
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
        {
            for(j=0;j<genotype->binding_sites_num[i];j++)
            {
                if(genotype->all_binding_sites[i][j].tf_id==N_SIGNAL_TF)                
                    regulated_by_signal[i]=1;
                else
                    regulated_by_signal[i]=0;
                break;                
            }
        }    
        
        
        i=0;   
        while(genotype->cisreg_cluster[i][0]!=-1) 
        {    
//            N_copies_reg_by_env=0;
//            N_copies_not_reg_by_env=0;
//            for(j=0;j<NGENES;j++)
//            {
//                copies_reg_by_env[j]=-1;
//                copies_not_reg_by_env[j]=-1;
//            }
           
            gene_id=genotype->cisreg_cluster[i][0];
            gene_id_copy=gene_id;
            protein_id=genotype->which_protein[gene_id];
            
            if(genotype->protein_identity[protein_id][1]!=NON_OUTPUT_PROTEIN) // is a effector gene
            {
                for(j=0;j<NPROTEINS;j++)
                {
                    repressors[j]=0;    
                    activators[j]=0;  
                }  
                /*effector genes can have the same cisreg_seq but differ in protein identity*/
                cluster_size=1; 
/*                cluster_size=0;
  *                  while(genotype->cisreg_cluster[i][cluster_size]!=-1)
   *                 cluster_size++;  
               */
                /*scan binding sites for tfs that regulate gene_id*/
                for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                {
                    protein_id=genotype->all_binding_sites[gene_id][j].tf_id;
                    if(genotype->which_protein[gene_id]!=protein_id)// do not look for self-regulation
                    {
                        if(genotype->locus_specific_TF_behavior[gene_id][protein_id]==ACTIVATOR ) // is a binding site of an activator                     
                            activators[protein_id]=1;
                        else
                            repressors[protein_id]=1;
                    }
                }
                /* move non-zeros entries in activators and repressors to the front. */
                k=0;
                j=0;
                N_activators=0;
                while(j<genotype->nproteins)
                {
                    if(activators[j]!=0)               
                    {
                        activators[k]=j;
                        if(k!=j)
                            activators[j]=0;
                        k++; 
                        N_activators++;
                    }    
                    j++;
                } 
                k=0;
                j=0;
                N_repressors=0;
                while(j<genotype->nproteins)
                {
                    if(repressors[j]!=0)               
                    {
                        repressors[k]=j;
                        if(k!=j)
                            repressors[j]=0;
                        k++; 
                        N_repressors++;
                    }    
                    j++;
                }            
                /*make lists of gene copies that are regulated by environmental signal and of those that are not*/             
//                for(j=0;j<N_repressors;j++)
//                {
//                    if(repressors[j]>=N_SIGNAL_TF)
//                    {  
//                        N_copies=genotype->protein_pool[repressors[j]][0][0]; 
//                        for(k=0;k<N_copies;k++)
//                        {
//                            gene_id=genotype->protein_pool[repressors[j]][1][k];
//                            found_bs=0;
//                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
//                            {
//                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1 )
//                                {
//                                    found_bs=1;
//                                    break;
//                                }
//                            } 
//                            if(found_bs)
//                            {
//                                copies_reg_by_env[N_copies_reg_by_env]=gene_id;
//                                N_copies_reg_by_env++;                                
//                            }
//                            else
//                            {
//                                copies_not_reg_by_env[N_copies_not_reg_by_env]=gene_id;
//                                N_copies_not_reg_by_env++;                                 
//                            }
//                        }                        
//                    }
//                } 
                /***count I1-ffls and NFB ***/  
                if(activators[0]==N_SIGNAL_TF-1) //effector is regulated by the signal
                {
                    protein_id=genotype->which_protein[gene_id_copy];
                    /*separate I1FFL, NFB, and I1FFL+NFB*/
                    for(j=0;j<N_repressors;j++)
                    {
                        for(k=0;k<genotype->protein_pool[repressors[j]][0][0];k++)
                        {
                            gene_id=genotype->protein_pool[repressors[j]][1][k];
                            if(regulated_by_signal[gene_id]==1) // gene_is regualted by the signal
                            {
                                found_bs=0;
                                for(l=0;l<genotype->binding_sites_num[gene_id];l++)
                                {
                                    if(genotype->all_binding_sites[gene_id][l].tf_id==protein_id) // gene_id is regulated by the effector
                                    {
                                        found_bs=1;
                                        break;
                                    }                                       
                                }
                                if(found_bs)
                                {                                    
                                    if(genotype->locus_specific_TF_behavior[gene_id][protein_id]==ACTIVATOR)
                                    {                                        
                                        if(genotype->locus_specific_TF_behavior[gene_id][repressors[j]]==ACTIVATOR)
                                        {
                                            found_bs=0;
                                            for(l=0;l<genotype->binding_sites_num[gene_id];l++)
                                            {
                                                if(genotype->all_binding_sites[gene_id][l].tf_id==protein_id) // is gene_id is self-regulated?
                                                {
                                                    found_bs=1;
                                                    break;
                                                }                                       
                                            }
                                            if(found_bs)
                                                genotype->N_motifs[1]++; //an I1-FFL+NBF+auto-activation
                                            else
                                                genotype->N_motifs[2]++; //an I1-FFL+NBF
                                        }
                                    }
                                }
                                else  // gene_id is not regulated by the effector
                                    genotype->N_motifs[0]++; // an I1FFL 
                            }
                            else // gene_id is not regulated by the effector
                            {
                                found_bs=0;
                                for(l=0;l<genotype->binding_sites_num[gene_id];l++)
                                {
                                    if(genotype->all_binding_sites[gene_id][l].tf_id==protein_id) // gene_id is regulated by the effector?
                                    {
                                        found_bs=1;
                                        break;
                                    }                                       
                                }
                                if(found_bs)
                                {                                    
                                    if(genotype->locus_specific_TF_behavior[gene_id][protein_id]==ACTIVATOR)
                                    {                                        
                                        if(genotype->locus_specific_TF_behavior[gene_id][repressors[j]]==ACTIVATOR)
                                        {
                                            found_bs=0;
                                            for(l=0;l<genotype->binding_sites_num[gene_id];l++)
                                            {
                                                if(genotype->all_binding_sites[gene_id][l].tf_id==protein_id) // gene_id is self-regulated? 
                                                {
                                                    found_bs=1;
                                                    break;
                                                }                                       
                                            }
                                            if(found_bs)
                                                genotype->N_motifs[3]++; //an NBF+auto-activation
                                            else
                                                genotype->N_motifs[4]++; //an NBF
                                        }
                                    }
                                }
                                //if the repressor is not regulated by the signal nor the effector, it is probably useless
                            }
                        }                            
                    }
                    
                    
                    
//                    if(genotype->protein_identity[protein_id][0]==1) // the effector is an activator
//                    {
//                        for(j=0;j<N_copies_not_reg_by_env;j++)
//                        {
//                            found_bs=0;
//                            gene_id=copies_not_reg_by_env[j];
//                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
//                            {
//                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id) 
//                                {
//                                    found_bs=1;
//                                    break;
//                                }
//                            }
//                            if(found_bs)//gene_id is regulated by the effector, we are looking at pure NFB                            
//                                genotype->N_motifs[0]++;                                                  
//                        }   
//                        for(j=0;j<N_copies_reg_by_env;j++)
//                        {
//                            found_bs=0;
//                            gene_id=copies_reg_by_env[j];
//                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
//                            {
//                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id) 
//                                {
//                                    found_bs=1;
//                                    break;
//                                }
//                            }
//                            if(found_bs)//gene_id is regulated by the effector, we are looking at I1FFL+NFB                            
//                                genotype->N_motifs[2]++;                         
//                            else // gene_id is not regulated by the effector, we are looking at I1FFL                           
//                                genotype->N_motifs[1]++;
//                        }  
//                    }
//                    else // the effector is an repressor
//                    {                        
//                        for(j=0;j<N_copies_reg_by_env;j++)
//                        {
//                            found_bs=0;
//                            gene_id=copies_reg_by_env[j];
//                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
//                            {
//                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id) 
//                                {
//                                    found_bs=1;
//                                    break;
//                                }
//                            }
//                            if(!found_bs)//gene_id is not regulated by the effector, we are looking at I1FFL
//                                genotype->N_motifs[1]++; 
//                        } 
//                        found_bs=0;
//                        for(site_id=0;site_id<genotype->binding_sites_num[gene_id_copy];site_id++)
//                        {
//                            if(genotype->all_binding_sites[gene_id_copy][site_id].tf_id==protein_id)
//                            {
//                                found_bs=1;
//                                break;
//                            }
//                        }
//                        if(found_bs)
//                            genotype->N_motifs[3]++; //auto-inhibition
//                    }
                }        
            }
            i++;
        }  
     
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(genotype->protein_identity[i][0]==1)
            {
                genotype->N_act_genes+=genotype->protein_pool[i][0][0];
                if(genotype->protein_identity[i][1]!=-1)
                    N_act_effector_genes+=genotype->protein_pool[i][0][0];
            }
            else
            {
                if(genotype->protein_identity[i][1]!=-1)
                    N_rep_non_effector_genes+=genotype->protein_pool[i][0][0];
                else
                {
                    N_rep_effector_genes_per_protein[N_rep_effector]=genotype->protein_pool[i][0][0];
                    N_rep_effector_genes+=genotype->protein_pool[i][0][0];
                    N_rep_effector++;
                }
            }
        }
        
        genotype->normalizer1=N_act_effector_genes*(N_rep_non_effector_genes+N_rep_effector_genes); //normalizer for NFB
        genotype->normalizer2=N_act_effector_genes*(N_rep_non_effector_genes+N_rep_effector_genes)+ //normalizer for I1ffl
                                N_rep_non_effector_genes*N_rep_effector_genes+
                                N_rep_effector_genes*(N_rep_effector_genes-1)/2;
        for(i=0;i<N_rep_effector;i++)
            genotype->normalizer2-=N_rep_effector_genes_per_protein[i]*(N_rep_effector_genes_per_protein[i]-1)/2;    
    }       
}

void tidy_output_files(char *file_genotype_summary, char *file_mutations)
{
    int i,replay_N_steps,N_tot_mutations;
    char buffer[600];    
    FILE *fp1,*fp2;
    
    fp1=fopen("saving_point.txt","r");
    fscanf(fp1,"%d %d",&replay_N_steps,&N_tot_mutations);
    fclose(fp1);
    
    /*Basically delete the last line of the file if it is not complete*/
    fp1=fopen(file_genotype_summary,"r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps+2;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove(file_genotype_summary);
    rename("temp",file_genotype_summary);
    
    fp1=fopen(file_mutations,"r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove(file_mutations);
    rename("temp",file_mutations);
    
//    fp1=fopen("proportion_c1ffl.txt","r");
//    fp2=fopen("temp","w");
//    for(i=0;i<replay_N_steps;i++)
//    {
//        fgets(buffer,600,fp1);
//        fputs(buffer,fp2);
//    }
//    fclose(fp1);
//    fclose(fp2);
//    remove("proportion_c1ffl.txt");
//    rename("temp","proportion_c1ffl.txt");
    
#if OUTPUT_MUTANT_DETAILS    
    fp1=fopen("MUT_Detail.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<N_tot_mutations;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("MUT_Detail.txt");
    rename("temp","MUT_Detail.txt");
    
    fp1=fopen("Mut_detail_fitness.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<N_tot_mutations;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("Mut_detail_fitness.txt");
    rename("temp","Mut_detail_fitness.txt");
#endif
    
    fp1=fopen("precise_fitness.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("precise_fitness.txt");
    rename("temp","precise_fitness.txt");
}

void calc_leaping_interval(Genotype *genotype, CellState *state, float *minimal_interval, float t_unreachable,int which_gene)
{
    int protein_id,j;
    float dt;
    float N_proteins_cause_change,Kd,N_at_end_of_simulation;
    float overall_rate;
    float P_binding;
    float t_remaining;
    float N_proteins_at_now[genotype->ngenes];
    int gene_ids[NGENES]; 
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