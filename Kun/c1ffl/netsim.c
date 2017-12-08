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

#define INITIALIZATION 2
#define UPDATE_PACT_OR_PREP 1
#define KD2APP_KD 1.8e10

int MAXELEMENTS=100; 
//int min_act_to_transcribe[N_THREADS][MAX_BINDING];
const float TRANSLATION_TIME=1.0; 
const float TRANSCRIPTION_TIME=1.0;
const float PROB_ACTIVATING=0.62;
const float TRANSCRIPTINIT=6.75; 
const float MAX_NUC_DISASSEMBLY=1.7;
const float BASE_NUC_DISASSEMBLY=0.15;
const float MAX_NUC_ASSEMBLY=6.7;
const float BASE_NUC_ASSEMBLY=0.67;
const float MAX_PIC_ASSEMBLY=3.3; 
const float BASE_PIC_ASSEMBLY=0.025;
const float MEAN_PROTEIN_DECAY_RATE=-4.58;
const float SD_PROTEIN_DECAY_RATE=0.67;
const float MEAN_TM_DISASSEMBLY_RATE=2.92;
const float SD_TM_DISASSEMBLY_RATE=0.52;
const float MEAN_MRNA_DECAY_RATE=-3.43;
const float SD_MRNA_DECAY_RATE=0.63;
const float MEAN_TRANSLATION_INIT_RATE=0.94;
const float SD_TRANSLATION_INIT_RATE=0.013;
const float MIN_Kd=1.0e-9;
const float MAX_Kd=1.0e-6;
const float log_MIN_Kd=-9.0;
const float log_MAX_Kd=-6.0;
const float NS_Kd=1.0e-5;
const float MIN_INTERVAL_TO_UPDATE_Pact_or_Prep=0.5; /* min */
const float MAX_CHANGE_IN_Pact_or_Prep=0.1;

/*Mutations*/
float SUBSTITUTION = 3.5e-10; /* susbstitution rate per site per cell division. Lynch 2008*/
float INDEL=0.0;
/*the following mutations are enabled after burn-in*/
float DUPLICATION = 0.0;   
float SILENCING = 0.0;      
float MUTKINETIC = 3.5e-9; /* 1% subs in a gene (~1kb, including introns) will change kinetic rates and binding seq */        
float proportion_mut_binding_seq = 1.0;
float proportion_mut_identity = 1.0;
float proportion_mut_koff=1.0;
float proportion_mut_kdis=1.0;
float proportion_mut_mRNA_decay=1.0;
float proportion_mut_protein_decay=1.0;
float proportion_mut_translation_rate=1.0;
float proportion_mut_cooperation=0.0;
//float max_inset=3.0;
//float max_delet=3.0;
const float sd_mutation_effect=1.78; 
//float bias_in_mut=1.15; 
const float mutational_regression_rate=0.5;
float miu_PIC_disassembly=5.22;
float miu_mRNA_decay=-1.13;
float miu_protein_decay=-4.58;
float miu_translation_init=-1.36;
float miu_Kd=-11.51;

/*fitness*/
float Ne_saturate = 100000.0;
float c_transl=2.0e-6;
float bmax=1.0; 
float duration_of_burn_in_growth_rate; /* allow cells to reach (possiblly) steady growth*/
float growth_rate_scaling = 1.0; /* set default growth rate scaling factor */
int recalc_new_fitness; /*calculate the growth rate of the new genotype four more times to increase accuracy*/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   


/*initial conditions*/
int init_TF_genes=6;
int init_N_act=3;
int init_N_rep=3;
int N_TF_GENES;
int N_EFFECTOR_GENES;

/*Environments*/
float background_signal_strength = 0.0; /* number of TF0 protein */
float signal_off_strength=0.0; /* number of signal protein at basal level */
float signal_profile_matrix[N_THREADS][200][15];
float env1_signal_strength;
float env2_signal_strength;
float env1_t_development;
float env2_t_development;
float env1_t_signal_on;    
float env1_t_signal_off;     
float env2_t_signal_on;
float env2_t_signal_off;
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
    int i, j;
    /* the first N_SIGNAL_TF genes encode the sensor TFs. The concentration of a sensor TF
     * is determined by certain environmental signal*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {  
        #if RANDOM_COOPERATION_LOGIC        
            genotype->min_act_to_transc[i]=RngStream_RandInt(RS,1,2); //if one activator is sufficient to induce expression, the gene is regualted by OR gate.
        #else
            genotype->min_act_to_transc[i]=1; 
            genotype->min_act_to_transc[genotype->ngenes-1]=2;
        #endif             
        /* tf affinity */
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
        /* mRNA decay */
        genotype->mRNAdecay[i] = exp(SD_MRNA_DECAY_RATE*gasdev(RS)+MEAN_MRNA_DECAY_RATE);
        /* protein decay */
        genotype->proteindecay[i] = exp(SD_PROTEIN_DECAY_RATE*gasdev(RS)+MEAN_PROTEIN_DECAY_RATE);        
        /* translation rate */
        genotype->translation[i] = exp(SD_TRANSLATION_INIT_RATE*gasdev(RS)+MEAN_TRANSLATION_INIT_RATE);      
        /*PIC disassembly rate*/
        genotype->pic_disassembly[i]=exp(SD_TM_DISASSEMBLY_RATE*gasdev(RS)+MEAN_TM_DISASSEMBLY_RATE); 
    }          
    /* assign tf identity*/
    genotype->N_act=0;
    genotype->N_rep=0;    
    if(init_N_rep==-1 && init_N_act==-1) /*randomly generate activators and repressors*/
    {
        for(i=N_SIGNAL_TF;i<genotype->ntfgenes;i++)
        {   
            if (RngStream_RandU01(RS)<PROB_ACTIVATING) 
            {
                genotype->N_act++; 
                genotype->activating[i] = 1;
            }
            else 
            {
                genotype->N_rep++;
                genotype->activating[i]= 0;
            }
        }
    }
    else
    {
        genotype->N_act=init_N_act;
        genotype->N_rep=init_N_rep;
        for(i=N_SIGNAL_TF;i<N_SIGNAL_TF+init_N_act;i++)            
            genotype->activating[i]=1;            
        for(i=N_SIGNAL_TF+init_N_act;i<genotype->ntfgenes;i++)
            genotype->activating[i]=0;
    }
    /* parameterize sensor TF*/
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        genotype->mRNAdecay[i]=0.0; // we assume environmental signal toggles the state of sensor TF between active and inactive 
        genotype->proteindecay[i]=0.0; // the concentration of sensor TF is constant.
        genotype->translation[i]=0.0;
        genotype->pic_disassembly[i]=0.0; 
        genotype->activating[i]=1; /*make sensor TF an activator*/
        genotype->N_act++;
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
    }
#if RANDOMIZE_SIGNAL2
    #if N_SIGNAL_TF==2
        if(RngStream_RandU01(RS)<=0.5) // we assume there is a background "on" signal, which is sensor TF 0, in the network.
            genotype->activating[1]=1; // Other sensor TFs can be either activators or repressors.
        else
        {
            genotype->activating[1]=0;
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

    genotype->ngenes=1+N_SIGNAL_TF+init_TF_genes; /*including the signal genes and 1 selection gene*/
    genotype->ntfgenes=N_SIGNAL_TF+init_TF_genes; /*including the signal genes*/
    genotype->nproteins=genotype->ngenes;  /*at initialization, each protein is encoded by one copy of gene*/    
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
    initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*NGENES, genotype->ngenes, RS);  // initialize cis-reg sequence
    initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*TFGENES, genotype->ntfgenes, RS);    //initialize binding sequence of TFs    
    /* We now generate the complementary sequence of BS that are on the non-template strand.
     * The complementary sequence is used to search for BS that on the non-template strand.  
     * We also assume that all the TFs can work on both strands, but cen induce expression in one direction.*/  
    for(i=0;i< genotype->ntfgenes;i++)
    {        
        for(k=0;k<TF_ELEMENT_LEN;k++)
        {
            switch (genotype->tf_seq[i][TF_ELEMENT_LEN-k-1])
//            switch(genotype->tf_seq[i][k])
            {
                case 'a': genotype->tf_seq_rc[i][k]='t'; break;
                case 't': genotype->tf_seq_rc[i][k]='a'; break;
                case 'c': genotype->tf_seq_rc[i][k]='g'; break;
                case 'g': genotype->tf_seq_rc[i][k]='c'; break;
            }
        }        
    }      
  //  if(!DIRECT_REG)
    //    initialization_add_regulation(genotype,RS);
    initialize_genotype_fixed(genotype, RS);     
    #if !SET_BS_MANUALLY    
        calc_all_binding_sites(genotype);
    #endif
}

void initialization_add_regulation(Genotype *genotype,RngStream RS)
{
    int i, j, k, gene_id;
    char *nucleotide;
    gene_id=RngStream_RandInt(RS,1,genotype->ntfgenes-1);
//     while(genotype->activating[gene_id]!=1)
//         gene_id=RngStream_RandInt(RS,1,genotype->ntfgenes-1);
    nucleotide=&(genotype->cisreg_seq[gene_id][0]);
    j=RngStream_RandInt(RS,HIND_LENGTH,CISREG_LEN-TF_ELEMENT_LEN-HIND_LENGTH-1);
    for(i=0;i<TF_ELEMENT_LEN;i++)
    {
        nucleotide[j]=genotype->tf_seq[0][i];
        j++;
    }
    nucleotide=&(genotype->cisreg_seq[genotype->ngenes-1][0]);
    j=RngStream_RandInt(RS,HIND_LENGTH,CISREG_LEN-TF_ELEMENT_LEN-HIND_LENGTH-1);
    k=j;
    for(i=0;i<TF_ELEMENT_LEN;i++)
    {
        nucleotide[j]=genotype->tf_seq[gene_id][i];
        j++;
    }
#if ADD_2_PATHWAYS
    k=gene_id;
    gene_id=RngStream_RandInt(RS,1,genotype->ntfgenes-1);
//     while(gene_id==k || genotype->activating[gene_id]!=1)
    while(gene_id==k)
        gene_id=RngStream_RandInt(RS,1,genotype->ntfgenes-1);
    nucleotide=&(genotype->cisreg_seq[gene_id][0]);
    j=RngStream_RandInt(RS,HIND_LENGTH,CISREG_LEN-TF_ELEMENT_LEN-HIND_LENGTH-1);
    for(i=0;i<TF_ELEMENT_LEN;i++)
    {
        nucleotide[j]=genotype->tf_seq[0][i];
        j++;
    }
    nucleotide=&(genotype->cisreg_seq[genotype->ngenes-1][0]);
    j=RngStream_RandInt(RS,HIND_LENGTH,CISREG_LEN-TF_ELEMENT_LEN-HIND_LENGTH-1);
    k=j;
    for(i=0;i<TF_ELEMENT_LEN;i++)
    {
        nucleotide[j]=genotype->tf_seq[gene_id][i];
        j++;
    }    
#else
    j=RngStream_RandInt(RS,HIND_LENGTH,CISREG_LEN-TF_ELEMENT_LEN-HIND_LENGTH-1);
    while(j<k+TF_ELEMENT_LEN+2*HIND_LENGTH && j>k-2*HIND_LENGTH-TF_ELEMENT_LEN)
        j=RngStream_RandInt(RS,HIND_LENGTH,CISREG_LEN-TF_ELEMENT_LEN-HIND_LENGTH-1);
    for(i=0;i<TF_ELEMENT_LEN;i++)
    {
        nucleotide[j]=genotype->tf_seq[gene_id][i];
        j++;
    }    
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
    char *tf_seq;
    char *cis_seq;
    char *tf_seq_rc; 
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
        
        for (k=start_TF; k < genotype->nproteins-1; k++) 
        {    
            tf_seq=&(genotype->tf_seq[k][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[k][0]);            
            /*find BS on the template strand*/
            match=0;
            for (j=i; j < i+TF_ELEMENT_LEN; j++) /*calculate the number of nucleotides that match in each [i,i+TF_ELEMENT_LEN] window. The window slides by 1 each time when scanning the promoter*/
                if (cis_seq[j] == tf_seq[j-i]) match++; 
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
                if(genotype->activating[k]==1) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand.*/
            {
                match_rc=0;
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
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
                    if(genotype->activating[k]==1) genotype->N_act_BS[gene_id]++;
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
        if(genotype->activating[genotype->all_binding_sites[gene_id][i].tf_id]==1)
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
                                int activating[NPROTEINS],
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
    if(activating[BS_info[0].tf_id]==1) // if a activator binds to BS 0   
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
        switch(activating[BS_info[m].tf_id])
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
	for(j=min_act_to_transcr;j<max_N_binding_act;j++)
		temp+=ratio_matrices[pos_next_record][0][j];
	*PaNr = (float)(temp / sum);
    
    *Pno_TF=(float)(ratio_matrices[pos_next_record][0][0]/sum);
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
                        float t_to_update_Pact_or_Prep,
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
    if(t==t_to_update_Pact_or_Prep)
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
void initialize_cell(   CellState *state,                        
                        int ngenes,
                        int nproteins,
                        int protein_pool[NPROTEINS][2][NGENES], 
                        float mRNAdecay[NGENES],
                        int init_mRNA_number[NGENES],
                        float init_protein_number[NPROTEINS], 
                        RngStream RS,
                        float tdevelopment)
{
    int i, j;
    /* start fitness at 1.0 */
    state->cumulative_fitness = 0.0;     
    state->cumulative_fitness_after_burn_in = 0.0;
    /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
    state->instantaneous_fitness = 0.0;    
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
        state->active[i]= NUC_NO_PIC;       
        state->mRNA_aft_transl_delay_num[i]=init_mRNA_number[i];
        state->mRNA_under_transl_delay_num[i]=0;
        state->mRNA_under_transc_num[i]=0;
        state->last_Pact[i]=0.0;
        state->last_Prep[i]=0.0;
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
    #if JUST_PLOTTING 
        float t;
        int N_data_points; 
        t=0.0+TIME_OFFSET;
        N_data_points=(int)tdevelopment;
        for(i=0;i<N_data_points;i++)
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
        flag='o'; 
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
                    float t,
                    int UPDATE_WHAT,
                    int thread_id)
{
    int i,cluster_id,too_much_error;
    float salphc;
    int concurrent;
    float t_to_update_Pact_or_Prep,interval_to_update_Pact_or_Prep;
    
    /*Update the rates of protein synthesis and decay*/
    /*These values are usually updated in update_protein_number_and_fitness, unless at initialization or 
     *forced to update*/
    if(UPDATE_WHAT==UPDATE_PACT_OR_PREP || UPDATE_WHAT==INITIALIZATION) 
    {
        for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
        {          
            salphc = state->mRNA_aft_transl_delay_num[i] * genotype->translation[i]/(genotype->proteindecay[i]);
            state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = genotype->proteindecay[i];
            state->konvalues[i][KON_SALPHAC_INDEX] = salphc;
        }
        for (i=0;i<N_SIGNAL_TF;i++)
        {
            state->konvalues[i][KON_PROTEIN_DECAY_INDEX]=0.0;
            state->konvalues[i][KON_SALPHAC_INDEX]=0.0;
        }
    }    
    /* reset rates */
    rates->mRNAdecay=0.0;
    rates->pic_disassembly=0.0;
    rates->nuc_disassembly=0.0;
    rates->nuc_assembly=0.0;
    rates->transcript_init=0;
    rates->pic_assembly=0.0;
    rates->subtotal=0.0;    
    for(i=0;i<genotype->ngenes;i++)
    {
        rates->nuc_disassembly_rate[i]=0.0;
        rates->nuc_assembly_rate[i]=0.0;
        rates->pic_assembly_rate[i]=0.0;
        rates->pic_disassembly_rate[i]=0.0;
        rates->mRNAdecay_rate[i]=0.0;
        rates->transcript_init_rate[i]=0;
        state->Pact[i]=0.0;
        state->Prep[i]=0.0;
        state->Pact_No_rep[i]=0.0;
        state->Pno_TF[i]=0.0;
    } 
    /* update mRNA decay rates*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        rates->mRNAdecay_rate[i] = genotype->mRNAdecay[i] * (state->mRNA_aft_transl_delay_num[i] + state->mRNA_under_transl_delay_num[i]);
        rates->mRNAdecay += rates->mRNAdecay_rate[i];        
    }
    /*update probability of binding configurations that activates expression
     * and use it to update other rates*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {    
        #if NO_REGULATION //if we manually turn off the expression of non-sensor TFs
            if(genotype->activating[genotype->which_protein[i]]==-1) // we only calculate the binding configurations for effector genes
            {
                cluster_id=genotype->which_cluster[i];        
                if(genotype->cisreg_cluster[cluster_id][0]!=i)  /*if this gene does not have a unique cis-reg sequence*/
                {                
                    state->Pact[i]=state->Pact[genotype->cisreg_cluster[cluster_id][0]]; /* copy TF distribution from elsewhere*/
                    state->Prep[i]=state->Prep[genotype->cisreg_cluster[cluster_id][0]];
                    state->Pact_No_rep[i]=state->Pact_No_rep[genotype->cisreg_cluster[cluster_id][0]];
                    state->Pno_TF[i]=state->Pno_TF[genotype->cisreg_cluster[cluster_id][0]];
                }
                else /* otherwise, we need to calc the ratio*/
                {
                    if(genotype->N_act_BS[i]!=0 || genotype->N_rep_BS[i]!=0)                
                        calc_TF_dist_from_all_BS(   genotype->all_binding_sites[i],
                                                    genotype->nproteins,
                                                    genotype->max_hindered_sites[i],
                                                    genotype->binding_sites_num[i],                                                                
                                                    genotype->activating, 
                                                    genotype->max_unhindered_sites[i],
                                                    state->protein_number,
                                                    genotype->min_act_to_transc[i],
                                                    &(state->Pact[i]),
                                                    &(state->Prep[i]),
													&(state->Pact_No_rep[i]),
                                                    &(state->Pno_TF[i]));
                    else
                    {
                        state->Pact[i]=0.0;     
                        state->Prep[i]=0.0;
						state->Pact_No_rep[i] = 0.0;
                        state->Pno_TF[i]=0.0;
                    }
                }
            }
            else
            {
               state->Pact[i]=0.0; 
               state->Prep[i]=0.0;
			   state->Pact_No_rep[i] = 0.0;
               state->Pno_TF[i]=0.0;
            }
        #else
            cluster_id=genotype->which_cluster[i];        
            if(genotype->cisreg_cluster[cluster_id][0]!=i)  /*if this gene does not have a unique cis-reg sequence*/
            {                
                state->Pact[i]=state->Pact[genotype->cisreg_cluster[cluster_id][0]]; /* copy TF distribution from elsewhere*/
                state->Prep[i]=state->Prep[genotype->cisreg_cluster[cluster_id][0]];
                state->Pact_No_rep[i]=state->Pact_No_rep[genotype->cisreg_cluster[cluster_id][0]];
                state->Pno_TF[i]=state->Pno_TF[genotype->cisreg_cluster[cluster_id][0]];
            }
            else /* otherwise, we need to calc the ratio*/
            {
                if(genotype->N_act_BS[i]!=0 || genotype->N_rep_BS[i]!=0)                
                    calc_TF_dist_from_all_BS(   genotype->all_binding_sites[i],
                                                genotype->nproteins,
                                                genotype->max_hindered_sites[i],
                                                genotype->binding_sites_num[i],                                                                
                                                genotype->activating, 
                                                genotype->max_unhindered_sites[i],
                                                state->protein_number,
                                                genotype->min_act_to_transc[i],
                                                &(state->Pact[i]),
                                                &(state->Prep[i]),
												&(state->Pact_No_rep[i]),
                                                &(state->Pno_TF[i]));
                else
                {
                    state->Pact[i]=0.0;     
                    state->Prep[i]=0.0;
					state->Pact_No_rep[i] = 0.0;
                    state->Pno_TF[i]=0.0;
                }
            }
        #endif
        /* calc other rates*/
        switch (state->active[i])
        {
            case NUC_NO_PIC:
                rates->nuc_disassembly_rate[i]=state->Pact[i]*(MAX_NUC_DISASSEMBLY-BASE_NUC_DISASSEMBLY)+BASE_NUC_DISASSEMBLY;             
                rates->nuc_disassembly+=rates->nuc_disassembly_rate[i];
                rates->nuc_assembly_rate[i]=0.0;
                rates->pic_assembly_rate[i]=0.0;
                rates->pic_disassembly_rate[i]=0.0;
                break;                
            case NO_NUC_NO_PIC:
                rates->nuc_assembly_rate[i]=state->Prep[i]*(MAX_NUC_ASSEMBLY-BASE_NUC_ASSEMBLY)+BASE_NUC_ASSEMBLY;
                rates->nuc_assembly+=rates->nuc_assembly_rate[i];
                rates->pic_assembly_rate[i]=MAX_PIC_ASSEMBLY*state->Pact_No_rep[i]+BASE_PIC_ASSEMBLY*state->Pno_TF[i];
                rates->pic_assembly+=rates->pic_assembly_rate[i]; 
                rates->pic_disassembly_rate[i]=0.0;
                rates->nuc_disassembly_rate[i]=0.0;
                break;                
            case PIC_NO_NUC: 
                rates->pic_disassembly_rate[i]=genotype->pic_disassembly[i];
                rates->pic_disassembly+=rates->pic_disassembly_rate[i]; 
                rates->transcript_init_rate[i]= 1;
                rates->transcript_init+=1;
                rates->nuc_assembly_rate[i]=0.0; 
                rates->nuc_disassembly_rate[i]=0.0;
                rates->pic_assembly_rate[i]=0.0;
                break;
        }
    }
    rates->subtotal+=rates->nuc_assembly;
    rates->subtotal+=rates->pic_assembly;
    rates->subtotal+=rates->nuc_disassembly;
    rates->subtotal+=rates->mRNAdecay;
    rates->subtotal+=rates->pic_disassembly;
    rates->subtotal+=(float)rates->transcript_init*TRANSCRIPTINIT;  
    
    /*Check if Pact needs to be updated more or less often*/
    if(UPDATE_WHAT!=INITIALIZATION)
    {
        if(UPDATE_WHAT==UPDATE_PACT_OR_PREP)
            interval_to_update_Pact_or_Prep=MIN_INTERVAL_TO_UPDATE_Pact_or_Prep;
        else
        {
            too_much_error=0;        
            for(i=N_SIGNAL_TF;i<genotype->ngenes;i++) //check if Pact changes too much for at least one gene
            {
                if(state->last_Pact[i]!=0.0)
                {
                    if(fabs(state->Pact[i]-state->last_Pact[i])/state->last_Pact[i]>MAX_CHANGE_IN_Pact_or_Prep)
                    {
                        too_much_error=1;  
                        break;                        
                    }   
                }
                if(state->last_Pact[i]==0.0 && state->Pact[i]!=0.0)
                {
                    too_much_error=1;
                    break;
                }
                if(state->last_Prep[i]!=0.0)
                {
                    if(fabs(state->Prep[i]-state->last_Prep[i])/state->last_Prep[i]>MAX_CHANGE_IN_Pact_or_Prep)
                    {
                        too_much_error=1;  
                        break;                        
                    }   
                }
                if(state->last_Prep[i]==0.0 && state->Prep[i]!=0.0)
                {
                    too_much_error=1;
                    break;
                }
				if (state->last_Pact_No_rep[i] != 0.0)
				{
					if (fabs(state->Pact_No_rep[i] - state->last_Pact_No_rep[i]) / state->last_Pact_No_rep[i]>MAX_CHANGE_IN_Pact_or_Prep)
					{
						too_much_error = 1;
						break;
					}
				}
				if (state->last_Pact_No_rep[i] == 0.0 && state->Pact_No_rep[i] != 0.0)
				{
					too_much_error = 1;
					break;
				}
                if (state->last_Pno_TF[i] != 0.0)
				{
					if (fabs(state->Pno_TF[i] - state->last_Pno_TF[i]) / state->last_Pno_TF[i]>MAX_CHANGE_IN_Pact_or_Prep)
					{
						too_much_error = 1;
						break;
					}
				}
				if (state->last_Pno_TF[i] == 0.0 && state->Pno_TF[i] != 0.0)
				{
					too_much_error = 1;
					break;
				}
            }
            if(!too_much_error)
                interval_to_update_Pact_or_Prep=(t-state->last_event_t)*2.0;//the change in Pact is moderate, double the update interval     
            else        
                interval_to_update_Pact_or_Prep=(t-state->last_event_t)*0.5;//otherwise, reduce it by a half        
            if(interval_to_update_Pact_or_Prep<MIN_INTERVAL_TO_UPDATE_Pact_or_Prep)
                interval_to_update_Pact_or_Prep=MIN_INTERVAL_TO_UPDATE_Pact_or_Prep;
        }    
        /*Update the next time that Pact will be updated mandatorily*/
        t_to_update_Pact_or_Prep=t+interval_to_update_Pact_or_Prep;
        concurrent=check_concurrence(   t_to_update_Pact_or_Prep,
                                        state->mRNA_transl_init_time_end_head,
                                        state->mRNA_transcr_time_end_head,
                                        state->signal_on_head,
                                        state->signal_off_head,
                                        state->burn_in_growth_rate_head,
                                        state->t_to_update_Pact_or_Prep,
                                        state->change_signal_strength_head);
        while(concurrent)//if the time to update overlaps with existing events, add a tiny offset
        {
            t_to_update_Pact_or_Prep+=TIME_OFFSET;
            concurrent=check_concurrence(   t_to_update_Pact_or_Prep,
                                            state->mRNA_transl_init_time_end_head,
                                            state->mRNA_transcr_time_end_head,
                                            state->signal_on_head,
                                            state->signal_off_head,
                                            state->burn_in_growth_rate_head,
                                            state->t_to_update_Pact_or_Prep,
                                            state->change_signal_strength_head);        
        }
        state->t_to_update_Pact_or_Prep=t_to_update_Pact_or_Prep;
    }
    /*Keep a copy of Pact and time for comparison at next time Pact is updated*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        state->last_Pact[i]=state->Pact[i];     
        state->last_Prep[i]=state->Prep[i];
        state->last_Pact_No_rep[i] = state->Pact_No_rep[i];
        state->last_Pno_TF[i]=state->Pno_TF[i];
    }
    state->last_event_t=t;  
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
    int retval=0;
    float t1;
    float t2;
    float t3;
    float t4;
    float t5;  
    float t6;    
    float t7;
    t1 = state->mRNA_transcr_time_end_head ? state->mRNA_transcr_time_end_head->time : TIME_INFINITY;
    t2 = state->mRNA_transl_init_time_end_head ? state->mRNA_transl_init_time_end_head->time : TIME_INFINITY;
    t3 = state->signal_off_head? state->signal_off_head->time : TIME_INFINITY;
    t4 = state->signal_on_head ? state->signal_on_head->time : TIME_INFINITY;
    t5 = state->burn_in_growth_rate_head ? state->burn_in_growth_rate_head->time : TIME_INFINITY;
    t6 = state->t_to_update_Pact_or_Prep;
    t7= state->change_signal_strength_head ? state->change_signal_strength_head->time : TIME_INFINITY;
    if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5) &&(t1<=t6) && (t1<=t7))
    {
        retval = 1;	
    }
    else if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5)&&(t2<=t6)&&(t2<=t7))
    { 
        retval = 2;
    }  
    else if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5)&&(t3<=t6) &&(t3<=t7)) 
    {
        retval = 3;
    }
    else if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5)&&(t4<=t6)&&(t4<=t7)) 
    {
        retval = 4;
    }               
    else if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4)&&(t5<=t6)&&(t5<=t7))
    {
        retval = 5;
    }             
    else if((t6 <= t1) && (t6 <= t) && (t6 <= t2) && (t6 <= t3) && (t6<=t4)&&(t6<=t5)&&(t6<=t7))
    {
        retval=6;
    }
    else if((t7 <= t1) && (t7 <= t) && (t7 <= t2) && (t7 <= t3) && (t7<=t4)&&(t7<=t5)&&(t7<=t6))
    {
        retval=7;
    }
    else
    {
        retval = 0;
    }
    return retval;
}

/*
 * end transcription: update the mRNAs ready for translation initiation
 * etc. accordingly and delete the event from the queue
 */
void fixed_event_end_transcription( float *dt,
                                    float t,
                                    CellState *state,
                                    GillespieRates *rates,
                                    Genotype *genotype,
                                    int *end_state,                        
                                    int mut_step,
                                    Mutation *mut_record,
                                    char *effect_of_effector)
{
    int gene_id;
    float concurrent;
    float endtime;
    /* recompute the delta-t based on difference between now and the time of transcription end */
    *dt = state->mRNA_transcr_time_end_head->time - t;   
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, mut_step, mut_record);
    /* get the gene which is ending transcription */
    gene_id = state->mRNA_transcr_time_end_head->event_id;    
    /* increase number of mRNAs that are initializing translation*/
    (state->mRNA_under_transl_delay_num[gene_id])++;
    /* decrease the number of mRNAs undergoing transcription */
    (state->mRNA_under_transc_num[gene_id])--;
    /* delete the fixed even which has just occurred */
    delete_fixed_event_from_head(&(state->mRNA_transcr_time_end_head), &(state->mRNA_transcr_time_end_tail));   
    /*add transcription initialization event*/  
    endtime=t+*dt+TRANSLATION_TIME;
    concurrent=check_concurrence(   endtime,
                                    state->mRNA_transl_init_time_end_head,
                                    state->mRNA_transcr_time_end_head,
                                    state->signal_on_head,
                                    state->signal_off_head,
                                    state->burn_in_growth_rate_head,
                                    state->t_to_update_Pact_or_Prep,
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
                                        state->t_to_update_Pact_or_Prep,
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
                                        float t,
                                        int *end_state,                          
                                        int mut_step,
                                        Mutation *mut_record,
                                        char *effect_of_effector)
{
    int gene_id;    
    /* calc the remaining time till translation initiation ends */
    *dt = state->mRNA_transl_init_time_end_head->time - t;         
    /* update cell growth and protein concentration during dt*/
    update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, mut_step, mut_record);
    /* get identity of gene that has just finished translating */
    gene_id=state->mRNA_transl_init_time_end_head->event_id; 
    /* there is one less mRNA that is initializing translation */
    (state->mRNA_under_transl_delay_num[gene_id])--;  
    /* delete the event that just happened */
    delete_fixed_event_from_head(&(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));    
    /* there is one more mRNA that produces protein */
    (state->mRNA_aft_transl_delay_num[gene_id])++;   
    /* update protein synthesis rate*/
    state->konvalues[gene_id][KON_SALPHAC_INDEX] = state->mRNA_aft_transl_delay_num[gene_id] * genotype->translation[gene_id] / state->konvalues[gene_id][KON_PROTEIN_DECAY_INDEX];
    
    if(genotype->which_protein[gene_id]==genotype->nproteins-1)//if the mRNA encodes a non-sensor TF, there could be a huge change in TF concentration
        return 0;
    else
        return UPDATE_PACT_OR_PREP; //so update Pact and reset updating interval to MIN_INTERVAL_TO_UPDATE_PACT
}

/*
 * compute t' factor used in the integration of growth rate
 * t' is the time the effector protein increases or decreases to saturation level
 */
float calc_tprime(Genotype *genotype, CellState* state, float *number_of_selection_protein_bf_dt, float dt, float RHS) 
{
    int n_copies;
    int i;          
    n_copies=genotype->protein_pool[genotype->nproteins-1][0][0];
    float alpha_s[n_copies],c[n_copies];
    for(i=0;i<n_copies;i++)
    {
        c[i]=state->konvalues[genotype->protein_pool[genotype->nproteins-1][1][i]][KON_PROTEIN_DECAY_INDEX];
        alpha_s[i]=c[i]*state->konvalues[genotype->protein_pool[genotype->nproteins-1][1][i]][KON_SALPHAC_INDEX];
    }       
    return rtsafe(  &calc_fx_dfx,
                    n_copies,
                    RHS,
                    number_of_selection_protein_bf_dt,
                    alpha_s,
                    c,
                    0.0,
                    dt,
                    1.0/60.0); //rtsafe is in numerical.c
}

/*
 * calculate f(x)-Pp_s and f'(x),
 * f(x) is the number of effector protein molecules at time x
 */
void calc_fx_dfx(float x, int n_copies, float RHS, float *p, float *alpha_s, float *c,float *fx, float *dfx)
{
    int i;    
    *fx=0;
    *dfx=0;    
    for(i=0;i<n_copies;i++)
    {
        *fx+=(p[i]-alpha_s[i]/c[i])*exp(-c[i]*x)+alpha_s[i]/c[i];
        *dfx+=(alpha_s[i]-c[i]*p[i])*exp(-c[i]*x);
    }
    *fx-=RHS;    
}

/*
 * calculate F(delta_t)/Pp_s. F(x) is the integral of f(x) over delta_t.
 * f(x) is the number of effector protein molecules at time x
 */
float calc_integral(Genotype *genotype, CellState *state, float *protein, float dt, float Pp)
{
    int i,n_copies,gene_ids[NPROTEINS];
    float integral=0.0,ect1;    
   
    n_copies=genotype->protein_pool[genotype->nproteins-1][0][0];
    for(i=0;i<n_copies;i++)
        gene_ids[i]=genotype->protein_pool[genotype->nproteins-1][1][i];
        
    for(i=0;i<n_copies;i++)
    {
        ect1=exp(-state->konvalues[gene_ids[i]][KON_PROTEIN_DECAY_INDEX]*dt)-1;    
        integral+=1.0/Pp*(state->konvalues[gene_ids[i]][KON_SALPHAC_INDEX]*ect1/state->konvalues[gene_ids[i]][KON_PROTEIN_DECAY_INDEX]-
                protein[i]*ect1/state->konvalues[gene_ids[i]][KON_PROTEIN_DECAY_INDEX]+ 
                state->konvalues[gene_ids[i]][KON_SALPHAC_INDEX]*dt);
    }
    return integral;
}

/*
 * return the instantaneous growth rate given the current cell state and environment,
 * also return the integrated growth rate as a pointer
 */
float calc_fitness(float *integrated_fitness,
                            Genotype *genotype,
                            CellState *state,
                            float* number_of_selection_protein_bf_dt,
                            float dt,
                            char effect_of_effector,
                            int *end_state,                         
                            int mut_step,
                            Mutation *mut_record)
{
    int i;
    float instantaneous_fitness=0.0;  /* this is returned from the function */
    float total_translation_rate = 0.0;    
    float dt_prime, dt_rest;
    float dt_fitness_flip;
    float Ne_flip;
    float penalty=bmax/Ne_saturate;   
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
            total_translation_rate += genotype->translation[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation[i]*(float)state->mRNA_under_transl_delay_num[i];
    } 
#else
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
    {     
        total_translation_rate += genotype->translation[i]*(float)state->mRNA_aft_transl_delay_num[i]+
                                    0.5*genotype->translation[i]*(float)state->mRNA_under_transl_delay_num[i];
    } 
#endif  
    cost_of_expression=total_translation_rate*c_transl;

#if ROUND_UP_NEGATIVE_FITNESS
    switch (effect_of_effector)
    {
        case 'b': /* effector is beneficial!*/
            /* choose the appropriate piecewise linear integral */           
            if ( (Ne >= Ne_saturate) && (Ne_next >= Ne_saturate) )  /* Ne > Ne_saturate throughout */
            {   
                if(bmax>cost_of_expression)
                    *integrated_fitness = bmax * dt-cost_of_expression*dt;
                else
                    *integrated_fitness =0.0;
            }
            else if ((Ne <= Ne_saturate) && (Ne_next <= Ne_saturate)) /* Ne < Ne_saturate throughout */
            {                    
                if(cost_of_expression<=bmax*Ne/Ne_saturate && cost_of_expression<=bmax*Ne_next/Ne_saturate) //the cost is always higher
                {
                    *integrated_fitness = bmax* calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                                -cost_of_expression*dt;                    
                }
                else if(cost_of_expression>=bmax*Ne/Ne_saturate && cost_of_expression>=bmax*Ne_next/Ne_saturate) //the cost is always lower
                {
                    *integrated_fitness =0.0;                                                    
                }
                else 
                {                    
                    Ne_flip=cost_of_expression*Ne_saturate/bmax; //This is the amount of effector that just canceled the cost of expression
                    dt_fitness_flip=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_flip);                    
                    if(Ne<=Ne_next) // net increase in effector protein. Instantaneous fitness is zero before dt_fitness_flip.  
                    {
                        *integrated_fitness=bmax*( calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                                        -calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_fitness_flip, Ne_saturate))
                                                -cost_of_expression*(dt-dt_fitness_flip);
                    }
                    else //net decrease in effector protein. Aft dt_fitness_flip, instantaneous fitness is zero. 
                        *integrated_fitness = bmax* calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_fitness_flip, Ne_saturate)-
                                                    cost_of_expression*dt_fitness_flip;                        
                }
            }
            else if ((Ne_saturate >= Ne) && (Ne_next >= Ne_saturate)) /* Ne < Ne_sat up until t' then Ne > Ne_sat */
            {                    
                if(cost_of_expression>=bmax) //the cost is always higher
                    *integrated_fitness =0.0;
                else if(cost_of_expression<=bmax*Ne/Ne_saturate) //the cost is always lower
                {  
                    dt_prime = calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate);                
                    dt_rest = dt - dt_prime;                
                    *integrated_fitness=bmax*( calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate)                
                                                    +dt_rest)
                                            -cost_of_expression*dt;
                }
                else                        
                {
                    Ne_flip=cost_of_expression*Ne_saturate/bmax; //This is the amount of effector that just canceled the cost of expression
                    dt_fitness_flip=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_flip); 
                    dt_prime = calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate);                
                    dt_rest = dt - dt_prime;   
                    *integrated_fitness=bmax*( calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate)
                                                    -calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_fitness_flip, Ne_saturate)
                                                    +dt_rest)
                                            -cost_of_expression*(dt-dt_fitness_flip);  
                }                
            }
            else // (Ne >= Ne_saturate) && (Ne_saturate >= Ne_next) i.e. Ne > Ne_sat up until t' then Ne < Ne_sat */
            {   
                if(cost_of_expression>=bmax) //the cost is always higher
                    *integrated_fitness =0.0;
                else if(cost_of_expression<=bmax*Ne_next/Ne_saturate) //the cost is always lower
                {                 
                    dt_prime = calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate);                 
                    dt_rest = dt - dt_prime;
                    *integrated_fitness=bmax*( dt_prime                
                                                    +calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, Ne_saturate)
                                                    -calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate))
                                            -cost_of_expression*dt;
                }
                else
                {
                    Ne_flip=cost_of_expression*Ne_saturate/bmax; //This is the amount of effector that just canceled the cost of expression
                    dt_fitness_flip=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_flip); 
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate);                
                    dt_rest=dt-dt_prime; 
                    *integrated_fitness=bmax*( dt_prime                
                                                    +calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_fitness_flip, Ne_saturate)
                                                    -calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, Ne_saturate))
                                            -cost_of_expression*dt_fitness_flip;                        
                }                
            }            
            /* compute instantaneous growth rate at t */
            if (Ne_next < Ne_saturate)
                instantaneous_fitness = bmax*Ne_next/Ne_saturate;
            else
                instantaneous_fitness = bmax;
            break;
    
        case 'd': /* effector is deleterious! */
            Ne_flip=(bmax-cost_of_expression)/penalty; //This is the amount of effector that just canceled the cost of expression            
            if(Ne>Ne_next)//decrease in effector protein
            {
                if(Ne_next>=Ne_flip) //too many effector throughout
                {
                    *integrated_fitness =0.0;
                }
                else if(Ne<=Ne_flip) // not enough effector throughout
                {
                    *integrated_fitness = bmax*dt-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)
                                                    -cost_of_expression*dt;
                }
                else // aft dt_fitness_flip, fitness is positive
                {
                    dt_fitness_flip=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_flip); 
                    *integrated_fitness = bmax*(dt-dt_fitness_flip)-penalty*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)-
                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_fitness_flip, 1.0))-
                                                    cost_of_expression*(dt-dt_fitness_flip);                    
                }                    
            }
            else // increase in effector protein
            {
                if(Ne>=Ne_flip) //too many effector throughout
                {
                    *integrated_fitness =0.0;
                }   
                else if(Ne_next<=Ne_flip)// not enough effector throughout
                {
                    *integrated_fitness = bmax*dt-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)
                                                    -cost_of_expression*dt;
                }
                else //Aft dt_fitness_flip, instantaneous fitness is zero.
                {
                    dt_fitness_flip=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_flip); 
                    *integrated_fitness = bmax*dt_fitness_flip-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_fitness_flip, 1.0)-
                                                    cost_of_expression*dt_fitness_flip;
                }                
            } 
            if(Ne_next<Ne_saturate)
                instantaneous_fitness = bmax - penalty*Ne_next;
            else
                instantaneous_fitness = 0.0;
            break;            
    } 
    /* and instantaneous integrated rate */
    instantaneous_fitness -= cost_of_expression;
    /* make sure growth rates can't be negative */  
    if (instantaneous_fitness < 0.0)
        instantaneous_fitness = 0.0;
#else
    switch (effect_of_effector)
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
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
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
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
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
                    *integrated_fitness = bmax*dt-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)
                                                    -cost_of_expression*dt;
                }
                else // aft dt_prime, the benefit becomes positive
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
                    *integrated_fitness = bmax*(dt-dt_prime)-penalty*(calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)-
                                                  calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0))-
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
                    *integrated_fitness = bmax*dt-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt, 1.0)
                                                    -cost_of_expression*dt;
                }
                else //Aft dt_prime, the benefit becomes zero
                {
                    dt_prime=calc_tprime(genotype,state,number_of_selection_protein_bf_dt,dt,Ne_saturate); 
                    *integrated_fitness = bmax*dt_prime-penalty*calc_integral(genotype, state, number_of_selection_protein_bf_dt, dt_prime, 1.0)-
                                                    cost_of_expression*dt;
                }                
            } 
            if(Ne_next<Ne_saturate)
                instantaneous_fitness = bmax - penalty*Ne_next;
            else
                instantaneous_fitness = 0.0;
            break; 
#else
            *integrated_fitness=(bmax-cost_of_expression)*dt;
            instantaneous_fitness=bmax;            
#endif
    } 
    /* and instantaneous integrated rate */
    instantaneous_fitness -= cost_of_expression;    
#endif    
    /* return the instantaneous growth rate */
    return instantaneous_fitness;
}


/* 
 * update both the protein concentration and current cell size *
 * 
 */
void update_protein_number_and_fitness( Genotype *genotype,
                                    CellState *state,                                   
                                    GillespieRates *rates,
                                    float dt,                                                                      
                                    float t,                                   
                                    char effect_of_effector,
                                    int *end_state,                                   
                                    int mut_step,
                                    Mutation *mut_record)
{
    int i,j;
    float ct, ect, ect1;
    float N_effector_molecules_bf_dt[genotype->protein_pool[genotype->nproteins-1][0][0]];
    float instantaneous_fitness = 0.0;
    float integrated_fitness = 0.0;
    FILE *fperror;
  
    /* store the numbers of the effector proteins encoded by each copy of gene before updating*/   
    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++)
        N_effector_molecules_bf_dt[i]=state->gene_specific_protein_number[genotype->protein_pool[genotype->nproteins-1][1][i]];
    /* update protein numbers*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {     
        ct = state->konvalues[i][KON_PROTEIN_DECAY_INDEX]*dt;
        ect = exp(-ct);
        if (fabs(ct)<EPSILON) ect1=ct;
        else ect1 = 1-ect;      
        /* get the new protein concentration for this gene */
        state->gene_specific_protein_number[i]=ect*state->gene_specific_protein_number[i]+state->konvalues[i][KON_SALPHAC_INDEX]*ect1 ;        
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
                                        dt,                                                                                                                     
                                        effect_of_effector,
                                        end_state,                                                   
                                        mut_step,
                                        mut_record);  
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
int Gillespie_event_mRNA_decay(GillespieRates *rates, CellState *state, Genotype *genotype, float t, RngStream RS)
{
    int gene_id;
    float x;
    int mRNA_id;
    while(1)/*in case of numerical error*/    
    {           
        x=RngStream_RandU01(RS)*rates->mRNAdecay;
        gene_id=N_SIGNAL_TF-1;
        /* loop through mRNA products, to choose the mRNA with the
        proportionally higher decay rate */
        while (gene_id < genotype->ngenes-1 && x > 0.0) 
        {
            gene_id++;
            x-= rates->mRNAdecay_rate[gene_id];
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
        state->konvalues[gene_id][KON_SALPHAC_INDEX] = state->mRNA_aft_transl_delay_num[gene_id] * genotype->translation[gene_id] / state->konvalues[gene_id][KON_PROTEIN_DECAY_INDEX];;
    } 
    else 
    {
        /* decay mRNA in process of translation initialization */       
        mRNA_id = RngStream_RandInt(RS,0,state->mRNA_under_transl_delay_num[gene_id]-1);
        /* delete this fixed event: this mRNA will never be translated */
        delete_fixed_event(gene_id, mRNA_id, &(state->mRNA_transl_init_time_end_head), &(state->mRNA_transl_init_time_end_tail));       
        /* remove the mRNA from the count */
        (state->mRNA_under_transl_delay_num[gene_id])--; 
    }
    
    if(genotype->which_protein[gene_id]==genotype->nproteins-1)
        return 0;
    else // an mRNA of transcription factor is degraded, which can cause fluctuation in transcription factor concentrations.
        return UPDATE_PACT_OR_PREP;
}

void Gillespie_event_histone_acteylation(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x= RngStream_RandU01(RS)*rates->nuc_disassembly;
        gene_id=N_SIGNAL_TF-1;
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->nuc_disassembly_rate[gene_id];
        }
        if(rates->nuc_disassembly_rate[gene_id]>0.0)
            break;
    }
    /* set state: eject nucleosome, but there is no PIC yet */
    state->active[gene_id] = NO_NUC_NO_PIC;
}

void Gillespie_event_histone_deacteylation(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id; 
    float x; 
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->nuc_assembly;
        gene_id=N_SIGNAL_TF-1;
        /* choose a particular gene copy to change state */
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->nuc_assembly_rate[gene_id];
        }
        if(rates->nuc_assembly_rate[gene_id]>0.0)
            break;
    }
    /* set state: nucleosome returns */
    state->active[gene_id] = NUC_NO_PIC;
}

void Gillespie_event_assemble_PIC(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    float x; 
    int gene_id;
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->pic_assembly;
        gene_id=N_SIGNAL_TF-1;
        /* choose a particular gene copy to change state */
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->pic_assembly_rate[gene_id];
        }
        if(rates->pic_assembly_rate[gene_id]>0.0)
            break;
    }
    /* turn gene fully on: ready for transcription and adjust rates */
    state->active[gene_id] = PIC_NO_NUC;  
}

void Gillespie_event_disassemble_PIC(Genotype *genotype, CellState *state,GillespieRates *rates, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x=RngStream_RandU01(RS)*rates->pic_disassembly;   
        gene_id=N_SIGNAL_TF-1;
        /* choose an appropriate gene copy to disassemble the PIC from */
        while(gene_id < genotype->ngenes-1 && x>0.0) 
        {
            gene_id++;       
            x -= rates->pic_disassembly_rate[gene_id];
        }
        if(rates->pic_disassembly_rate[gene_id]>0.0)
            break;
    }
    state->active[gene_id]=NO_NUC_NO_PIC;   
}

void Gillespie_event_transcription_init(GillespieRates *rates, CellState *state, Genotype *genotype, float dt, float t, RngStream RS)
{
    int gene_id;  
    int x;
    float candidate_t;
    int concurrent;    
    gene_id=N_SIGNAL_TF-1;    
    x=RngStream_RandInt(RS,1,rates->transcript_init);
    while(gene_id<genotype->ngenes-1 && x>0)
    {
        gene_id++;
        x-=rates->transcript_init_rate[gene_id];
    }     
    /* now that transcription of gene has been initiated, 
     * we add the timepoint at which the transcription ends, 
     * which is dt+time-of-transcription from now */
    candidate_t=t+dt+TRANSCRIPTION_TIME;
    concurrent=check_concurrence(   candidate_t,
                                    state->mRNA_transl_init_time_end_head,
                                    state->mRNA_transcr_time_end_head,
                                    state->signal_on_head,
                                    state->signal_off_head,
                                    state->burn_in_growth_rate_head,
                                    state->t_to_update_Pact_or_Prep,
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
                                        state->t_to_update_Pact_or_Prep,                                        
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
        genotype_clone->min_act_to_transc[i]=MAX_BINDING+1;
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
    /* copy binding sites' sequences*/  
    for(i=0; i < genotype_templet->ntfgenes; i++) 
    {          
        for(j=0;j<TF_ELEMENT_LEN;j++)
        {    
            genotype_clone->tf_seq[i][j]=genotype_templet->tf_seq[i][j];
            genotype_clone->tf_seq_rc[i][j]=genotype_templet->tf_seq_rc[i][j];
        }
    }
    /*copy kinetic constants*/
    for(i=0; i < genotype_templet->ngenes; i++) 
    {            
        genotype_clone->mRNAdecay[i]=genotype_templet->mRNAdecay[i];
        genotype_clone->proteindecay[i]=genotype_templet->proteindecay[i];
        genotype_clone->translation[i]=genotype_templet->translation[i];            
        genotype_clone->pic_disassembly[i]=genotype_templet->pic_disassembly[i];
        genotype_clone->which_protein[i]=genotype_templet->which_protein[i];
        genotype_clone->min_act_to_transc[i]=genotype_templet->min_act_to_transc[i];
    } 
    /* copy TF information*/
    for(i=0;i<NPROTEINS;i++)
    {
        genotype_clone->activating[i]=genotype_templet->activating[i];
        genotype_clone->Kd[i]=genotype_templet->Kd[i];
    }    
    /* copy gene and protein numbers*/
    genotype_clone->ngenes=genotype_templet->ngenes;
    genotype_clone->ntfgenes=genotype_templet->ntfgenes;
    genotype_clone->nproteins=genotype_templet->nproteins;
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
                        float *t,                          
                        char *effect_of_effector, 
                        float signal_on_strength,
                        float duration_signal_on,
                        float duration_signal_off,
                        char init_effector_effect,
                        int fixed_effector_effect,
                        RngStream RS,                       
                        int mut_step,
                        Mutation *mut_record,
                        int *end_state,
                        int *did_burn_in,
                        int thread_id,
                        float tdevelopment,
                        float *signal_profile) 
{    
    int event, UPDATE_WHAT;     
    float fixed_time; 
    float dt;
    float x;
    FILE *fperror;

    x = expdev(RS);        /* draw random number */
    dt = x/rates->subtotal;
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
    fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;
    event = does_fixed_event_end(state, fixed_time);
    while(event!=0)
    {           
        /*after doing fixed event, return a flag to indicate whether mandatorily update Pact or Prep*/
        UPDATE_WHAT=do_fixed_event( genotype, 
                                    state, 
                                    rates, 
                                    &dt, 
                                    *t, 
                                    event, 
                                    signal_on_strength,
                                    duration_signal_on, 
                                    duration_signal_off, 
                                    effect_of_effector, 
                                    init_effector_effect,
                                    fixed_effector_effect, 
                                    end_state, 
                                    mut_step, 
                                    mut_record,
                                    signal_profile);        
        if(*end_state==0) 
            return; 
        *t += dt;  /* advance time by the dt */        
#if CAUTIOUS //watch for overlapping fixed events -- one of which will be skipped and cause error
        float dt_copy,x_copy,rate_subtotal_copy,t_copy;  
        int event2;
        if(event==5)
        {
            *did_burn_in=1;
            state->cell_size_after_burn_in=state->cell_size;
        }          
        fixed_time=*t+dt;
        event2=does_fixed_event_end(state, fixed_time);
        if(event2!=0)
        {
            fperror=fopen(error_file,"a+");
            LOG("overlapping fixed events at mut_step %d", mut_step);
            fclose(fp);
        }                       
        t_copy=*t; 
        dt_copy=dt;
        x_copy=x; 
        rate_subtotal_copy=rates->subtotal; 
#endif        
        x -= dt*rates->subtotal;  /* we've been running with rates->subtotal for dt, so substract it from x*/        
        calc_all_rates(genotype, state, rates, *t, UPDATE_WHAT,thread_id);  /* update rates->subtotal and re-compute a new dt */      
        dt = x/rates->subtotal;
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
        fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment; 
        /* check to see there aren't more fixed events to do */
        event = does_fixed_event_end(state, fixed_time);                                    
    } 
  /* no remaining fixed events to do in dt, now do stochastic events */  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (*t+dt < tdevelopment)
    {
        update_protein_number_and_fitness(genotype, state, rates, dt, *t, *effect_of_effector, end_state, mut_step, mut_record);        
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, *t, RS, end_state, mut_step, mut_record);
        if(*end_state==0)
            return;  
        /* Gillespie step: advance time to next event at dt */
        *t += dt;
        calc_all_rates(genotype,state,rates,*t,UPDATE_WHAT,thread_id);        
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - *t;
        /* final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, dt, *t, *effect_of_effector, end_state, mut_step, mut_record); 
        if(*end_state==0)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        *t = tdevelopment;
    }
}

/* while there are fixed events
     occuring in current t->dt window */
int do_fixed_event(Genotype *genotype, 
                    CellState *state, 
                    GillespieRates *rates, 
                    float *dt,
                    float t,        
                    int event,
                    float signal_on_strength,
                    float duration_signal_on,
                    float duration_signal_off,
                    char *effect_of_effector,
                    char init_effector_effect,
                    int fixed_effector_effect,                
                    int *end_state,                  
                    int mut_step,
                    Mutation *mut_record,
                    float *signal_profile)
{     
    int return_value;
    return_value=0;
    switch (event) 
    {
        case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, t, state, rates, genotype,end_state,mut_step,mut_record,effect_of_effector); 
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, dt, t,end_state,mut_step,mut_record,effect_of_effector);
            break;
        case 3:     /* turn signal off*/ 
            *dt = state->signal_off_head->time - t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, mut_step, mut_record); 
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));
            if(fixed_effector_effect)
                *effect_of_effector=init_effector_effect;
            else
                *effect_of_effector='d';
            state->protein_number[N_SIGNAL_TF-1]=signal_off_strength;
            break;
        case 4:     /*turn signal on*/
            *dt = state->signal_on_head->time - t;   
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, mut_step, mut_record);  
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
            state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;
            if(fixed_effector_effect)                               
                *effect_of_effector=init_effector_effect;            
            else                 
                *effect_of_effector='b'; 
            break;	
        case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, mut_step, mut_record);
            state->cumulative_fitness_after_burn_in=state->cumulative_fitness;           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6: /* mandatorily updating Pact and Prep*/
            *dt=state->t_to_update_Pact_or_Prep-t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, mut_step, mut_record);          
            break;
        case 7: /* update signal strength */
            *dt=state->change_signal_strength_head->time-t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, mut_step, mut_record);
            state->protein_number[N_SIGNAL_TF-1]=signal_profile[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
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
                        float t,
                        RngStream RS,
                        int *end_state,                        
                        int mut_step,
                        Mutation *mut_record)
{
    float x; 
    int return_value;
    return_value=0;
    FILE *fperrors; 
    while(1)
    {    
        x=RngStream_RandU01(RS)*(rates->subtotal);          
        if (x <= rates->mRNAdecay)  /*STOCHASTIC EVENT: an mRNA decay event */
        {    
            return_value=Gillespie_event_mRNA_decay(rates, state, genotype, t, RS);                
            break;
        } 
        else 
        {               
            x -= rates->mRNAdecay;             
            if (x <= rates->pic_disassembly) /* STOCHASTIC EVENT: PIC disassembly*/
            {                   
                Gillespie_event_disassemble_PIC(genotype, state, rates, RS);
                break;
            } 
            else 
            {  
                x -= rates->pic_disassembly;  
                if (x <= rates->nuc_disassembly)  /* acetylation*/
                {    
                    Gillespie_event_histone_acteylation(rates, state, genotype, RS);
                    break;
                } 
                else 
                { 
                    x-= rates->nuc_disassembly;                  
                    if (x <= rates->nuc_assembly)/* STOCHASTIC EVENT: histone deacetylation */ 
                    {   
                        Gillespie_event_histone_deacteylation(rates, state, genotype, RS);
                        break;
                    } 
                    else 
                    {
                        x -= rates->nuc_assembly;                        
                        if (x <= rates->pic_assembly)/* STOCHASTIC EVENT: PIC assembly*/
                        {  
                            Gillespie_event_assemble_PIC(rates, state, genotype, RS); 
                            break;
                        } 
                        else 
                        {
                            x -= rates->pic_assembly;                            
                            if (x <= rates->transcript_init * TRANSCRIPTINIT) /* STOCHASTIC EVENT: transcription initiation */
                            {   
                                Gillespie_event_transcription_init(rates, state, genotype, dt, t, RS);
                                break;
                            }                           
                            else 
                            {                               
                              /* shouldn't get here, unless rounding error*/
                                fperrors=fopen(error_file,"a+");
                                LOG("at time %f in mut step %d\n", t, mut_step)                            
                                fclose(fperrors); 
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
}

/**
 *Calculate the fintess of a given genotype.
 *Essentially calling do_single_timestep until tdevelopment and calculate 
 *average growth rate over tdevelopment.
 */
void calc_avg_growth_rate(  Genotype *genotype,
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
        float t;
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];       
        int did_burn_in;  
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
            genotype_clone.ntfgenes=genotype->ntfgenes;
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
    #if !SET_BS_MANUALLY
        calc_all_binding_sites(&genotype_clone);
    #endif        
        /* now calc growth rate under two environments*/
        for(i=0;i<N_replicates_per_thread;i++) /* env 1, usually a constant signal that matches env*/
        {	 
            effect_of_effector=env1_initial_effect_of_effector; // whether the effector is beneficial or deleterious
            end_state=1; //flag to show whether do_single_time_step encounters an error. 1 means no error, 0 means error
            did_burn_in=0; //if burn-in growth is enabled, and correctly executed, this flag will be set to 1 by do_single_time_step.
#if EXTERNAL_SIGNAL
            j=RngStream_RandInt(RS_parallel[thread_ID],0,99);
            signal_profile=&(signal_profile_matrix[thread_ID][j][0]);
#else
            signal_profile=NULL;             
#endif
            /*Set initial mRNA and protein number using given values*/
            for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)  //expression of the sensor TF is not considered in the model           
                mRNA[j] = init_mRNA_clone[j];                       
            for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
            {
                for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                    protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0]; //split the initial protein number equally to different copies
                                                                                                                                            //this is to make sure all proteins have equal initial numbers
            }
            /*initialize mRNA and protein numbers, and gene states, in cell*/
            initialize_cell(&state_clone,                             
                            genotype_clone.ngenes, 
                            genotype_clone.nproteins,
                            genotype_clone.protein_pool,
                            genotype_clone.mRNAdecay, 
                            mRNA, 
                            protein, 
                            RS_parallel[thread_ID],
                            env1_t_development);
            /*set how the environment signal should change during simulation*/
            set_signal(&state_clone,
                        env1_t_signal_on,
                        env1_t_signal_off,
                        signal_profile,
                        env1_t_development,
                        env1_signal_strength);          
            /*growth starts at 0 min*/
            t = 0.0;
            /*calcualte the rates of cellular activity based on the initial cellular state*/
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, t, INITIALIZATION,thread_ID); 
            /*run growth simulation until tdevelopment or encounter an error*/
            while(t<env1_t_development && end_state==1) 
            {
                 do_single_timestep(&genotype_clone, 
                                    &state_clone, 
                                    &rate_clone, 
                                    &t,                                                                                   
                                    &effect_of_effector,
                                    env1_signal_strength,
                                    env1_t_signal_on,
                                    env1_t_signal_off,
                                    env1_initial_effect_of_effector,
                                    env1_fixed_effector_effect,
                                    RS_parallel[thread_ID],                                   
                                    mut_step,  
                                    mut_record,
                                    &end_state,
                                    &did_burn_in,
                                    thread_ID,
                                    env1_t_development,
                                    signal_profile);
            } 
            if(end_state==1) // no error
            {	
                /*calculate average growth rate*/
                gr1[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/(env1_t_development-duration_of_burn_in_growth_rate);              
            #if CAUTIOUS                
                if(gr1[i]<0.0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("negative growth rate at replicate %d in test 1 thread %d at mut step %d\n", i, thread_ID, mut_step);
                    fclose(fp);
                }
                if(did_burn_in==0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("burn_in did not engage at replicate %d in test 1 thread %d at mut step %d\n", i, thread_ID, mut_step);
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
            free_fixedevent(&state_clone);           
        }            
        /*******************env2*********************/
        for(i=0;i<N_replicates_per_thread;i++) 
        {	 
            effect_of_effector=env2_initial_effect_of_effector;  
            end_state=1;
            did_burn_in=0;
            signal_profile=NULL;  
            for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)             
                mRNA[j] = init_mRNA_clone[j];                       
            for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
            {
                for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                    protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0];                
            }
            initialize_cell(&state_clone, 
                            genotype_clone.ngenes, 
                            genotype_clone.nproteins,
                            genotype_clone.protein_pool,
                            genotype_clone.mRNAdecay, 
                            mRNA, 
                            protein, 
                            RS_parallel[thread_ID],
                            env2_t_development);
            set_signal( &state_clone,
                        env2_t_signal_on,   
                        env2_t_signal_off,
                        signal_profile,
                        env2_t_development,
                        env2_signal_strength);  
            t = 0.0;
            calc_all_rates( &genotype_clone, 
                            &state_clone, 
                            &rate_clone, 
                            t, 
                            INITIALIZATION,thread_ID);	
            while(t<env2_t_development && end_state==1)
            {
                do_single_timestep(&genotype_clone, 
                                    &state_clone, 
                                    &rate_clone, 
                                    &t,                                                                                 
                                    &effect_of_effector,
                                    env2_signal_strength,
                                    env2_t_signal_on,
                                    env2_t_signal_off,
                                    env2_initial_effect_of_effector,
                                    env2_fixed_effector_effect,                                  
                                    RS_parallel[thread_ID],                                  
                                    mut_step,
                                    mut_record,
                                    &end_state,
                                    &did_burn_in,
                                    thread_ID,
                                    env2_t_development,
                                    signal_profile);
            } 
            if(end_state==1)
            {
                gr2[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/(env2_t_development-duration_of_burn_in_growth_rate);             
#if CAUTIOUS
                FILE *fp;
                if(gr2[i]<0.0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("negative growth rate at replicate %d in test 2 thread %d at mut step %d\n", i, thread_ID, mut_step);
                    fclose(fp);
                }
                if(did_burn_in==0)
                {
                    fperror=fopen(error_file,"a+");
                    LOG("burn_in did not engage at replicate %d in test 2 thread %d at mut step %d\n", i, thread_ID, mut_step);
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
    t7 = state->t_to_update_Pact_or_Prep;
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
                                float t,        
                                int event, 
                                float duration_signal_on,
                                float duration_signal_off,
                                char *effect_of_effector,
                                char init_effector_effect,
                                int fixed_effector_effect,                              
                                int *end_state,
                                float *signal_profile)
{      
    int return_value;
    return_value=0;
    switch (event) 
    {        
         case 1:     /* a transcription event ends */
            fixed_event_end_transcription(dt, t, state, rates, genotype,end_state,0,NULL,effect_of_effector); 
            break;
        case 2:     /* a translation initialization event ends */ 
            return_value=fixed_event_end_translation_init(genotype, state, rates, dt, t,end_state,0,NULL,effect_of_effector);
            break;
        case 3:     /* turn off signal*/ 
            *dt = state->signal_off_head->time - t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, 0,NULL); 
            delete_fixed_event_from_head(&(state->signal_off_head),&(state->signal_off_tail));            
            state->protein_number[N_SIGNAL_TF-1]=signal_off_strength;
            if(fixed_effector_effect)
                *effect_of_effector=init_effector_effect;
            else
                *effect_of_effector='d';
            state->protein_number[N_SIGNAL_TF-1]=signal_off_strength;
            break;
        case 4:     /* turn on signal*/
            *dt = state->signal_on_head->time - t;   
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, 0,NULL);  
            delete_fixed_event_from_head(&(state->signal_on_head),&(state->signal_on_tail));
//            state->protein_number[N_SIGNAL_TF-1]=signal_on_strength;
             if(fixed_effector_effect)                               
                *effect_of_effector=init_effector_effect;            
            else                 
                *effect_of_effector='b'; 
            break;	
        case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-t;     
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, 0,NULL);
            state->cumulative_fitness_after_burn_in=state->cumulative_fitness;           
            delete_fixed_event_from_head(&(state->burn_in_growth_rate_head),&(state->burn_in_growth_rate_tail));
            break;
        case 6:
            *dt=state->sampling_point_end_head->time-t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, 0, NULL);
            delete_fixed_event_from_head(&(state->sampling_point_end_head),&(state->sampling_point_end_tail));
            break;
        case 7: /* update force to update Pact and Prep*/
            *dt=state->t_to_update_Pact_or_Prep-t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, 0,NULL);          
            break;  
        case 8: /* update signal strength */
            *dt=state->change_signal_strength_head->time-t;
            update_protein_number_and_fitness(genotype, state, rates, *dt, t, *effect_of_effector, end_state, 0, NULL);
            state->protein_number[N_SIGNAL_TF-1]=signal_profile[state->change_signal_strength_head->event_id];
            delete_fixed_event_from_head(&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            break;
    }    
    return return_value;
}

void do_single_timestep_plotting(   Genotype *genotype, 
                                    CellState *state,                         
                                    GillespieRates *rates, 
                                    float *t,                         
                                    char *effect_of_effector,                                    
                                    float duration_signal_on,
                                    float duration_signal_off,
                                    int fixed_effector_effect,
                                    char initial_effect_of_effector,    
                                    float (*phenotype)[90],                                    
                                    float fitness[90],                                    
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
    dt = x/rates->subtotal;    
    fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;
    event = does_fixed_event_end_plotting(state, fixed_time);
    while(event!=0)
    {           
        UPDATE_WHAT=do_fixed_event_plotting(genotype, 
                                            state, 
                                            rates, 
                                            &dt, 
                                            *t, 
                                            event, 
                                            duration_signal_on, 
                                            duration_signal_off, 
                                            effect_of_effector,
                                            initial_effect_of_effector,
                                            fixed_effector_effect, 
                                            end_state,
                                            signal_profile);
        if(*end_state==0)
            return;        
        if(event==6)/*time to take a snapshot of fitness and protein number*/
        {    
            for(i=0;i<genotype->nproteins;i++)
                phenotype[i][*timepoint]=state->protein_number[i];            
            fitness[*timepoint]=state->instantaneous_fitness;          
            (*timepoint)++;    
        }              
        fixed_time=*t+dt; 
        *t += dt;                     
        x -= dt*rates->subtotal;       
        calc_all_rates(genotype, state, rates, *t, UPDATE_WHAT,thread_id);  
        dt = x/rates->subtotal;	
        if(dt<0.0)	
            dt=TIME_OFFSET;         
        /* check to see there aren't more fixed events to do */
        fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;         
        event = does_fixed_event_end_plotting(state, fixed_time);                                    
    } 
    /* no remaining fixed events to do in dt, now do stochastic events */  
    /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (*t+dt < tdevelopment)
    { 
        update_protein_number_and_fitness(genotype, state, rates, dt, *t, *effect_of_effector, end_state, 0, NULL);  
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, *t, RS, end_state, 0, NULL);
        if(*end_state==0)
            return;        
        calc_all_rates(genotype,state,rates,*t,UPDATE_WHAT,thread_id);       
        /* Gillespie step: advance time to next event at dt */
        *t += dt;
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - *t;
        /* final update of protein concentration */
        update_protein_number_and_fitness(genotype, state, rates, dt, *t, *effect_of_effector, end_state, 0, NULL); 
       
        if(*end_state==0)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        *t = tdevelopment;
    }
}

void calc_avg_growth_rate_plotting( Genotype *genotype,
                                    int init_mRNA[NGENES],
                                    float init_protein_number[NPROTEINS],
                                    RngStream RS[N_THREADS])
{     
    float phenotypeA[N_REPLICATES][genotype->nproteins][90];
    float phenotypeB[N_REPLICATES][genotype->nproteins][90];
    float fitnessA[N_REPLICATES][90];
    float fitnessB[N_REPLICATES][90];
    int l,m,n;
    char filename1[32],filename2[32];
    FILE *fp1,*fp2;
    
    for(l=0;l<N_REPLICATES;l++)
    {
        for(m=0;m<genotype->nproteins;m++)
        {
            for(n=0;n<90;n++)
            {
                phenotypeA[l][m][n]=0.0;
                phenotypeB[l][m][n]=0.0;               
            }
        }
    }    
    for(l=0;l<N_REPLICATES;l++)
    {   
        for(n=0;n<90;n++)
            {
                fitnessA[l][n]=0.0;
                fitnessB[l][n]=0.0;               
            }        
    }
   
//    #pragma omp parallel num_threads(N_THREADS)
    {
//        int thread_ID=omp_get_thread_num();
        int thread_ID=0;
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
//        #pragma omp critical
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
        for(i=0;i<N_replicates_per_thread;i++)        
        {             
            effect_of_effector=env1_initial_effect_of_effector;
            end_state=1;
            signal_profile=NULL;
            for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)             
                mRNA[j] = init_mRNA_clone[j];                       
            for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
            {
                for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                    protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0];                
            }
            initialize_cell(&state_clone, genotype_clone.ngenes, genotype_clone.nproteins,
                            genotype_clone.protein_pool,genotype_clone.mRNAdecay, mRNA, protein, RS[thread_ID],env1_t_development);
            set_signal(&state_clone,env1_t_signal_on,env1_t_signal_off,signal_profile,env1_t_development,env1_signal_strength);        
            t = 0.0;
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, t, INITIALIZATION,thread_ID);
            timepoint=0;
            while(t<env1_t_development && end_state==1)
            {    
                do_single_timestep_plotting(&genotype_clone, 
                                            &state_clone, 
                                            &rate_clone, 
                                            &t,                                                                                   
                                            &effect_of_effector,
                                            env1_t_signal_on,
                                            env1_t_signal_off,
                                            env1_fixed_effector_effect,  
                                            env1_initial_effect_of_effector,
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
            free_fixedevent(&state_clone);
        }
        for(i=0;i<N_replicates_per_thread;i++)        
        {              
            effect_of_effector=env2_initial_effect_of_effector; 
            end_state=1;
            signal_profile=NULL;            
            for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)             
                mRNA[j] = init_mRNA_clone[j];                       
            for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
            {
                for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                    protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0];                
            }
            initialize_cell(&state_clone, genotype_clone.ngenes, genotype_clone.nproteins,
                            genotype_clone.protein_pool,genotype_clone.mRNAdecay, mRNA, protein, RS[thread_ID],env2_t_development);
            set_signal(&state_clone,env2_t_signal_on,env2_t_signal_off,signal_profile,env2_t_development,env2_signal_strength);        
            t = 0.0;
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, t, INITIALIZATION,thread_ID);        
            timepoint=0;
            while(t<env2_t_development && end_state==1)
            {
                do_single_timestep_plotting(&genotype_clone, 
                                            &state_clone, 
                                            &rate_clone, 
                                            &t,                                                                                          
                                            &effect_of_effector,
                                            env2_t_signal_on,
                                            env2_t_signal_off,
                                            env2_fixed_effector_effect,  
                                            env2_initial_effect_of_effector,
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
            free_fixedevent(&state_clone);
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
        for(m=0;m<90;m++)
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
    for(m=0;m<90;m++)
    {
        for(n=0;n<N_REPLICATES;n++)
        {
            fprintf(fp1,"%f ",fitnessA[n][m]);
            fprintf(fp2,"%f ",fitnessB[n][m]);
        }
        fprintf(fp1,"\n");
        fprintf(fp2,"\n");
    }  
    fclose(fp1);
    fclose(fp2);
}
#endif /*End of plotting functions*/


/*
 ************ begin of mutation functions **************
 */
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
 *small (1-3 nt) insertion and deletion to cis-reg are no longer modeled,
 *but they can be turned on. HOWEVER, the code for comparing binding sites
 *were not added to these functions. 
 */
void mut_insertion(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
//    int inset_size=0;
//    int which_gene,which_nucleotide,i;
//    char n; 
//
//    inset_size=RngStream_RandInt(RS,1,max_inset);
//    which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
//    which_nucleotide=RngStream_RandInt(RS,0,CISREG_LEN-1);                /* at which new seq will be inserted*/
//    
//    if (which_nucleotide+inset_size>CISREG_LEN) inset_size=CISREG_LEN-which_nucleotide;
//
//    for(i=1;i<=CISREG_LEN-inset_size-which_nucleotide;i++)    
//        genotype->cisreg_seq[which_gene][CISREG_LEN-i]=genotype->cisreg_seq[which_gene][CISREG_LEN-inset_size-i];
//    for(i=which_nucleotide;i<which_nucleotide+inset_size;i++)
//    {        					
//        n=set_base_pair(RngStream_RandU01(RS));	
//        genotype->cisreg_seq[which_gene][i]=n;
//        /*record mutation info*/
//        mut_record->nuc_diff[i-which_nucleotide]=n;
//    }
//    
//    genotype->recalc_TFBS[which_gene]=1;  /*recalc TFBS*/
//    genotype->clone_info[which_gene]=1;   /*restore info in clone_cell*/
//    update_cisreg_cluster(genotype,which_gene,'i',NULL,-1,-1);    
//    
//    /*record mutation info*/
//    mut_record->which_nucleotide=which_nucleotide;
//    mut_record->which_gene=which_gene;    
}

void reproduce_insertion(Genotype *genotype, Mutation *mut_record)
{
//    int inset_size;    
//    int which_gene,which_nucleotide,i;    
//    
//    inset_size=(int)strlen(mut_record->nuc_diff);
//    which_gene=mut_record->which_gene;
//    which_nucleotide=mut_record->which_nucleotide;               /* at which new seq will be inserted*/
//    
//    if (which_nucleotide+inset_size>CISREG_LEN) inset_size=CISREG_LEN-which_nucleotide;
//
//    for(i=1;i<=CISREG_LEN-inset_size-which_nucleotide;i++)    
//        genotype->cisreg_seq[which_gene][CISREG_LEN-i]=genotype->cisreg_seq[which_gene][CISREG_LEN-inset_size-i];
//    for(i=which_nucleotide;i<which_nucleotide+inset_size;i++)      	
//        genotype->cisreg_seq[which_gene][i]=mut_record->nuc_diff[i-which_nucleotide];    
//    genotype->recalc_TFBS[which_gene]=1;  /*recalc TFBS*/
//    genotype->clone_info[which_gene]=1;   /*restore info in clone_cell*/
//    update_cisreg_cluster(genotype,which_gene,'i',NULL,-1,-1);   
}

void mut_partial_deletion(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
//    int delet_size=0;
//    int which_gene, which_nucleotide, i;
//    char n;    
//
//    delet_size = RngStream_RandInt(RS,1,max_delet);              
//    which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);/* from which a seq will be deleted */
//    which_nucleotide=RngStream_RandInt(RS,0,CISREG_LEN-1);
//    
//    /*record mutation info*/
//    mut_record->which_gene=which_gene;
//    mut_record->which_nucleotide=which_nucleotide;
//
//    if (which_nucleotide+delet_size>CISREG_LEN)	/* if only the tail is deleted*/
//    {
//        delet_size=CISREG_LEN-which_nucleotide;
//        
//        for(i=which_nucleotide;i<which_nucleotide+delet_size;i++)
//        {					
//            n=set_base_pair(RngStream_RandU01(RS));
//            genotype->cisreg_seq[which_gene][i]=n;
//             /*record mutation info*/
//            mut_record->nuc_diff[i-which_nucleotide]=n;
//        }				
//    }
//    else /* otherwise, join the two fragments aside the deletion */
//    {
//        for(i=which_nucleotide;i<CISREG_LEN-delet_size;i++)
//        {
//            genotype->cisreg_seq[which_gene][i]=genotype->cisreg_seq[which_gene][i+delet_size];
//        }                				
//        for(i=CISREG_LEN-delet_size;i<CISREG_LEN;i++) /* and fill the gab by generating new seq */
//        {				
//            n= set_base_pair(RngStream_RandU01(RS));
//            genotype->cisreg_seq[which_gene][i]=n;
//             /*record mutation info*/
//            mut_record->nuc_diff[i-(CISREG_LEN-delet_size)]=n;
//        }						
//    }	
//    
//    genotype->recalc_TFBS[which_gene]=1;  /*recalc TFBS*/
//    genotype->clone_info[which_gene]=1;   /*restore info in clone_cell*/
//    update_cisreg_cluster(genotype,which_gene,'p',NULL,-1,-1);  
}

void reproduce_partial_deletion(Genotype *genotype, Mutation *mut_record)
{
//    int delet_size;   
//    int which_gene, which_nucleotide, i;    
//    
//    delet_size=(int)strlen(mut_record->nuc_diff);		
//    which_gene=mut_record->which_gene;    		
//    which_nucleotide=mut_record->which_nucleotide;                /* from which a seq will be deleted */
//   
//    if (which_nucleotide+delet_size>CISREG_LEN)	/* if only the tail is deleted*/
//    {
//        delet_size=CISREG_LEN-which_nucleotide;
//        
//        for(i=which_nucleotide;i<which_nucleotide+delet_size;i++)
//        {
//            genotype->cisreg_seq[which_gene][i]=mut_record->nuc_diff[i-which_nucleotide];            
//        }				
//    }
//    else /* else, join the two fragments aside the deletion */
//    {
//        for(i=which_nucleotide;i<CISREG_LEN-delet_size;i++)
//        {
//            genotype->cisreg_seq[which_gene][i]=genotype->cisreg_seq[which_gene][i+delet_size];
//        }                				
//        for(i++;i<CISREG_LEN;i++) /* and fill the gab by generating new seq */
//        {          
//            genotype->cisreg_seq[which_gene][i]=mut_record->nuc_diff[i-which_nucleotide];           
//        }						
//    }	
//    
//    genotype->recalc_TFBS[which_gene]=1;  /*recalc TFBS*/
//    genotype->clone_info[which_gene]=1;   /*restore info in clone_cell*/
//    update_cisreg_cluster(genotype,which_gene,'p',NULL,-1,-1); 
}

/**
 *Deletion whole cis-reg sequence
 */
void mut_whole_gene_deletion(Genotype *genotype, Mutation *mut_record, RngStream RS) // any gene can be deleted
{
    int i,j;
    int which_gene, offset, cisreg_seq_cluster_id,cisreg_seq_cluster_id_copy;
    char *temp;
    int protein_id, N_gene_of_selection_protein, puppet_gene,position_of_puppet_gene;    
    /* check which genes can be deleted*/   
    N_gene_of_selection_protein=genotype->protein_pool[genotype->nproteins-1][0][0]; //does the effector gene have multipe copies   
    if(genotype->ntfgenes-N_SIGNAL_TF>1)//if there are more than one copies of genes of non-sensor TF  
    {
        if(N_gene_of_selection_protein==1)//if the effector have only one copies
        {
            protein_id=genotype->nproteins-1;
            while(protein_id==genotype->nproteins-1) //the last copy of the effector cannot be deleted. choose a different gene
            {   
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                protein_id=genotype->which_protein[which_gene];
            }
        }
        else //otherwise any gene can be deleted
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    }
    else //if there's only one copy of non-sensor TF gene
    {        
        while(1)
        {
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
            protein_id=genotype->which_protein[which_gene];
            if(protein_id==genotype->nproteins-1)//obviously we can only delete a copy of the effector gene
                break;                           //it's impossible for deletion when the non-sensor TF and the effector each has only one copy
                                                 //function draw_mutation would have set deletion rate to 0 in such a case
        }
    }
    /*record mutation info*/
    mut_record->which_gene=which_gene;    
    /*To preserve data structure, the last gene has to be an effector gene. If this gene is selected for deletion, 
     * we overwrite it with one of its copy. Then we deleted the original copy. 
     * For example, assuming ngenes=10, gene 9 is an effector gene and is marked for deletion. 
     * We know gene 5 was duplicated from gene 9, so we overwrite gene 9 with all the info of gene 5, 
     * then we delete gene 5 (the puppet gene).*/
    if(which_gene==genotype->ngenes-1)
    {
        cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
        /*check if the effector gene has a copy that has the same cis-reg sequence*/
        if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=-1)//there is at least one copy that meets the requirement
        {            
            puppet_gene=genotype->cisreg_cluster[cisreg_seq_cluster_id][0];//genes in cisreg_cluster is ordered acendingly
        }
        else// otherwise, this effector gene has a unique cis-reg sequence. We have to use a copy that has a different cis-reg sequence.
            // We also need to remove the cis-reg cluster of the effector gene to be deleted. 
        {     
            /*find such a copy -- call it puppet_gene*/
            protein_id=genotype->which_protein[which_gene];
            position_of_puppet_gene=0;
            if(genotype->protein_pool[protein_id][1][position_of_puppet_gene]==which_gene)position_of_puppet_gene++;
            puppet_gene=genotype->protein_pool[protein_id][1][position_of_puppet_gene];
            /*Add the effector gene to the cis_reg cluster of the puppet_gene, so that we can edit the cluster with standard procedure later*/
            cisreg_seq_cluster_id_copy=genotype->which_cluster[puppet_gene];            
            j=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=-1)j++;//find an empty slot in cisreg_cluster
            genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=which_gene;            
            genotype->which_cluster[which_gene]=cisreg_seq_cluster_id_copy;
            /*then delete which_gene's original cisreg cluster by shifting cisreg_cluster*/
            cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id;            
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=-1)
            {
                /*reset cluster=cisreg_seq_cluster_id*/
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=-1)
                {
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=-1;
                    j++;
                }
                /* copy cisreg_seq_cluster_id+1 to cisreg_seq_cluster_id, and update which_cluster*/                    
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]!=-1)
                {                        
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j];
                    genotype->which_cluster[genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]]--;
                    j++;
                }
                cisreg_seq_cluster_id_copy++;
            }
        }       
        /*overwrite which_gene with the info of puppet_gene*/
        for(j=0;j<CISREG_LEN;j++)        
            genotype->cisreg_seq[which_gene][j]=genotype->cisreg_seq[puppet_gene][j];
        genotype->min_act_to_transc[which_gene]=genotype->min_act_to_transc[puppet_gene]; 
        genotype->pic_disassembly[which_gene]=genotype->pic_disassembly[puppet_gene];           
        genotype->mRNAdecay[which_gene]=genotype->mRNAdecay[puppet_gene];
        genotype->proteindecay[which_gene]=genotype->proteindecay[puppet_gene];
        genotype->translation[which_gene]=genotype->translation[puppet_gene]; 
        genotype->recalc_TFBS[which_gene]=1;             
        /*now the effector gene is same as the puppet_gene, so delete the puppet_gene*/
        which_gene=puppet_gene;        
    }
    /*points the first nucleotide of the cis-reg sequence to be deleted*/
    temp = &genotype->cisreg_seq[which_gene][0];			
    offset=CISREG_LEN;
    /* shift the cisreg array to overwrite the cis_seq to be deleted */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-which_gene-1);i++) 
    {				
        *temp=*(temp+offset);         
        temp++;				
    }    
    /* if the to-be-deleted gene is a tf gene*/
    protein_id=genotype->which_protein[which_gene];   
    if(protein_id<genotype->nproteins-1)
    {  
        /*if this tf has only one copy of gene, then we'll delete the binding seq of this tf*/ 
        if(genotype->protein_pool[protein_id][0][0]==1) 
        {            
            /* shift the tf_reg array to overwrite the binding sequence to be deleted */
            temp=&genotype->tf_seq[protein_id][0];
            offset=TF_ELEMENT_LEN;
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }
            /* shift the tf_reg_rc array to overwrite the binding sequence to be deleted */    
            temp=&genotype->tf_seq_rc[protein_id][0];
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }   
        }
    } 
    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
    for(i=which_gene;i<genotype->ngenes;i++)
    {
        genotype->min_act_to_transc[i]=genotype->min_act_to_transc[i+1];
        genotype->pic_disassembly[i]=genotype->pic_disassembly[i+1];             
        genotype->mRNAdecay[i]=genotype->mRNAdecay[i+1];
        genotype->proteindecay[i]=genotype->proteindecay[i+1];
        genotype->translation[i]=genotype->translation[i+1];         
        genotype->recalc_TFBS[i]=1;        
    }
    /* if the to-be-deleted gene is a tf gene, change ntfgenes as well*/
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes--;    
    /* now change protein_pool and cisreg_cluster*/   
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);  
    genotype->ngenes--;
}

void reproduce_whole_gene_deletion(Genotype *genotype, Mutation *mut_record) // any gene can be deleted
{    
    int which_gene, offset, i,j, cisreg_seq_cluster_id, cisreg_seq_cluster_id_copy,puppet_gene,position_of_puppet_gene;
    char *temp;    
    int protein_id;       
    which_gene=mut_record->which_gene;    
    if(which_gene==genotype->ngenes-1)
    {
        cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
        if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=-1)
        {            
            puppet_gene=genotype->cisreg_cluster[cisreg_seq_cluster_id][0];
        }
        else 
        {
            protein_id=genotype->which_protein[which_gene];
            position_of_puppet_gene=0;
            if(genotype->protein_pool[protein_id][1][position_of_puppet_gene]==which_gene)position_of_puppet_gene++;
            puppet_gene=genotype->protein_pool[protein_id][1][position_of_puppet_gene];            
            cisreg_seq_cluster_id_copy=genotype->which_cluster[puppet_gene];            
            j=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=-1)j++;
            genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=which_gene;            
            genotype->which_cluster[which_gene]=cisreg_seq_cluster_id_copy;            
            cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id;            
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=-1)
            {                
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]!=-1)
                {
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=-1;
                    j++;
                }                                    
                j=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]!=-1)
                {                        
                    genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][j]=genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j];
                    genotype->which_cluster[genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][j]]--;
                    j++;
                }
                cisreg_seq_cluster_id_copy++;
            }
        }
        for(j=0;j<CISREG_LEN;j++)        
            genotype->cisreg_seq[which_gene][j]=genotype->cisreg_seq[puppet_gene][j];
        genotype->min_act_to_transc[which_gene]=genotype->min_act_to_transc[puppet_gene];  
        genotype->pic_disassembly[which_gene]=genotype->pic_disassembly[puppet_gene];           
        genotype->mRNAdecay[which_gene]=genotype->mRNAdecay[puppet_gene];
        genotype->proteindecay[which_gene]=genotype->proteindecay[puppet_gene];
        genotype->translation[which_gene]=genotype->translation[puppet_gene]; 
        genotype->recalc_TFBS[which_gene]=1;        
        which_gene=puppet_gene;        
    }    
    temp = &genotype->cisreg_seq[which_gene][0];			
    offset=CISREG_LEN;    
    for(i=0;i<CISREG_LEN*(genotype->ngenes-which_gene-1);i++) 
    {				
        *temp=*(temp+offset);         
        temp++;				
    }
    protein_id=genotype->which_protein[which_gene];   
    if(protein_id<genotype->nproteins-1)
    {      
        if(genotype->protein_pool[protein_id][0][0]==1) 
        {  
            temp=&genotype->tf_seq[protein_id][0];
            offset=TF_ELEMENT_LEN;
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }            
            temp=&genotype->tf_seq_rc[protein_id][0];
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp=*(temp+offset);
                temp++;
            }   
        }
    }
    for(i=which_gene;i<genotype->ngenes;i++)
    {
        genotype->min_act_to_transc[i]=genotype->min_act_to_transc[i+1];
        genotype->pic_disassembly[i]=genotype->pic_disassembly[i+1];            
        genotype->mRNAdecay[i]=genotype->mRNAdecay[i+1];
        genotype->proteindecay[i]=genotype->proteindecay[i+1];
        genotype->translation[i]=genotype->translation[i+1];         
        genotype->recalc_TFBS[i]=1;
    }
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes--;
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);
    genotype->ngenes--;   
}

/*
 * Duplicate a whole cis-reg sequence
 */
void mut_duplication(Genotype *genotype, Mutation *mut_record, RngStream RS) 
{
    int which_gene, which_gene_copy, i, protein_id;
    char *temp1, *temp2;    
    if(genotype->protein_pool[genotype->nproteins-1][0][0]>=N_EFFECTOR_GENES) // too many effector gene copies
    {     
        //note that it's not possible to have too many effector gene copies and too many tf gene copies at the same time
        //because that will make duplication rate 0.
        which_gene=genotype->ngenes-1;
        while(genotype->which_protein[which_gene]==genotype->nproteins-1)
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    }
    else
    {
        if(genotype->ntfgenes>=N_TF_GENES) // too many tf gene copies
        {
            while(1)
            {
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                if(genotype->which_protein[which_gene]==genotype->nproteins-1)
                    break;
            }    
        }            
        else //any gene can be duplicated
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1); 
    }   
    /*record mutation info*/
    mut_record->which_gene=which_gene;  
    /*To preserve data structure, the last gene has to be an effector gene.*/
    /*We want to shift the effector gene, and put the duplicated gene to the second last place.
     *But if it is the effector gene that gets duplicated, things are a little different */
    if(which_gene==genotype->ngenes-1)
        which_gene_copy=which_gene+1; 
    else
        which_gene_copy=which_gene;
    /* copy the promoter*/
    temp1=&genotype->cisreg_seq[which_gene_copy][0]; /* points to the gene to be duplicated*/
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][CISREG_LEN-1]; /* points to the end of the effector gene */
    /* shift the sequences of the effector gene CISREG_LEN bp */
    for(i=0;i<CISREG_LEN;i++) 
    {
        *(temp2+CISREG_LEN)=*temp2;
        temp2--;
    }  
    /* point temp2 to the start of the original position of the selection gene */
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0];
    /* put the duplicated gene at the original place of the selection gene */
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;
    /* shift the effector genes to make a slot*/
    genotype->min_act_to_transc[genotype->ngenes]=genotype->min_act_to_transc[genotype->ngenes-1];
    genotype->pic_disassembly[genotype->ngenes]=genotype->pic_disassembly[genotype->ngenes-1]; 
    genotype->mRNAdecay[genotype->ngenes]=genotype->mRNAdecay[genotype->ngenes-1];
    genotype->proteindecay[genotype->ngenes]=genotype->proteindecay[genotype->ngenes-1];
    genotype->translation[genotype->ngenes]=genotype->translation[genotype->ngenes-1];
    genotype->recalc_TFBS[genotype->ngenes]=1;        
    /* copy and paste info to the slot*/
    genotype->min_act_to_transc[genotype->ngenes-1]=genotype->min_act_to_transc[which_gene_copy];
    genotype->pic_disassembly[genotype->ngenes-1]=genotype->pic_disassembly[which_gene_copy];
    genotype->mRNAdecay[genotype->ngenes-1]=genotype->mRNAdecay[which_gene_copy];
    genotype->proteindecay[genotype->ngenes-1]=genotype->proteindecay[which_gene_copy];
    genotype->translation[genotype->ngenes-1]=genotype->translation[which_gene_copy];
    genotype->recalc_TFBS[genotype->ngenes-1]=1;
    /* update protein_pool*/
    protein_id=genotype->which_protein[which_gene];    
    update_protein_pool(genotype,protein_id,which_gene,'d'); 
    /* update cisreg_cluster*/    
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    /* update gene numbers*/  
    if(protein_id<genotype->nproteins-1)//note duplication do not change nproteins           
        genotype->ntfgenes++;     
    genotype->ngenes++;
}

void reproduce_gene_duplication(Genotype *genotype, Mutation *mut_record) //any gene can be duplicated
{   
    int which_gene,which_gene_copy, i, protein_id;
    char *temp1, *temp2; 
    /*get the gene to be duplicated from record*/
    which_gene=mut_record->which_gene;    
    if(which_gene==genotype->ngenes-1)
        which_gene_copy=which_gene+1; 
    else
        which_gene_copy=which_gene; 
    temp1=&genotype->cisreg_seq[which_gene_copy][0]; 
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][CISREG_LEN-1];    
    for(i=0;i<CISREG_LEN;i++) 
    {
        *(temp2+CISREG_LEN)=*temp2; 
        temp2--;
    }    
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0];  
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;
    genotype->min_act_to_transc[genotype->ngenes]=genotype->min_act_to_transc[genotype->ngenes-1];
    genotype->pic_disassembly[genotype->ngenes]=genotype->pic_disassembly[genotype->ngenes-1]; 
    genotype->mRNAdecay[genotype->ngenes]=genotype->mRNAdecay[genotype->ngenes-1];
    genotype->proteindecay[genotype->ngenes]=genotype->proteindecay[genotype->ngenes-1];
    genotype->translation[genotype->ngenes]=genotype->translation[genotype->ngenes-1];
    genotype->recalc_TFBS[genotype->ngenes]=1;
    genotype->min_act_to_transc[genotype->ngenes-1]=genotype->min_act_to_transc[which_gene_copy];
    genotype->pic_disassembly[genotype->ngenes-1]=genotype->pic_disassembly[which_gene_copy];
    genotype->mRNAdecay[genotype->ngenes-1]=genotype->mRNAdecay[which_gene_copy];
    genotype->proteindecay[genotype->ngenes-1]=genotype->proteindecay[which_gene_copy];
    genotype->translation[genotype->ngenes-1]=genotype->translation[which_gene_copy]; 
    genotype->recalc_TFBS[genotype->ngenes-1]=1;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'d');     
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes++;      
    genotype->ngenes++;     
}

/*
 * Mutation to the binding sequence of a TF gene
 */
void mut_binding_sequence(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene, which_nucleotide, protein_id, i;
    char nucleotide;      
    char *tf_seq, *tf_seq_rc,*temp1,*temp2;
    /*get which gene to mutate*/
    which_gene=RngStream_RandInt(RS,0,genotype->ngenes-2); // ngenes-1 is the selection gene, so tf genes must be 0 to ngenes-2
    while(genotype->which_protein[which_gene]==genotype->nproteins-1)
        which_gene=RngStream_RandInt(RS,0,genotype->ngenes-2);
    /*if this TF has more than one copies of gene, then the mutation adds a new tf
     *which requires a new slot in tf_seq and tf_seq_rc to store the new binding seq*/
    protein_id=genotype->which_protein[which_gene];    
    if(genotype->protein_pool[protein_id][0][0]>1)
    {
        /*points to the first nucleotide of the binding sequence*/
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        /*points to an empty slot*/
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        /*copy the binding sequences to empty slots*/
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }
        /*point tf_seq and tf_seq_rc to the new slot so that we can apply mutation later*/
        tf_seq=&genotype->tf_seq[genotype->nproteins-1][0];
        tf_seq_rc=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        /*update protein pool*/
        update_protein_pool(genotype,protein_id,which_gene,'c');  
    }    
    else /*if this tf has only one copy of gene, no new slot is required*/
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
    }
    /*mutate the binding sequence*/
    /*mutation only changes one nucleotide in the binding sequence*/
    which_nucleotide=RngStream_RandInt(RS,0,TF_ELEMENT_LEN-1);
    nucleotide=set_base_pair(RngStream_RandU01(RS));        
    while (nucleotide == tf_seq[which_nucleotide])
        nucleotide=set_base_pair(RngStream_RandU01(RS));
    tf_seq[which_nucleotide]=nucleotide;    
    /*record mutation info*/
    mut_record->which_gene=which_gene;
    mut_record->which_nucleotide=which_nucleotide;
    mut_record->nuc_diff[0]=nucleotide;    
    /* update the reverse complement sequence*/
    switch (nucleotide)
    {
        case 'g':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
    }  
    /* The binding sites on every promoter needs recalculation */
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=1;
    /*decide whether to update cisreg clusters. Mutation to binding seq may differ bs distributions among genes in a cluster*/       
    int new_clusters[NGENES][NGENES]; // 
    int genes_in_cluster[NGENES];
    int N_genes_in_cluster,no_difference,reference_gene,gene_to_be_sorted;
    int N_new_clusters,N_genes_in_new_cluster,j,k;
    calc_all_binding_sites(genotype); 
    i=N_SIGNAL_TF;
    while(genotype->cisreg_cluster[i][0]!=-1)/*check each cisreg cluster*/
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
        /*the while loop below sort genes in a cluster into groups based on whether they have the same BS distributions*/ 
        /*We use one gene in the original cluster as a reference, and sort all genes 
         *that have the same binding sites as the reference gene into a new cluster. 
         *After all genes are sorted, we check if the original cluster turns into multiple
         *new cluster.*/
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
                if(no_difference)//if the gene under comparison has the same binding sites as the reference gene
                {
                    /*put the gene into the new cluster*/
                    new_clusters[N_new_clusters][N_genes_in_new_cluster]=genes_in_cluster[which_gene];
                    N_genes_in_new_cluster++;
                    /*shift to remove the gene from the original cluster*/
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
    char *tf_seq, *tf_seq_rc, *temp1,*temp2;     
    which_gene=mut_record->which_gene;
    protein_id=genotype->which_protein[which_gene];  
    if(genotype->protein_pool[protein_id][0][0]>1)
    {        
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }      
        tf_seq=&genotype->tf_seq[genotype->nproteins-1][0];
        tf_seq_rc=&genotype->tf_seq_rc[genotype->nproteins-1][0];        
        update_protein_pool(genotype,protein_id,which_gene,'c');  
    }    
    else 
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
    }        			
    which_nucleotide=mut_record->which_nucleotide;        
    tf_seq[which_nucleotide]=mut_record->nuc_diff[1];
    switch (tf_seq[which_nucleotide])
    {
//        case 'g':
//            tf_seq_rc[which_nucleotide]='c'; break;
//        case 'c':
//            tf_seq_rc[which_nucleotide]='g'; break;
//        case 'a':
//            tf_seq_rc[which_nucleotide]='t'; break;
//        case 't':
//            tf_seq_rc[which_nucleotide]='a'; break;
        case 'g':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
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
    float random1,temp;    
    int which_gene;
    float total_mut_flux=   proportion_mut_kdis+
                            proportion_mut_mRNA_decay+
                            proportion_mut_protein_decay+
                            proportion_mut_translation_rate+
                            proportion_mut_cooperation;
    
    random1=RngStream_RandU01(RS);
    which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);    
    /*record mutation info*/
    mut_record->which_gene=which_gene;
    /********************************** mut kdis *******************************/
    if(random1<=proportion_mut_kdis/total_mut_flux) 
    {    
    #if REGRESSIVE_MUTATION
        genotype->pic_disassembly[which_gene]=mut_make_new_value(genotype->pic_disassembly[which_gene],miu_PIC_disassembly,RS);
    #else
        temp=genotype->pic_disassembly[which_gene]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);  //otherwise apply a fold change
        while(temp<MIN_PICDISASSEMBLE || temp>MAX_PICDISASSEMBLE || temp==genotype->pic_disassembly[which_gene])
            temp=genotype->pic_disassembly[which_gene]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);
        genotype->pic_disassembly[which_gene]=temp;
    #endif        
        /*record mutation info*/
        mut_record->kinetic_type=0;
        mut_record->kinetic_diff=genotype->pic_disassembly[which_gene];
    }
    else /******************************** mut mRNAdecay **********************/
    {
        random1-=proportion_mut_kdis/total_mut_flux;
        if(random1<=proportion_mut_mRNA_decay/total_mut_flux) 
        {
        #if REGRESSIVE_MUTATION
            genotype->mRNAdecay[which_gene]=mut_make_new_value(genotype->mRNAdecay[which_gene],miu_mRNA_decay,RS);
        #else            
            temp=genotype->mRNAdecay[which_gene]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut); 
            while(temp>MAX_mRNA_decay || temp<MIN_mRNA_decay || temp==genotype->mRNAdecay[which_gene])
                temp=genotype->mRNAdecay[which_gene]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);
            genotype->mRNAdecay[which_gene]=temp;
        #endif
            /*record mutation info*/
            mut_record->kinetic_type=1;
            mut_record->kinetic_diff=genotype->mRNAdecay[which_gene];
        }
        else /*************************** mut translation *********************/
        {
            random1-=proportion_mut_mRNA_decay/total_mut_flux;
            if(random1<=proportion_mut_translation_rate/total_mut_flux) 
            {
        #if REGRESSIVE_MUTATION
                genotype->translation[which_gene]=mut_make_new_value(genotype->translation[which_gene],miu_translation_init,RS);
        #else
                temp=genotype->translation[which_gene]*exp(sd_mut_effect*gasdev(RS)-bias_in_mut);
                while(temp>MAX_TRANSLATION || temp< MIN_TRANSLATION || temp==genotype->translation[which_gene])
                    temp=genotype->translation[which_gene]*exp(sd_mut_effect*gasdev(RS)-bias_in_mut);
                genotype->translation[which_gene]=temp;
        #endif
                mut_record->kinetic_type=2;
                mut_record->kinetic_diff=genotype->translation[which_gene];
            }
            else /********************* mut protein decay **********************/
            {  
                random1-=proportion_mut_translation_rate/total_mut_flux;
                if(random1<=proportion_mut_protein_decay/total_mut_flux)
                {
                    #if REGRESSIVE_MUTATION
                        genotype->proteindecay[which_gene]=mut_make_new_value(genotype->proteindecay[which_gene],miu_protein_decay,RS);
                    #else
                        temp=genotype->proteindecay[which_gene]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);  
                        while(temp>MAX_protein_decay || temp<MIN_protein_decay || temp==genotype->proteindecay[which_gene])
                            temp=genotype->proteindecay[which_gene]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut); 
                        genotype->proteindecay[which_gene]=temp;
                    #endif
                        mut_record->kinetic_type=3;
                        mut_record->kinetic_diff=genotype->proteindecay[which_gene];
                }
                else /****************** mut cooperation ***********************/
                {
                    while(genotype->which_protein[which_gene]!=genotype->nproteins-1)
                        which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                    mut_record->which_gene=which_gene;
                    if(genotype->min_act_to_transc[which_gene]==1)
                        genotype->min_act_to_transc[which_gene]=2;
                    else
                        genotype->min_act_to_transc[which_gene]=1;
                    mut_record->kinetic_type=4;
                    mut_record->kinetic_diff=(float)genotype->min_act_to_transc[which_gene];
                    update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
                }
            } 
        }
    }
}

float mut_make_new_value(float old_val, float miu, RngStream RS)
{
    float new_val;
    new_val=old_val;
    while(new_val==old_val)
    {
        new_val=old_val*exp(sd_mutation_effect*gasdev(RS)+mutational_regression_rate*(miu-log(old_val)));
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
            genotype->pic_disassembly[which_gene]=mut_record->kinetic_diff;
            break;
        case 1: /* mut mRNAdecay */ 
            genotype->mRNAdecay[which_gene]=mut_record->kinetic_diff;
            break;
        case 2: /* mut translation */                      
            genotype->translation[which_gene]=mut_record->kinetic_diff;           
            break;        
        case 3: /* mut protein decay */
            genotype->proteindecay[which_gene]=mut_record->kinetic_diff;            
            break; 
        case 4:
            genotype->min_act_to_transc[which_gene]=(int)mut_record->kinetic_diff;
            update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
}

void mut_identity(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int tf_id,protein_id,i;   
    char *tf_seq,*tf_seq_rc,*temp1,*temp2; 
    /*which tf gene to mutate*/
    tf_id = RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-2); // the first sensor tf must be an activator, therefore is not subject to mutation
    protein_id=genotype->which_protein[tf_id];    
    while(protein_id==genotype->nproteins-1) // copies of effector gene are not subjected to this mutation
    {
        tf_id = RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-2); 
        protein_id=genotype->which_protein[tf_id];
    }    
    /*save record*/
    mut_record->which_gene=tf_id;       
    /* if this tf gene has more than one copies, the mutation adds a new protein*/
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        /*give the new protein its own binding sequence*/
        tf_seq=&genotype->tf_seq[protein_id][0]; 
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        } 
        update_protein_pool(genotype,protein_id,tf_id,'e');
        /* update_protein_pool put the new protein at nproteins-1 and then increases nproteins by 1, 
        * so the new protein is at nproteins-2 now. Note that N_act and N_rep is updated in update_protein_pool*/ 
        genotype->activating[genotype->nproteins-2]=(genotype->activating[genotype->nproteins-2]==1)?0:1;         
    }
    else //otherwise we just flip the property of an exisiting TF
    {
        genotype->activating[protein_id]=(genotype->activating[protein_id]==1)?0:1; 
        if(genotype->activating[protein_id]==1)//mutate to activator
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else //otherwise to repressor
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
    }    
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=1; /* recalculate binding sites on every promoter */        
}

void reproduce_mut_identity(Genotype *genotype, Mutation *mut_record)
{
    int tf_id, protein_id,i;
    char *tf_seq,*tf_seq_rc,*temp1,*temp2;    
    tf_id = mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }     
        update_protein_pool(genotype,protein_id,tf_id,'e');  
        genotype->activating[genotype->nproteins-2]=(genotype->activating[genotype->nproteins-2]==1)?0:1;     
    }
    else
    {
        genotype->activating[protein_id]=(genotype->activating[protein_id]==1)?0:1; 
        if(genotype->activating[protein_id]==1)
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
    }    
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=1;   
}

/*
 * mutate affinity of TF
 */
void mut_koff(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int i;   
    float new_Kd; 
    int tf_id,protein_id;
    char *tf_seq,*tf_seq_rc,*temp1,*temp2;
    /*which TF to mutate*/
    tf_id=RngStream_RandInt(RS,0,genotype->ngenes-2);
    protein_id=genotype->which_protein[tf_id];
    while(protein_id==genotype->nproteins-1)
    {
        tf_id=RngStream_RandInt(RS,0,genotype->ngenes-2);
        protein_id=genotype->which_protein[tf_id];
    }
    /*generate a new koff */ 
#if REGRESSIVE_MUTATION    
    new_Kd=mut_make_new_value(genotype->Kd[protein_id],miu_Kd,RS);  
    while(new_Kd>NS_Kd)
        new_Kd=mut_make_new_value(genotype->Kd[protein_id],miu_Kd,RS);
#else 
    new_Kd=genotype->Kd[protein_id]*exp(gasdev(RS)*sd_mut_effect+bias_in_mut);
    while(new_Kd<MIN_Kd || new_Kd>MAX_Kd || new_Kd==genotype->Kd[protein_id])
        new_Kd=genotype->Kd[protein_id]*exp(gasdev(RS)*sd_mut_effect+bias_in_mut);
#endif    
    /* if this tf gene has more than one copies, the mutation adds a new protein*/    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }
        update_protein_pool(genotype,protein_id,tf_id,'f');
        /* update_protein_pool put the new protein at nproteins-1 and then increases nproteins by 1, 
         * so the new protein is at nproteins-2 now.*/
        genotype->Kd[genotype->nproteins-2]=new_Kd;  
    }                                                                                                           
    else
        genotype->Kd[protein_id]=new_Kd; 
#if SET_BS_MANUALLY 
    int j;
    for(i=1;i<genotype->ngenes-1;i++)
    {
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {
            if(genotype->all_binding_sites[i][j].tf_id==tf_id)
                genotype->all_binding_sites[i][j].Kd=genotype->Kd[protein_id]*
                        pow(NS_Kd/genotype->Kd[protein_id],(float)genotype->all_binding_sites[i][j].mis_match/(TF_ELEMENT_LEN-NMIN+1));
            
        }
    }
#endif    
    /*record mutation*/
    mut_record->which_gene=tf_id;
    mut_record->kinetic_diff=new_Kd; 
    /* recalculate binding sites on every promoter */
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=1;
}

void reproduce_mut_koff(Genotype *genotype, Mutation *mut_record)
{   
    int tf_id,protein_id,i;
    char *tf_seq, *tf_seq_rc, *temp1, *temp2;    
    tf_id=mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }    
        update_protein_pool(genotype,protein_id,tf_id,'f');  
        genotype->Kd[genotype->nproteins-2]=mut_record->kinetic_diff;  
    }    
    else
        genotype->Kd[protein_id]=mut_record->kinetic_diff;
#if SET_BS_MANUALLY 
    int j;
    for(i=1;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {
            if(genotype->all_binding_sites[i][j].tf_id==mut_record->which_gene)
                genotype->all_binding_sites[i][j].Kd=genotype->Kd[genotype->which_protein[mut_record->which_gene]]*
                        pow(NS_Kd/genotype->Kd[protein_id],(float)genotype->all_binding_sites[i][j].mis_match/(TF_ELEMENT_LEN-NMIN+1));
        }
    }
#endif
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=1;   
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
        case 'i': // insertion in cis-reg       
            mut_insertion(genotype,mut_record,RS);
            break;
        case 'p': // small deletion in cis-reg       
            mut_partial_deletion(genotype,mut_record,RS);
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
            mut_koff(genotype,mut_record,RS);            
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
        case 'i': // insertion        
            reproduce_insertion(genotype,mut_record);
            break;
        case 'p': // partial deletion        
            reproduce_partial_deletion(genotype,mut_record);
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
            reproduce_mut_koff(genotype,mut_record);
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
    float tot_subs_rate, tot_indel_rate, tot_dup_rate, tot_sil_rate, tot_mut_kin_rate, tot_mut_identity_rate, tot_mut_binding_seq_rate, tot_mut_koff_rate; 
    int N_genes_to_be_duplicated; 
    /* duplication rate*/    
    tot_dup_rate=0.0;
    N_genes_to_be_duplicated=genotype->ngenes-N_SIGNAL_TF;//NA to sensor TF gene
    if(genotype->ntfgenes>=N_TF_GENES)//too many non-sensor TF genes.
        N_genes_to_be_duplicated-=genotype->ntfgenes-N_SIGNAL_TF; //do not duplicate non-sensor TF gene anymore
    if(genotype->protein_pool[genotype->nproteins-1][0][0]>=N_EFFECTOR_GENES)//too many effector gene
        N_genes_to_be_duplicated-=genotype->protein_pool[genotype->nproteins-1][0][0];//do not duplicate effector gene anymore
    tot_dup_rate=N_genes_to_be_duplicated*DUPLICATION;
    tot_mut_rate+=tot_dup_rate;   
    /* silencing rate*/ 
    if(genotype->ntfgenes-N_SIGNAL_TF>1)//if there are more than one copy of non-sensor TF gene
        tot_sil_rate=(genotype->ntfgenes-N_SIGNAL_TF)*SILENCING;    
    else
        tot_sil_rate=0.0; //the last non-sensor TF cannot be deleted
    if(genotype->protein_pool[genotype->nproteins-1][0][0]>1)//if there are more than one effector gene
        tot_sil_rate+=genotype->protein_pool[genotype->nproteins-1][0][0]*SILENCING;
    else
        tot_sil_rate+=0.0; // otherwise last copy effector gene cannot be deleted either*/
    tot_mut_rate+=tot_sil_rate;    
    /* calc total susbtitution rate*/
    tot_subs_rate=(genotype->ngenes-N_SIGNAL_TF)*CISREG_LEN*SUBSTITUTION; //NA to the sensor TF gene
    tot_mut_rate+=tot_subs_rate;    
    /* indel rate*/
    tot_indel_rate=(genotype->ngenes-N_SIGNAL_TF)*CISREG_LEN*INDEL; //NA to the sensor TF gene
    tot_mut_rate+=tot_indel_rate;    
    /* mut in kinetic constants */    
    tot_mut_kin_rate=(genotype->ngenes-N_SIGNAL_TF)*MUTKINETIC*(proportion_mut_kdis+proportion_mut_mRNA_decay+proportion_mut_protein_decay+proportion_mut_translation_rate); // NA to the sensor TF gene
    tot_mut_kin_rate+=genotype->protein_pool[genotype->nproteins-1][0][0]*MUTKINETIC*proportion_mut_cooperation;
    tot_mut_rate+=tot_mut_kin_rate;
    /* mut in binding seq*/
    tot_mut_binding_seq_rate=genotype->ntfgenes*MUTKINETIC*proportion_mut_binding_seq; // NA to the effector genes
    tot_mut_rate+=tot_mut_binding_seq_rate;    
    /* mut in identity*/
    tot_mut_identity_rate=(genotype->ntfgenes-N_SIGNAL_TF)*MUTKINETIC*proportion_mut_identity; // NA to the sensor TF gene
    tot_mut_rate+=tot_mut_identity_rate;
    /* mut in koff*/
    tot_mut_koff_rate=genotype->ntfgenes*MUTKINETIC*proportion_mut_koff;  
    tot_mut_rate+=tot_mut_koff_rate;    
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
                if(random<=tot_indel_rate/tot_mut_rate)
                {
                    if(RngStream_RandU01(RS)<=0.5)
                        *mut_type='i';                      /* insertion*/
                    else
                        *mut_type='p';                      /* partial deletion*/
                }
                else
                {
                    random-=tot_indel_rate/tot_mut_rate;
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
    int i, j, protein_id, gene_id; 
    /*Protein_pool stores the numbers of gene copies that encode a given protein, and the ids of these gene copies.
     *One important thing is that the genes encoding a given protein are not stored by the order of their ids in protein_pool.
     *To delete a gene, which might be the only gene, encoding a given protein, we shift protein_pool to overwrite the to-be-deleted gene
     *We need to update the ids of the remaining genes and proteins
     *For gene duplication, the new gene is always add to the end of the list of genes encoding a given protein.
     *A new protein is also add to the end of protein_pool
     *which_protein can be updated easily. Changing which protein a gene encodes is always easy. For deletion, 
     *we just shift the array to overwrite the to-be-deleted gene and update the ids of the remaining genes.*/ 
    switch (mut_type)
    {
        case 'w':/*a whole gene deletion*/        
            if(genotype->protein_pool[which_protein][0][0]==1) /* if this is the only gene copy,we also need to delete a protein*/
            {   
                /*
                 * UPDATE protein_pool for protein>=which_protein
                 */
                protein_id=which_protein;
                /*shift protein>which_protein to overwrite the to-be-deleted protein*/
                for(i=0;i<genotype->nproteins-which_protein;i++)  
                {   
                    /*reset the portion of protein_pool to be overwritten*/
                    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                        genotype->protein_pool[protein_id][1][j]=-1;
                    /*overwrite*/
                    genotype->protein_pool[protein_id][0][0]=genotype->protein_pool[protein_id+1][0][0]; //this is number of gene copies encoding a protein
                    for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
                    {
                        gene_id=genotype->protein_pool[protein_id+1][1][j];//these are the gene copies encoding a protein
                        /*note that deletion changes the ids of the remaining genes!!! Any gene that is greater than which_gene is reduced by one*/
                        genotype->protein_pool[protein_id][1][j]=(gene_id>which_gene)?gene_id-1:gene_id;
                    }            
                    protein_id++;
                }
                /*
                 * SPECIAL CASE
                 */
                /*if a tf PROTEIN is deleted, we need to change the number of TF proteins*/
                if(which_protein<genotype->nproteins-1) 
                {   
                    /* reduce the number of activators or that of repressors */ 
                    if(genotype->activating[which_protein]) 
                        genotype->N_act--;
                    else
                        genotype->N_rep--;
                    /* also remove it from activating and koff */
                    protein_id=which_protein;
                    for(i=0;i<genotype->nproteins-which_protein-1;i++)
                    {
                        genotype->activating[protein_id]=genotype->activating[protein_id+1];
                        genotype->Kd[protein_id]=genotype->Kd[protein_id+1];
                        protein_id++;
                    }                   
                    /* in the case, all genes need to recalc binding sites*/
                    for(i=N_SIGNAL_TF;i<which_gene;i++)                    
                        genotype->recalc_TFBS[i]=1; /* recalc BS */                                       
                }  
                /*
                 * UPDATE which_protein
                 */
                /* update which_protein for gene<which_gene in which_protein*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i]<which_protein)?genotype->which_protein[i]:genotype->which_protein[i]-1;//the deletion also changes the ids of proteins
                /* shift and update which_protein for gene>=which_gene in which_protein*/                
                for(i=which_gene;i<genotype->ngenes;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i+1]>which_protein)?genotype->which_protein[i+1]-1:genotype->which_protein[i+1];                   
                /*one less protein*/
                genotype->nproteins--;
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
                for(;i<genotype->protein_pool[which_protein][0][0];i++)
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
                for(i=which_gene;i<genotype->ngenes;i++)
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
            genotype->protein_pool[which_protein][1][genotype->protein_pool[which_protein][0][0]]=genotype->ngenes-1; //the newly duplicated gene takes the original place of the effector gene
            genotype->protein_pool[which_protein][0][0]++; 
            /*update the id of the original effector gene stored in protein_pool*/           
            j=0;
            while(genotype->protein_pool[genotype->nproteins-1][1][j]!=genotype->ngenes-1)j++;
            genotype->protein_pool[genotype->nproteins-1][1][j]++; 
            /*update which_protein*/          
            genotype->which_protein[genotype->ngenes]=genotype->nproteins-1;//points to the effector protein       
            genotype->which_protein[genotype->ngenes-1]=which_protein;            
            break;
        case 'c': /*mutation in tf binding seq, creating a new tf*/
            /* remove this copy of gene from the original protein_pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            /*shift to delete which_gene*/
            for(;i<genotype->protein_pool[which_protein][0][0];i++) 
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            /*one less gene copy to encoding which_protein*/
            genotype->protein_pool[which_protein][0][0]--;
            /* increase the protein id of the effector protein to make room for the new tf*/ 
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j]; //shift the protein pool of effector gene
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;//update which_protein
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1; //reset the original protein pool of the effector 
            }                       
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins-1; //put the new protein to the original place of the effector protein
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=which_gene;
            /* update activator or repressor numbers, and activating*/
            if(genotype->activating[which_protein]) //mutation to binding seq does not change the identity of a tf
                genotype->N_act++;
            else
                genotype->N_rep++;
            genotype->activating[genotype->nproteins-1]=genotype->activating[which_protein];
            /* update Kd*/
            genotype->Kd[genotype->nproteins-1]=genotype->Kd[which_protein];
            /* finally, update protein numbers*/
            genotype->nproteins++;            
            break;            
        case 'e': /*mutation in the identity of a TF, creating a new tf*/
            /* remove this copy of gene from the original protein*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0];i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            }
            genotype->protein_pool[which_protein][0][0]--;
            /* increase the protein id of the effector protein*/               
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j];
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1;
            }            
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins-1; 
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=which_gene;
            /* update activating*/
            if(genotype->activating[which_protein]) 
                genotype->N_rep++;  /* an activator turns into a repressor */
            else
                genotype->N_act++;
            genotype->activating[genotype->nproteins-1]=genotype->activating[which_protein];
            /* update Kd*/
            genotype->Kd[genotype->nproteins-1]=genotype->Kd[which_protein];
            /* finally, update protein numbers*/
            genotype->nproteins++;           
            break;            
        case 'f': /*mutation in tf koff, creating a new tf*/
            /* remove this copy of gene from the original protein pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0];i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; /* rearrange data array */
            }
            genotype->protein_pool[which_protein][0][0]--;
            /* increase the protein id of the effector protein*/               
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j];
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1;
            }            
            /* create a new protein and link it to this gene*/
            genotype->which_protein[which_gene]=genotype->nproteins-1; 
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=which_gene;            
            /* update activating*/
            if(genotype->activating[which_protein]) 
                genotype->N_act++;  
            else
                genotype->N_rep++;
            genotype->activating[genotype->nproteins-1]=genotype->activating[which_protein];  
            /* finally, update protein numbers*/
            genotype->nproteins++;
            /* NOTE: this mutation does not change the number of genes*/
            break;
    }
}


/*To reduce the amount of calculation on the probability of binding distributions, we group gene copies
 *that are created by whole gene duplication. We call such a group a cis-reg cluster because gene copies 
 *in the group should have the same cis-reg sequence. For each cis-reg cluster we only need to calculate
 *the probability of binding distributions once. However, substitutions is cis-reg sequence can create/remove
 *binding sites, therefore we need to check whether a gene copy is still in the original cis-reg cluster 
 *after mutation.We use cisreg_cluster and which_cluster to track the bi-way relation between a gene and 
 *a cis-reg cluster.*/
void update_cisreg_cluster(Genotype *genotype, int which_gene, char mut_type, int new_clusters[NGENES][NGENES], int N_new_clusters, int original_cluster_id)
{
    /*In a cis-reg cluster, gene copies are ordered ascendingly by their ids. There are no empty slots in the list 
     *of gene copies. Empty slots after the list are marked by -1. We do not track the number of gene copies in a
     *cluster. In order to count the number of gene copies correctly, marking the empty slots accurately become important.
     *Therefore, we always set all the slots in a cluster to -1, when deleting or overwrting the cluster. 
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
        if(mut_type!='w') /*not a gene deletion. Substitution (or small indel) just kicked one gene copy out of a cluster*/
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
                for(i=which_gene;i<genotype->ngenes;i++)                              
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
                for(i=which_gene;i<genotype->ngenes;i++)                                
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
        /*Assuming the last gene is in cis-reg cluster X, check whether the duplicated
         *gene is also in cluster X*/
        if(genotype->which_cluster[which_gene]!=genotype->which_cluster[genotype->ngenes-1]) //No
        {    
            /*find an empty slot in cisreg_seq_cluster_id*/
            i=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=-1)i++;
            /*the duplicated gene is always add to ngene-1. Note that ngenes has not been 
             *updated when update_cisreg_cluster is called.*/
            genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->ngenes-1;            
            /* now find the cis-reg cluster of the effector gene at ngenes-1*/           
            cisreg_seq_cluster_id_copy=genotype->which_cluster[genotype->ngenes-1];
            /*find this effector gene in the cluster*/
            i=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=genotype->ngenes-1)i++;
            /*increase the gene id by 1*/
            genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=genotype->ngenes;
        }
        else //the last gene is duplicated
        {   
            /*simply add ngenes to the cis-reg cluster of the last gene*/
            i=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=-1)i++; //find an empty slot
            genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->ngenes; //add ngenes to the slot
            genotype->cisreg_cluster[cisreg_seq_cluster_id][i-1]=genotype->ngenes-1; //this is the duplicated gene                     
        }
        /*update which_cluster*/
        genotype->which_cluster[genotype->ngenes]=genotype->which_cluster[genotype->ngenes-1]; 
        genotype->which_cluster[genotype->ngenes-1]=cisreg_seq_cluster_id;
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
        genotype->min_act_to_transc[j]=MAX_BINDING+1; /*by default a gene cannot be turned on. 
                                                       *MAX_BINDING is the maximum number of tf that 
                                                       *can bind to a cis-reg sequence.*/
        genotype->Kd[j]=-1.0;
        for(k=0;k<NGENES;k++)        
            genotype->cisreg_cluster[j][k]=-1;
    }    
    for(j=0;j<NGENES;j++)
        genotype->cisreg_cluster[NGENES][j]=-1;
    /* initialize variables that applies to protein */
    for(j=0;j<NPROTEINS;j++)
    {
        genotype->protein_pool[j][0][0]=0;
        for(k=0;k<NGENES;k++)        
            genotype->protein_pool[j][1][k]=-1;        
        genotype->activating[j]=-1;
    }
    /* alloc space for binding sites*/
    genotype->N_allocated_elements=MAXELEMENTS;
    for(j=0;j<NGENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(MAXELEMENTS*sizeof(AllTFBindingSites));
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
    float score;
#if RULE_OF_REPLACEMENT==0
    float ref,p;
    score=(mutant->fitness-resident->fitness)/sqrt(resident->sq_SE_fitness+mutant->sq_SE_fitness);
    p=pnorm(score);
    ref=RngStream_RandU01(RS);
    if(ref<p)
        *fixation=1;
    else
        *fixation=0;
#elif RULE_OF_REPLACEMENT==1
    float p,mean_rank_resident,mean_rank_mutant;    
    p=Wilcoxon_test(resident, mutant, N_measurement_resident, N_measurement_mutant,&score,&mean_rank_resident,&mean_rank_mutant);    
    if(p<ALPHA)
    {
        if(mean_rank_resident < mean_rank_mutant)
            *fixation=1;
        else
            *fixation=0;   
    }
    else
        *fixation=0;   
#elif RULE_OF_REPLACEMENT==2
    float ref;
    if(resident->fitness < mutant->fitness)
        *fixation=1;
    else
    {
        if(resident->fitness==mutant->fitness)
        {
            ref=RngStream_RandU01(RS);
            if(ref>0.5)
                *fixation=1;
            else
                *fixation=0;
        } 
        else
            *fixation=0;
    }
    score=-1.0;
#elif RULE_OF_REPLACEMENT==3
    float sd,margin,ref;
    sd=sqrt(resident->sq_SE_fitness+mutant->sq_SE_fitness);
    margin=qnorm7(1-ALPHA/2.0,sd);
    if(resident->fitness+margin < mutant->fitness)
        *fixation=1;
    else
    {
        if(resident->fitness+margin==mutant->fitness)
        {
            ref=RngStream_RandU01(RS);
            if(ref>0.5)
                *fixation=1;
            else
                *fixation=0;
        } 
        else
            *fixation=0;
    }
    score=-1.0;
#else //RULE 4: s>minimal_selection_coefficient
    float s, ref;
    s=mutant->fitness/resident->fitness-1.0;
    if(s>=minimal_selection_coefficient)
        *fixation=1;
    else          
        *fixation=0;    
    score=-1.0;
#endif
    return score;
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
    
    remove("proportion_c1ffl.txt"); 
    remove("summary_BS.txt");
    
    for(i=0;i<replay_N_steps;i++)
    {
        calc_all_binding_sites(genotype_ori);
        if(i%OUTPUT_INTERVAL==0)
            summarize_binding_sites(genotype_ori,i); 
        find_ffl(genotype_ori);
        print_core_c1ffls(genotype_ori);
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
    summarize_binding_sites(genotype_ori,i); 
    find_ffl(genotype_ori);
    print_core_c1ffls(genotype_ori);
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
        while(1)
        {
            clone_genotype(genotype_ori,genotype_ori_copy);
            mutate(genotype_ori_copy,RS_main,mut_record);         
            if(RngStream_RandU01(RS_main)<=0.5)
            {
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
        clone_genotype(genotype_ori_copy,genotype_ori);   
        i++;
        calc_all_binding_sites(genotype_ori);
        find_ffl(genotype_ori); 
        print_core_c1ffls(genotype_ori);        
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
         calc_all_binding_sites(genotype_ori);  
    #endif
    summarize_binding_sites(genotype_ori,1);   
//    exit(0);
    /* conditions under which the phenotype and fitness is measured */    
    env1_t_development=89.9;
    env2_t_development=89.9;
    env1_signal_strength=10000.0;
    env2_signal_strength=10000.0;
    duration_of_burn_in_growth_rate = 0.0; 
    env1_t_signal_on=200.0;     
    env1_t_signal_off=0.0;
    env2_t_signal_on=10.0;
    env2_t_signal_off=130.0;
    env1_initial_effect_of_effector='b';    
    env2_initial_effect_of_effector='d'; 
    env1_fixed_effector_effect=0;
    env2_fixed_effector_effect=1;
    recalc_new_fitness=1;
    env1_occurence=0.5;
    env2_occurence=0.5;            
    calc_avg_growth_rate_plotting(  genotype_ori, 
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
                calc_avg_growth_rate(   genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,                                        
                                        0,
                                        GR1[j],
                                        GR2[j],
                                        mut_record); 
            }

            calc_fitness_stats( genotype_ori,
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
        env1_t_development=89.9;
        env2_t_development=89.9;                 // global variable
        duration_of_burn_in_growth_rate=0.0;// global variable
        env1_signal_strength=10000.0;
        env2_signal_strength=10000.0;
        env1_t_signal_on=200.0;    
        env1_t_signal_off=0.0;     
        env2_t_signal_on=10.0;
        env2_t_signal_off=130.0;
        env1_initial_effect_of_effector='b';
        env2_initial_effect_of_effector='d';
        env1_fixed_effector_effect=0;    
        env2_fixed_effector_effect=1;            // global variable
        recalc_new_fitness=5;               // global variable, make sure its value is smaller than MAX_RECALC_FITNESS
        env1_occurence=0.67;                 // global variable
        env2_occurence=0.33;                 // global variable
        DUPLICATION=5.25e-9;                 
        SILENCING =5.25e-9;
        N_EFFECTOR_GENES=2;
        N_TF_GENES=9; 
        miu_PIC_disassembly=0.62;
        miu_mRNA_decay=-3.43;
        miu_protein_decay=-4.58;
        miu_translation_init=0.94;
      
        fp=fopen(RuntimeSumm,"a+");
        fprintf(fp,"**********Burn-in_I conditions**********\n");
        fprintf(fp,"BURN_IN_I=%d\n",BURN_IN_I);                
        fprintf(fp,"N_REPLICATES=%d\n",N_REPLICATES);        
        fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
        fprintf(fp,"env1_t_development=%f, env2_t_development=%f\n",env1_t_development,env2_t_development);
        fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);        
        fprintf(fp,"env1: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env1_t_signal_on, env1_t_signal_off, env1_initial_effect_of_effector, env1_fixed_effector_effect, env1_occurence);
        fprintf(fp,"env2: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env2_t_signal_on, env2_t_signal_off, env2_initial_effect_of_effector, env2_fixed_effector_effect, env2_occurence);
        fprintf(fp,"env1 init effecto effect %c, fixed? %d. env2 init effecto effect %c, fixed? %d.\n",env1_initial_effect_of_effector,env1_fixed_effector_effect,env2_initial_effect_of_effector,env2_fixed_effector_effect);
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
//         if(first_step==BURN_IN_I+1)
//             return;        
        /*2nd stage of burn-in*/
        if(BURN_IN_II)
        {
            end_state=0;
            run_burn_in=1;      
            max_mut_steps=BURN_IN_II;    
            env1_t_development=89.9;
            env2_t_development=89.9;                 // global variable
            duration_of_burn_in_growth_rate=0.0;// global variable
            env1_signal_strength=10000.0;
            env2_signal_strength=0.0;
            env1_t_signal_on=200.0;    
            env1_t_signal_off=0.0;     
            env2_t_signal_on=100.0;
            env2_t_signal_off=130.0;
            env1_initial_effect_of_effector='b';
            env2_initial_effect_of_effector='d';
            env1_fixed_effector_effect=0;    
            env2_fixed_effector_effect=1;            // global variable
            recalc_new_fitness=5;               // global variable, make sure its value is smaller than MAX_RECALC_FITNESS
            env1_occurence=0.33;                 // global variable
            env2_occurence=0.67;                 // global variable
            DUPLICATION=5.25e-9;                 
            SILENCING =5.25e-9;
            N_EFFECTOR_GENES=2;
            N_TF_GENES=9;

            fp=fopen(RuntimeSumm,"a+");
            fprintf(fp,"**********Burn-in_II conditions**********\n");
            fprintf(fp,"BURN_IN_II=%d\n",BURN_IN_II);                
            fprintf(fp,"N_REPLICATES=%d\n",N_REPLICATES);        
            fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
            fprintf(fp,"env1_t_development=%f, env2_t_development=%f\n",env1_t_development,env2_t_development);
            fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);        
            fprintf(fp,"env1: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env1_t_signal_on, env1_t_signal_off, env1_initial_effect_of_effector, env1_fixed_effector_effect, env1_occurence);
            fprintf(fp,"env2: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env2_t_signal_on, env2_t_signal_off, env2_initial_effect_of_effector, env2_fixed_effector_effect, env2_occurence);
            fprintf(fp,"env1 init effecto effect %c, fixed? %d. env2 init effecto effect %c, fixed? %d.\n",env1_initial_effect_of_effector,env1_fixed_effector_effect,env2_initial_effect_of_effector,env2_fixed_effector_effect);
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
            fprintf(fp,"Burn_in completes after the %dth step.\n",first_step);
            fclose(fp);
            
            if(end_state==-1)
                return;           
        }        
    }    
    
    /* post-burn-in simulations*/
    run_burn_in=0;
    max_mut_steps=MAX_MUT_STEP;    
    env1_t_development=89.9;
    env2_t_development=89.9;                     // global variable
    duration_of_burn_in_growth_rate=0.0;    // global variable    
    env1_signal_strength=10000.0;
    env2_signal_strength=10000.0;
    env1_t_signal_on=200.0;    
    env1_t_signal_off=0.0;     
    env2_t_signal_on=10.0;
    env2_t_signal_off=130.0;
    env1_initial_effect_of_effector='b';
    env2_initial_effect_of_effector='d';
    env1_fixed_effector_effect=0;    
    env2_fixed_effector_effect=1;                 // global variable
    recalc_new_fitness=5;                   // global variable, make sure its value is smaller than MAX_RECALC_FITNESS        
    env1_occurence=0.33;                     // global variable
    env2_occurence=0.67;                     // global variable
    DUPLICATION=1.5e-7;                 
    SILENCING = 1.5e-7;
    N_EFFECTOR_GENES=EFFECTOR_GENES;
    N_TF_GENES=TFGENES;
    miu_PIC_disassembly=5.22;
    miu_mRNA_decay=-1.13;
    miu_protein_decay=-4.58;
    miu_translation_init=-1.36;
    
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
    fprintf(fp,"env1: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env1_t_signal_on, env1_t_signal_off, env1_initial_effect_of_effector, env1_fixed_effector_effect, env1_occurence);
    fprintf(fp,"env2: signal on duration=%f min, signal off duration=%f min, initial effector effect=%c, always_deleterious_effector:%d occurrence=%f\n",env2_t_signal_on, env2_t_signal_off, env2_initial_effect_of_effector, env2_fixed_effector_effect, env2_occurence);       
    fprintf(fp,"Background signal strength=%f\n",background_signal_strength);
    fprintf(fp,"Signal off strength=%f, env1 signal on strength=%f, env2 signal on strength=%f \n",signal_off_strength,env1_signal_strength,env2_signal_strength);
    fprintf(fp,"env1 init effecto effect %c, fixed? %d. env2 init effecto effect %c, fixed? %d.\n",env1_initial_effect_of_effector,env1_fixed_effector_effect,env2_initial_effect_of_effector,env2_fixed_effector_effect);
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

void calc_fitness_stats(Genotype *genotype,
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
#if ADJUST_FITNESS
    float diff;
    diff=fabs((float)(ADJUST_FITNESS_CONDITION-genotype->ntfgenes+N_SIGNAL_TF));
//    if(diff>0.0)    
        genotype->fitness=genotype->fitness*(1.0-ADJUSTMENT*diff);    
#endif
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
            #if !SET_BS_MANUALLY		
                calc_all_binding_sites(genotype_ori_copy);
            #endif
            MAXELEMENTS=genotype_ori_copy->N_allocated_elements;
            /*calculate the fitness of the mutant at low resolution*/
            calc_avg_growth_rate(   genotype_ori_copy,
                                    init_mRNA,
                                    init_protein_number,
                                    RS_parallel,                                   
                                    i,
                                    GR1[0],
                                    GR2[0],
                                    mut_record);
            calc_fitness_stats( genotype_ori_copy,
                                &(GR1[0]),
                                &(GR2[0]),
                                1); // N_recalc_fitness=1 implicitly when calculating the fitness of a mutant.
            /*record low-resolution fitness*/
            fp=fopen("Mut_detail_fitness.txt","a+");
            fprintf(fp,"%.10f %.10f\n",
                    genotype_ori_copy->fitness,
                    genotype_ori_copy->sq_SE_fitness);
            fclose(fp);        
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
        for(j=1;j<recalc_new_fitness;j++)                    
            calc_avg_growth_rate(   genotype_ori, 
                                    init_mRNA,
                                    init_protein_number,
                                    RS_parallel,                                    
                                    0,
                                    GR1[j],
                                    GR2[j],
                                    mut_record);  
        calc_fitness_stats( genotype_ori,
                            &(GR1[0]),
                            &(GR2[0]),
                            recalc_new_fitness); 
        #if !SET_BS_MANUALLY
            calc_all_binding_sites(genotype_ori); 
        #endif            
        /*calculate the number of c1-ffls every step*/
        find_ffl(genotype_ori); 
        print_core_c1ffls(genotype_ori); 
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
//    /* initialize min_act_to_transcribe*/
//    for(i=0;i<N_THREADS;i++)
//    {
//        for(j=0;j<MAX_BINDING;j++)
//        {
//            min_act_to_transcribe[i][j]=(int)ceil(fabs((j-0.31)/0.33)); // determine the number of activators needs for transcription 
//                                                                    // based the number of repressors, j
//        }
//    }
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
    genotype_ori_copy.ntfgenes=genotype_ori.ntfgenes;
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
//                                 RS_parallel[2*N_THREADS],
//                                 &(RS_parallel[2*N_THREADS+1])
                                RS_main,
                                RS_parallel
                                );                                                                    
            
        }  
    }
    else /* otherwise the simulation starts over from beginning*/
    {   
        /* record the initial network topology*/
        summarize_binding_sites(&genotype_ori,init_step); /*snapshot of the initial (0) distribution binding sites */   
        find_ffl(&genotype_ori); 
        print_core_c1ffls(&genotype_ori);

    #if JUST_PLOTTING 
        fp=fopen("MUT_8.txt","r");    
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
                            60000);
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
            env1_t_development=89.9;
            env2_t_development=89.9;
            duration_of_burn_in_growth_rate=0.0;
            env1_signal_strength=10000.0;
            env2_signal_strength=10000.0;
            env1_t_signal_on=200.0;    
            env1_t_signal_off=0.0;     
            env2_t_signal_on=10.0;
            env2_t_signal_off=130.0;
            env1_initial_effect_of_effector='b';
            env2_initial_effect_of_effector='d';
            env1_fixed_effector_effect=0;    
            env2_fixed_effector_effect=1; 
            recalc_new_fitness=5; // make sure its value is smaller than MAX_RECALC_FITNESS
            env1_occurence=0.33;
            env2_occurence=0.67;    
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
                calc_avg_growth_rate(   &genotype_ori, 
                                        init_mRNA,
                                        init_protein_number,
                                        RS_parallel,     
//                                         &(RS_parallel[2*N_THREADS+1]),
                                        0,
                                        GR1[i],
                                        GR2[i],
                                        &mut_record);                
            }
                                        
            calc_fitness_stats(&genotype_ori,&(GR1[0]),&(GR2[0]),recalc_new_fitness);
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
//                         RS_parallel[2*N_THREADS],
//                         &(RS_parallel[2*N_THREADS+1])
                        RS_main,
                        RS_parallel
                                );
    
        #else // no selection at all, just randomly shuffle network topology and kinetic parameters            
        int max_mut_steps=MAX_MUT_STEP;
        evolve_neutrally(   &genotype_ori,
                            &genotype_ori_copy,                          
                            &mut_record,
                            max_mut_steps,
                            kdis,
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
//                    fprintf(OUTPUT,"%c",genotype->tf_seq[k][j]);
//            fprintf(OUTPUT,"\n");
//            fclose(OUTPUT);
//            snprintf(filename,sizeof(char)*32,"TF_r_%i.txt",k);
//            OUTPUT=fopen(filename,"a+");
//            fprintf(OUTPUT,"%d",step_i);
//            for(j=0;j<TF_ELEMENT_LEN;j++)
//                    fprintf(OUTPUT,"%c",genotype->tf_seq_rc[k][j]);
//            fprintf(OUTPUT,"\n");
//            fclose(OUTPUT);
//    }
}

void print_core_c1ffls(Genotype *genotype)
{
    FILE *fp; 
    fp=fopen("proportion_c1ffl.txt","a+");
    fprintf(fp,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                                                                                                genotype->N_motifs[0],
                                                                                                genotype->N_motifs[1],
                                                                                                genotype->N_motifs[2],
                                                                                                genotype->N_motifs[3],
                                                                                                genotype->N_motifs[4],
                                                                                                genotype->N_motifs[5],
                                                                                                genotype->N_motifs[6],
                                                                                                genotype->N_motifs[7],
                                                                                                genotype->N_motifs[8],
                                                                                                genotype->N_motifs[9],
                                                                                                genotype->N_motifs[10],
                                                                                                genotype->N_motifs[11],
                                                                                                genotype->N_motifs[12],
                                                                                                genotype->N_motifs[13],
                                                                                                genotype->N_motifs[14],
                                                                                                genotype->N_motifs[15],
                                                                                                genotype->N_motifs[16],
                                                                                                genotype->N_motifs[17],
                                                                                                genotype->N_motifs[18],
                                                                                                genotype->N_motifs[19],
                                                                                                genotype->N_motifs[20],
                                                                                                genotype->N_motifs[21],
                                                                                                genotype->N_motifs[22],
                                                                                                genotype->N_motifs[23],
                                                                                                genotype->N_motifs[24],
                                                                                                genotype->N_motifs[25],
                                                                                                genotype->N_motifs[26],
                                                                                                genotype->N_act_genes,
                                                                                                genotype->N_act_genes_reg_by_env,
                                                                                                genotype->N_act_genes_not_reg_by_env,
                                                                                                genotype->protein_pool[genotype->nproteins-1][0][0]);
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
    int i,j;
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
    for(i=0;i<genotype->nproteins-1;i++)
    {
        if(genotype->activating[i]==1)
        {
            fprintf(OUTPUT1," A%d ",i);
//            N_act++;
        }
        if(genotype->activating[i]==0)
        {
            fprintf(OUTPUT1," R%d ",i);
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
        
        for(j=0;j<genotype->nproteins-1;j++)
        {
            if(table[i][j]<10)
                fprintf(OUTPUT1," %d  ",table[i][j]);
            else
                fprintf(OUTPUT1," %d ",table[i][j]);
        }
        if(genotype->which_protein[i]==genotype->nproteins-1)
            fprintf(OUTPUT1,"S");
        else
        {
            if(genotype->activating[genotype->which_protein[i]]==1)
                fprintf(OUTPUT1,"A%d",genotype->which_protein[i]); 
            if(genotype->activating[genotype->which_protein[i]]==0)
                fprintf(OUTPUT1,"R%d",genotype->which_protein[i]);
        }
        fprintf(OUTPUT1," %d \n",genotype->min_act_to_transc[i]);
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
//        if(genotype->activating[i]==1)
//        {
//            fprintf(OUTPUT1," A%d ",i);
////            N_act++;
//        }
//        if(genotype->activating[i]==0)
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
//            if(genotype->activating[genotype->which_protein[i]]==1)
//                fprintf(OUTPUT1,"A%d",genotype->which_protein[i]); 
//            if(genotype->activating[genotype->which_protein[i]]==0)
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
        fprintf(fp,"0 %d %d %d\n",genotype->all_binding_sites[which_gene][i].tf_id,genotype->all_binding_sites[which_gene][i].BS_pos+8,genotype->activating[genotype->all_binding_sites[which_gene][i].tf_id]);
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
    
    for(i=0;i<genotype->ntfgenes;i++)
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

void find_ffl(Genotype *genotype)
{
    int i,j,k,l,cluster_size;
    int found_bs;
    int gene_id,gene_id_copy,site_id,protein_id,N_copies,N_activators;
    int master_TF,aux_TF;
//    int repressors[NPROTEINS];
    int activators[NPROTEINS];
    int copies_reg_by_env[NGENES],copies_not_reg_by_env[NGENES],N_copies_reg_by_env,N_copies_not_reg_by_env;
    int hindrance[NPROTEINS][NPROTEINS];
    int pos_binding_sites_of_j[MAXELEMENTS],pos_binding_sites_of_k[MAXELEMENTS],N_binding_sites_of_j,N_binding_sites_of_k;
    int N_all_motifs;
    int N_motifs;
    
#if FORCE_OR_GATE   
    for(i=0;i<NGENES;i++)
    {
        genotype->gene_in_core_C1ffl[i]=0;
        for(j=0;j<NPROTEINS;j++)
            genotype->TF_in_core_C1ffl[i][j]=0;
    }
#endif  
    
    if(genotype->N_act>N_SIGNAL_TF)
    {        
        for(i=0;i<27;i++)
            genotype->N_motifs[i]=0; 
        i=0;   
        while(genotype->cisreg_cluster[i][0]!=-1) 
        {    
            N_copies_reg_by_env=0;
            N_copies_not_reg_by_env=0;
            for(j=0;j<NGENES;j++)
            {
                copies_reg_by_env[j]=-1;
                copies_not_reg_by_env[j]=-1;
            }
           
            gene_id=genotype->cisreg_cluster[i][0];
            gene_id_copy=gene_id;
            if(genotype->which_protein[gene_id]==genotype->nproteins-1) // is a effector gene
            {
                for(j=0;j<NPROTEINS;j++)
                {
//                    repressors[j]=0;    
                    activators[j]=0;  
                }
                cluster_size=0;
                while(genotype->cisreg_cluster[i][cluster_size]!=-1)
                    cluster_size++;        

                /*scan binding sites for tfs that regulate gene_id*/
                for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                {
                    protein_id=genotype->all_binding_sites[gene_id][j].tf_id;
                    if(genotype->activating[protein_id]==1) // is a binding site of an activator 
                        activators[protein_id]=1;
//                    else
//                        repressors[protein_id]=1;
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
//                k=0;
//                j=0;
//                N_repressors=0;
//                while(j<genotype->nproteins)
//                {
//                    if(repressors[j]!=0)               
//                    {
//                        repressors[k]=j;
//                        if(k!=j)
//                            repressors[j]=0;
//                        k++; 
//                        N_repressors++;
//                    }    
//                    j++;
//                } 
                
                /*build a table to show hindrance between binding sites on effector gene*/
                /*if every binding site of j can hinder all the binding site of k, denote hindrance[j][k]=1. Otherwise 0*/
                /*For hindrance[j][j]=1, this means at most one binding site of j can be bound by TF j */
                /*A strict AND gate between j and K should have H[j][j]=H[k][k]=1, and H[j][k]=0*/
                for(j=0;j<NPROTEINS;j++)
                {
                    for(l=0;l<NPROTEINS;l++)
                        hindrance[j][l]=1;
                }
                for(j=0;j<N_activators;j++)
                {
                    /*list the positions of all the binding sites of j*/
                    for(l=0;l<MAXELEMENTS;l++)
                        pos_binding_sites_of_j[l]=-CISREG_LEN; 
                    N_binding_sites_of_j=0;
                    gene_id=genotype->cisreg_cluster[i][0];
                    for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                    {
                        if(genotype->all_binding_sites[gene_id][site_id].tf_id==activators[j])
                        {
                            pos_binding_sites_of_j[N_binding_sites_of_j]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                            N_binding_sites_of_j++;
                        }
                    } 
                    if(N_binding_sites_of_j==1 && genotype->min_act_to_transc[gene_id_copy]==2)                        
                    {
                        hindrance[activators[j]][activators[j]]=1;
                    }
                    else
                    {
                        /*if the first bs hinders the last bs, we can be sure at most one binding site of j can be bound by TF j */
                        if(pos_binding_sites_of_j[N_binding_sites_of_j-1]-pos_binding_sites_of_j[0]<TF_ELEMENT_LEN+2*HIND_LENGTH)                     
                            hindrance[activators[j]][activators[j]]=1;
                        else
                            hindrance[activators[j]][activators[j]]=0;
                    }
                    
                    /*list the positions of all the binding sites of k*/                  
                    N_binding_sites_of_k=0;
                    for(k=j+1;k<N_activators;k++)
                    {
                        for(l=0;l<MAXELEMENTS;l++)
                            pos_binding_sites_of_k[l]=-CISREG_LEN;
                        N_binding_sites_of_k=0;
                        gene_id=genotype->cisreg_cluster[i][0];                     
                        for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                        {
                            if(genotype->all_binding_sites[gene_id][site_id].tf_id==activators[k])
                            {
                                pos_binding_sites_of_k[N_binding_sites_of_k]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                N_binding_sites_of_k++;
                            } 
                        }                    
                        if(pos_binding_sites_of_j[N_binding_sites_of_j-1]-pos_binding_sites_of_k[0]>= TF_ELEMENT_LEN+2*HIND_LENGTH ||
                            pos_binding_sites_of_k[N_binding_sites_of_k-1]-pos_binding_sites_of_j[0]>= TF_ELEMENT_LEN+2*HIND_LENGTH    )
                        {
                            hindrance[activators[j]][activators[k]]=0; 
                            hindrance[activators[k]][activators[j]]=0;
                        }
                    }
                }
                
                /*make lists of gene copies that are regulated by environmental signal and of those that are not*/
                for(j=0;j<N_activators;j++)
                {
                    if(activators[j]>=N_SIGNAL_TF)
                    {  
                        N_copies=genotype->protein_pool[activators[j]][0][0]; 
                        for(k=0;k<N_copies;k++)
                        {
                            gene_id=genotype->protein_pool[activators[j]][1][k];
                            found_bs=0;
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1)
                                {
                                    found_bs=1;
                                    break;
                                }
                            } 
                            if(found_bs)
                            {
                                copies_reg_by_env[N_copies_reg_by_env]=gene_id;
                                N_copies_reg_by_env++;                                
                            }
                            else
                            {
                                copies_not_reg_by_env[N_copies_not_reg_by_env]=gene_id;
                                N_copies_not_reg_by_env++;                                 
                            }
                        }                        
                    }
                } 

                /*******************************count c1-ffls ************************/
                #if DIRECT_REG
                    /*count c1-ffls formed by environment tf and one env-regulated copy */
                    if(activators[0]==N_SIGNAL_TF-1 || activators[1]==N_SIGNAL_TF-1)
                    {
                        N_motifs=cluster_size*N_copies_reg_by_env;
                        genotype->N_motifs[0]+=N_motifs; 
                #if FORCE_OR_GATE
                        if(N_copies_reg_by_env!=0)
                            genotype->gene_in_core_C1ffl[gene_id_copy]=1;
                #endif
               
                        for(j=0;j<N_copies_reg_by_env;j++)
                        {
                            protein_id=genotype->which_protein[copies_reg_by_env[j]];
                            if(hindrance[N_SIGNAL_TF-1][protein_id])
                            {
                                if(hindrance[N_SIGNAL_TF-1][N_SIGNAL_TF-1])
                                {
                                    if(hindrance[protein_id][protein_id])
                                        genotype->N_motifs[1]+=cluster_size;
                                    else
                                        genotype->N_motifs[2]+=cluster_size;
                                }   
                                else
                                {
                                    if(hindrance[protein_id][protein_id])
                                        genotype->N_motifs[3]+=cluster_size;
                                    else
                                        genotype->N_motifs[4]+=cluster_size;
                                }
                            }
                            else
                            {
                                if(hindrance[N_SIGNAL_TF-1][N_SIGNAL_TF-1])
                                {
                                    if(hindrance[protein_id][protein_id])
                                        genotype->N_motifs[5]+=cluster_size;
                                    else
                                        genotype->N_motifs[6]+=cluster_size;
                                }   
                                else
                                {
                                    if(hindrance[protein_id][protein_id])
                                        genotype->N_motifs[7]+=cluster_size;
                                    else
                                        genotype->N_motifs[8]+=cluster_size;
                                }
                            }
                        }
                    } 
                #endif

                /*count c1-ffl formed by one env-regulated copy and one unregulated copy*/
                 for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=0;k<N_copies_not_reg_by_env;k++)
                    {
                        /*search bs of k on j*/
                        gene_id=copies_reg_by_env[j];
                        found_bs=0;
                        protein_id=genotype->which_protein[copies_not_reg_by_env[k]];
                        for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                        {
                            if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id)
                            {
                                found_bs=1;
                                break;
                            }
                        }
                        if(!found_bs) // j is not regulated by k  
                        {
                            /*search bs of j on k*/
                            gene_id=copies_not_reg_by_env[k];
                            found_bs=0;
                            protein_id=genotype->which_protein[copies_reg_by_env[j]];
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id)
                                {
                                    found_bs=1;
                                    break;
                                }
                            }
                            if(found_bs)                            
                            {    
                                genotype->N_motifs[9]+=cluster_size; 
                                master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                aux_TF=genotype->which_protein[copies_not_reg_by_env[k]];
                                if(hindrance[master_TF][aux_TF])
                                {
                                    if(hindrance[master_TF][master_TF])
                                    {
                                        if(hindrance[aux_TF][aux_TF])
                                            genotype->N_motifs[10]+=cluster_size;
                                        else  
                                            genotype->N_motifs[11]+=cluster_size;
                                    }
                                    else
                                    {
                                        if(hindrance[aux_TF][aux_TF])
                                            genotype->N_motifs[12]+=cluster_size;
                                        else  
                                            genotype->N_motifs[13]+=cluster_size;
                                    }    
                                }
                                else
                                {
                                    if(hindrance[master_TF][master_TF])
                                    {
                                        if(hindrance[aux_TF][aux_TF])
                                            genotype->N_motifs[14]+=cluster_size;
                                        else  
                                            genotype->N_motifs[15]+=cluster_size;
                                    }
                                    else
                                    {
                                        if(hindrance[aux_TF][aux_TF])
                                            genotype->N_motifs[16]+=cluster_size;
                                        else  
                                            genotype->N_motifs[17]+=cluster_size;
                                    }                                     
                                } 
                            }
                        }
                    }
                }
                /*count c1-ffls formed by two env-regulated copies*/
                /*also count parallel structures*/
                for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=j+1;k<N_copies_reg_by_env;k++)
                    {
                        if(k!=j)
                        {                    
                            /*search bs of k on j*/
                            gene_id=copies_reg_by_env[j];
                            found_bs=0;
                            protein_id=genotype->which_protein[copies_reg_by_env[k]];
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id)
                                {
                                    found_bs=1;
                                    break;
                                }
                            }
                            if(!found_bs) // j is not regulated by k
                            {
                                /*search bs of j on k*/
                                gene_id=copies_reg_by_env[k];
                                found_bs=0;
                                protein_id=genotype->which_protein[copies_reg_by_env[j]];
                                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(found_bs)  // k is regulated by j                          
                                {                                    
                                    genotype->N_motifs[18]+=cluster_size; 
                                    master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    aux_TF=genotype->which_protein[copies_reg_by_env[k]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[19]+=cluster_size;
                                            else  
                                                genotype->N_motifs[20]+=cluster_size;
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[21]+=cluster_size;
                                            else  
                                                genotype->N_motifs[22]+=cluster_size;
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[23]+=cluster_size;
                                            else  
                                                genotype->N_motifs[24]+=cluster_size;
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[25]+=cluster_size;
                                            else  
                                                genotype->N_motifs[26]+=cluster_size;
                                        }                                     
                                    }
                                }
                            }
                            else // j is regulated by k
                            {
                                /*search bs of j on k*/
                                gene_id=copies_reg_by_env[k];
                                found_bs=0;
                                protein_id=genotype->which_protein[copies_reg_by_env[j]];
                                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(!found_bs) // k is not regulated by j                           
                                {
                                    genotype->N_motifs[18]+=cluster_size; 
                                    master_TF=genotype->which_protein[copies_reg_by_env[k]];
                                    aux_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[19]+=cluster_size;
                                            else  
                                                genotype->N_motifs[20]+=cluster_size;
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[21]+=cluster_size;
                                            else  
                                                genotype->N_motifs[22]+=cluster_size;
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[23]+=cluster_size;
                                            else  
                                                genotype->N_motifs[24]+=cluster_size;
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[25]+=cluster_size;
                                            else  
                                                genotype->N_motifs[26]+=cluster_size;
                                        }                                     
                                    }                                   
                                }  
                            }
                        }
                    }
                }               
            }
            i++;
        }  
        genotype->N_act_genes=0;
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(genotype->activating[i]==1)
            {
                genotype->N_act_genes+=genotype->protein_pool[i][0][0];
            }
        }
//         #if DIRECT_REG    
//         for(i=0;i<9;i++)
//             genotype->proportion_motifs[i]=(float)genotype->N_motifs[i]/(genotype->N_act_genes*genotype->protein_pool[genotype->nproteins-1][0][0]);
//         #endif

        genotype->N_act_genes_reg_by_env=0;
        genotype->N_act_genes_not_reg_by_env=0;
        for(gene_id=N_SIGNAL_TF;gene_id<genotype->ngenes;gene_id++)
        {
            protein_id=genotype->which_protein[gene_id];
            if(genotype->activating[protein_id]==1)
            {
                found_bs=0;
                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                {
                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1)
                    {
                        found_bs=1;
                        break;
                    }
                }
                if(found_bs)
                    genotype->N_act_genes_reg_by_env++;
                else
                    genotype->N_act_genes_not_reg_by_env++;                            
            }
        }
//         N_all_motifs=genotype->protein_pool[genotype->nproteins-1][0][0]*(genotype->N_act_genes*(genotype->N_act_genes-1)/2);   
//         for(i=9;i<18;i++)
//             genotype->proportion_motifs[i]=(float)genotype->N_motifs[i]/N_all_motifs; 
    }    
    else
    {
        for(i=0;i<27;i++)
        {
//             genotype->proportion_motifs[i]=-1.0; 
            genotype->N_motifs[i]=0;            
        }
        genotype->N_act_genes=0; 
        genotype->N_act_genes_not_reg_by_env=0;
        genotype->N_act_genes_reg_by_env=0;
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
    
    fp1=fopen("proportion_c1ffl.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("proportion_c1ffl.txt");
    rename("temp","proportion_c1ffl.txt");
    
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

