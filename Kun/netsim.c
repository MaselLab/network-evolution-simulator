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
#include "random.h"
#include "lib.h"
#include "netsim.h"
#include "RngStream.h"

#define INITIALIZATION 2
#define UPDATE_PACT 1

enum {COPY_ALL=1,MUT_KCONST=2,MUT_CISSEQ=3,MUT_TFSEQ=4,MUT_N_GENE=5}; /*tags used in clone_cell*/

int MAXELEMENTS=100; 
/* start by allocating maxelements when initializing a genotype, double as needed, reduce at end */
const int MAXBOUND=500;
const int MAXALLOC=10;
const float KRNA=618.0;
const float TTRANSLATION=1.0; 
const float TTRANSCRIPTION=1.0;
const float PROB_ACTIVATING=0.62;
const float TRANSCRIPTINIT=5.0; /* replace betaon and betaoff */
const float DEACETYLATE=0.667;
const float ACETYLATE=0.167;
const float PICASSEMBLY=0.04; //revised from 0.0277
const float STARTNUCLEUS=0.1;
const float MIN_PICDISASSEMBLE = 0.1;
const float MAX_PICDISASSEMBLE = 1.0;
const float MAX_TRANSLATION = 600.0;   /* molecules per min. Based on 3 codons/s, 30 codons apart*/
const float MIN_TRANSLATION = 60.0 ;
const float MAX_mRNA_decay = 0.69;          /* 1 min half-life */
const float MIN_mRNA_decay = 5.8e-3;        /* 120 min */
const float MAX_protein_decay = 0.069;  /* Wang, 2002 PNAS 99:5860 */
const float MIN_protein_decay = 5.8e-3; 
const float KR=10.0;        /*this is the ratio of the affinity of a perfect binding site to that of a non-specific site*/
                            /*It controls the reduction in affinity per mismatch*/
const float lumped_kon=1.0e-2*0.9/6.02e23*1.0e6/3.0e-15;   /* use molecules of a tf to calculate on rate*/ 
                                                            /* kon=1.0e-2 uM^-1*s^-1 */
                                                            /* vol_nuc = 3.0e-15 l^3 */
                                                            /* NA = 6.02e23/1.0e6 convert number of molecules to umol*/
                                                            /* assuming 90% of tf molecules are translocated to nuclear */ 
const float koff=1.0e-2; /*s^-1, this is the initial value of koff*/
const float MIN_koff= 1.0e-5; /*s^-1, Fordecy 2010 Nat Biotech 28:970 reported Kd of 10 nM. N
                               * alefski 2006 Biochem 45:13794 reported 5.4nM*/
const float MAX_koff= 1.0; /*s^-1*/
float Koff[TF_ELEMENT_LEN-NMIN+1];
float protein_aging = 1e-4;/* protein aging term: used when c=c'+g=0, set to 1e-4 < mean-3*sd of
   Belle et al. (2006) and small with respect to most growth rates  */
const float MAX_RATE_TO_UPDATE_PACT=2.0;
const float MAX_CHANGE_IN_PACT=0.1;
/* end default options */ 

/*Mutations*/
float SUBSTITUTION = 0.33e-9; /* susbstitution rate per site per cell division. Lynch 2008*/
float INDEL = 0.02e-9;      /* indel rate per site per cell division Lynch 2008*/
/*the following mutations are enabled after burn-in*/
float DUPLICATION = 0.0;   
float SILENCING = 0.0;      
float MUTKINETIC = 2.6e-7*0.33; /* 50% subs and indel in a gene (~1.5kb, including introns) will change kinetic rates and binding seq */        
float proportion_mut_binding_seq = 0.1;
float proportion_mut_identity = 0.05;
float proportion_mut_koff=0.2;
float proportion_mut_kdis=0.05;
float proportion_mut_mRNA_decay=0.2;
float proportion_mut_protein_decay=0.2;
float proportion_mut_translation_rate=0.2;
float max_inset=3.0;
float max_delet=3.0;
float sd_mut_effect=0.1;
float bias_in_mut=0.0; /*mutations to kinetic constants are x*e^N(bias_in_mut,sd_mut_effect). */


/*Cell growth and development*/
int avg_protein_number = 12000;
float penalty = 1.0e-5;
float Pact_scaling = 1.0;
float tdevelopment;/* default  development time: can be changed at runtime */
float duration_of_burn_in_growth_rate; /* allow cells to reach (possiblly) steady growth*/
float growth_rate_scaling = 1.0; /* set default growth rate scaling factor */
int N_replicates;
int recalc_new_fitness; /*calculate the growth rate of the new genotype four more times to increase accuracy*/
float cost_term=1.0e-4;      /* this determined the cost of translation */                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
float penalty_of_extra_copies=1.0e-4 ;       /* this is the penalty for having more copies than the MAX_COPIES.
                                            * extra copies reduce growth rate. The amount of reducation per extra
                                            * copy is this number times gpeak (defined in initialized_growth_rate_parameters).
                                            * we use this penalty as a soft restriction to copy numbers. Note that if  
					    * the value of this term is similar to that of cost_term, then the reduction
                                            * in growth rate by an extra copy is numerically the same as the cost of expressing
                                            * the extra copy with a expression level of avg_protein_number*/
float reduction_in_P_dup=100.0;     /* this is another way to restrict copy numbers, by reducing the probabity
                                             * of duplication of a gene x folds, once its copy number exceed the upperlim.
                                             * reducing the probabity prevents even calculating the growth rate of the
                                             * new genotype*/

/*initial conditions*/
int init_TF_genes=6;
int init_N_act=3;
int init_N_rep=3;
int min_act_to_transcr_selection_protein=2;

/*Environments*/
float background_signal_strength = 10.0; /* number of TF0 protein */
float signal_strength=12000.0; /* number of each signal protein */
float basal_signal_strength=10.0; /* number of signal protein at basal level */
float env1_t_signalA;           /* duration of signal A in environment 1 */
float env1_t_signalB;           /* duration of signal B in environment 1 */
float env2_t_signalA;
float env2_t_signalB;
int env1_signalA_as_noise;      
int env2_signalA_as_noise;
float env1_occurence;             /* one environment can occur more frequently than*/                                
float env2_occurence;             /* the other*/ 

/* file output parameters */
char *output_directory = "output";   /* default output directory */
int verbose = 0;                     /* don't log verbosely by default */ 

/* initialize the growth rate parameters: 
 * do computations here so that we can easily change the scaling factor and Pp */
void initialize_growth_rate_parameters() {
  float hc, Ltf;
  gpeak = 0.005776*growth_rate_scaling;  /* in min^-1 based on doubling time of 120 min: ln(2)/(120 min)=0.005776 */
//  gpeak_b = 0.005776*growth_rate_scaling;
  Pp_s = 20000;
//  Pp_b = 20000;              /* mean gene expression of all proteins is 12064.28 */
  Ltf= 1418;               /* mean gene expression of only TFs is 1418 */
  hc = (gpeak/avg_protein_number)*(1-(log(2-2*cost_term)/log(2)));      /* in min^-1 cost of doubling gene expression, based on Wagner (2005) 
  * using {s=0.2, N=500} matches {s=10^-5, N=10^7} combination (both Ns=100). Geiler-Samerotte (2011) reported that 47,000 YFP 
  * molecules at equilibrium reduces yeast growth by 1.4%. Assuming the half-life of YFP equals to the doubling time of yeast (120 min), 
  * then the expression rate of YFP is about 268 molecules/min. This means when translation increases by 1 molecule/min, growth rate is 
  * reduced by 3e-7, i.e. h should equal to 3e-7, which is equivalent to setting the cost_term equals to 0.01.*/ 
  h = hc/0.023;            /* using c=0.023/min from mean of distribution from Belle et al (2006)*/
  gmax_a = gpeak + hc*(Pp_s+(10*Ltf));    /* compute the gmax coefficient based on gpeak and other parameters */
                                            /* Here we assume there will be 10 TFs protein in the end of the simulation. This can 
                                            * generate an interesting question: what is the ratio of protein content (or investment) 
                                            * between regulatory proteins and proteins that directly promote growth.*/
  gmax_b = gpeak+hc*10*Ltf;  
  h_extra_copy = gmax_a*penalty_of_extra_copies; /* this is the reduction in GR per extra copy of gene*/
  
 }

char set_base_pair(float x) {
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

void initialize_genotype_fixed( Genotype *genotype, 
                                float kdis[],
                                RngStream RS)
{
    int i, j;
    /* the first two genes are "signal" that take the form of tfs, i.e. having a binding seq. 
     * They have no promoter, no mRNA, no cost of expression, but constant protein conc*/
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {
        genotype->min_act_to_transc[i]=1;
    #if RANDOM_INIT_KINETIC_CONST
        /* tf affinity */
        genotype->koff[i]=koff*exp(gasdev(RS)*sd_mut_effect);
        if(genotype->koff[i]<MIN_koff || genotype->koff[i]>MAX_koff)
            genotype->koff[i]=koff*exp(gasdev(RS)*sd_mut_effect);
        /* mRNA decay */
        genotype->mRNAdecay[i] = exp(0.4909*gasdev(RS)-3.20304);
        while (genotype->mRNAdecay[i]<0.0)
            genotype->mRNAdecay[i] = exp(0.4909*gasdev(RS)-3.20304);
        /* protein decay */
        genotype->proteindecay[i]=-1.0;
        while (genotype->proteindecay[i] < 0.0) 
        {
            if (RngStream_RandU01(RS) < 0.08421)
                genotype->proteindecay[i] = EPSILON; /* was 0.0. changed to non-zero for stability*/
            else 
                genotype->proteindecay[i] = exp(0.7874*gasdev(RS)-3.7665);
        }
        /* translation rate */
        genotype->translation[i] = exp(0.7406*gasdev(RS)+4.56);
        while (genotype->translation[i] < 0.0)
            genotype->translation[i] = exp(0.7406*gasdev(RS)+4.56);
    #else
        genotype->koff[i]=koff;        
        genotype->mRNAdecay[i]=0.1386; //half life is about 5min
        genotype->proteindecay[i]=0.1386; 
        genotype->translation[i]=360; // 1 protein per 6s 
    #endif 
        j = RngStream_RandInt(RS,0,NUM_K_DISASSEMBLY-1);
        genotype->pic_disassembly[i] = kdis[j];
    }    
    /* make the activations the same in each copy */
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
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        genotype->mRNAdecay[i]=0.0;
        genotype->proteindecay[i]=0.0;
        genotype->translation[i]=0.0;
        genotype->pic_disassembly[i]=0.0; 
        genotype->activating[i]=1; /*make this signal an activator*/
        genotype->N_act++;   
        genotype->koff[i]=koff;
    }
#if RANDOMIZE_SIGNAL2
#if N_SIGNAL_TF==2
    if(RngStream_RandU01(RS)<=0.5)
        genotype->activating[1]=1;
    else
    {
        genotype->activating[1]=0;
        genotype->N_act--;
        genotype->N_rep++;
    }
#endif
#endif
    genotype->min_act_to_transc[genotype->ngenes-1]=min_act_to_transcr_selection_protein; 
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, but the same random TF sequence,
 * hindrance positions, replication times, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype,                        
                         float kdis[],
                         RngStream RS)
{ 
    int i,k;

    genotype->ngenes=1+N_SIGNAL_TF+init_TF_genes; /*including the signal genes and selection genes*/
    genotype->ntfgenes=N_SIGNAL_TF+init_TF_genes; /*including the signal genes*/
    genotype->nproteins=genotype->ngenes;  /*including the signal genes and selection genes*/
    
    for(i=0;i<genotype->ngenes;i++)
    {    
        genotype->which_cluster[i]=i;
        genotype->cisreg_cluster[i][0]=i;
    }
     
    /* initially, each protein has only one copy of gene and mRNA*/    
    for(i=0;i<genotype->nproteins;i++)
    {
        genotype->protein_pool[i][0][0]=1;
        genotype->protein_pool[i][1][0]=i;
        genotype->which_protein[i]=i;
    }
    
    genotype->ncopies_under_penalty=0;

    /* the first two genes are "signals" that take the forms of tf. They have no promoter, no mRNA, but constant protein conc*/
    initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*NGENES, genotype->ngenes, RS);
    initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*TFGENES, genotype->ntfgenes, RS); 
  /* We now generate the complementary sequence of BS that are on the non-template strand.
   * The complementary sequence is used to search for BS that on the non-template strand.  
   * We also assume that all the TFs have strong orientation preference.*/  
    for(i=0;i< genotype->ntfgenes;i++)
    {        
        for(k=0;k<TF_ELEMENT_LEN;k++)
        {
            switch (genotype->tf_seq[i][k])
            {
                case 'a': genotype->tf_seq_rc[i][k]='t'; break;
                case 't': genotype->tf_seq_rc[i][k]='a'; break;
                case 'c': genotype->tf_seq_rc[i][k]='g'; break;
                case 'g': genotype->tf_seq_rc[i][k]='c'; break;
            }
        }        
    }       
    initialize_genotype_fixed(genotype, kdis, RS); 
    
#if !SET_BS_MANUALLY    
    calc_all_binding_sites(genotype);
#endif
}

/*
 * Manually set binding sites
 *
 */
void set_binding_sites(Genotype *genotype,int *N_BS_per_gene, int (*tf_id_per_site)[3], int (*hindrance_per_site)[3], float (*koff_per_site)[3])
{
    int i,j;
    
    for(i=2;i<genotype->ngenes;i++)
    {
        genotype->binding_sites_num[i]=N_BS_per_gene[i-2];
        genotype->max_hindered_sites[i]=0;
        for(j=0;j<N_BS_per_gene[i-2];j++)
        {
            genotype->all_binding_sites[i][j].tf_id=tf_id_per_site[i-2][j];
            genotype->all_binding_sites[i][j].N_hindered=hindrance_per_site[i-2][j];
            genotype->all_binding_sites[i][j].Koff=koff_per_site[i-2][j];
            genotype->max_hindered_sites[i]=(genotype->max_hindered_sites[i]>hindrance_per_site[i-2][j])?
                                            genotype->max_hindered_sites[i]:hindrance_per_site[i-2][j];
            genotype->all_binding_sites[i][j].mis_match=0;
        }       
    }
}

/*
 * compute the list binding sites for specified gene and gene copy
 */
/*considering computing all the possible configurations here*/
void calc_all_binding_sites_copy(Genotype *genotype, int gene_id)
{
    int i, j, k, match,match_rc;
    int N_hindered_BS=0;
    int N_hindered_BS_copy;
    int FOUND_AT_THE_SAME_POS;
    int N_binding_sites=0;
    int start_TF;
    FILE *fp;
    genotype->N_act_BS[gene_id]=0;
    genotype->N_rep_BS[gene_id]=0;
    genotype->max_hindered_sites[gene_id]=0;  
    //some helper pointer 
    char *tf_seq;
    char *cis_seq;
    char *tf_seq_rc; 
    cis_seq=&(genotype->cisreg_seq[gene_id][0]); 
  
    for(i=0; i < CISREG_LEN-TF_ELEMENT_LEN; i++) /* scan promoter */
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
        /* set flags and copy*/
        N_hindered_BS_copy=N_hindered_BS;
        FOUND_AT_THE_SAME_POS=0;
        /* loop through TFs */
        /* NOTE: because of duplication, ntfgenes can be bigger than the number of actual tfs*/
        if(genotype->which_protein[gene_id]==genotype->nproteins-1) // the environmental signals cannot directly regulate the selection gene
            start_TF=2;
        else
            start_TF=0;
        
        for (k=start_TF; k < genotype->nproteins-1; k++) 
        {    
            tf_seq=&(genotype->tf_seq[k][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[k][0]);            
            /*find BS on the template strand*/
            match=0;
            for (j=i; j < i+TF_ELEMENT_LEN; j++) 
                if (cis_seq[j] == tf_seq[j-i]) match++; 
            if (match >= NMIN)
            {  
                if (N_binding_sites + 1 >= MAXELEMENTS) 
                {         
                    MAXELEMENTS+=MAXELEMENTS;
                    genotype->all_binding_sites[gene_id] = realloc(genotype->all_binding_sites[gene_id], MAXELEMENTS*sizeof(AllTFBindingSites));
                    if (!genotype->all_binding_sites[gene_id]) 
                    {
                        fp=fopen("output.txt","a+");
                        fprintf(fp,"error in calc_all_binding_sites_copy\n");
                        fclose(fp);
                        exit(1);
                    }
                    else 
                    {
                        fp=fopen("output.txt","a+");
                        fprintf(fp,"alloc more space for binding sites!\n");
                        fclose(fp);                        
                    }                        
                }                
                genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;
                genotype->all_binding_sites[gene_id][N_binding_sites].Koff = genotype->koff[k]*pow(KR,(float)((TF_ELEMENT_LEN-match)/(TF_ELEMENT_LEN-NMIN+1)));                
                genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ;  
//                genotype->all_binding_sites[gene_id][N_binding_sites].compressed = 0;
                
                if(!FOUND_AT_THE_SAME_POS) /* first BS found at this pos */                   
                {
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                    N_hindered_BS_copy++;
                    FOUND_AT_THE_SAME_POS=1;
                }
                else /* found BS of another TF at this pos again*/
                {
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS_copy;
                    N_hindered_BS_copy++;
                }
                N_binding_sites++;
                if(genotype->activating[k]==1) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand*/
            {
                match_rc=0;
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
                if (match_rc >= NMIN)
                {
                    /**********************************************************************/     
                    if (N_binding_sites + 1 >= MAXELEMENTS) 
                    {   
                        MAXELEMENTS+=MAXELEMENTS;
                        genotype->all_binding_sites[gene_id] = realloc(genotype->all_binding_sites[gene_id], MAXELEMENTS*sizeof(AllTFBindingSites));

                        if (!genotype->all_binding_sites[gene_id]) 
                        {
                            fp=fopen("output.txt","a+");
                            fprintf(fp,"error in calc_all_binding_sites_copy\n");
                            fclose(fp);                           
                            exit(1);
                        }
                        else
                        {
                            fp=fopen("output.txt","a+");
                            fprintf(fp,"alloc more space for binding sites!\n");
                            fclose(fp);                        
                        }      
                    }
                    /************************************************************************************************************/
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;
                    genotype->all_binding_sites[gene_id][N_binding_sites].Koff = genotype->koff[k]*pow(KR,(float)((TF_ELEMENT_LEN-match_rc)/(TF_ELEMENT_LEN-NMIN+1)));                   
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
//                    genotype->all_binding_sites[gene_id][N_binding_sites].compressed = 0;
                   
                    if(!FOUND_AT_THE_SAME_POS) /* first BS found at this pos */                   
                    {
                        genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                        N_hindered_BS_copy++;
                        FOUND_AT_THE_SAME_POS=1;
                    }
                    else /* found BS of another TF at this pos again*/
                    {
                        genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS_copy;
                        N_hindered_BS_copy++;
                    }                    
                    N_binding_sites++;                    
                    if(genotype->activating[k]==1) genotype->N_act_BS[gene_id]++;
                }
            } 
        }/* looping through TFs ends */
    }/*end of promoter scanning*/ 
    genotype->binding_sites_num[gene_id]=N_binding_sites;  
    genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);
    /* calc max_hindered_sites */
    for(i=0;i<N_binding_sites;i++)
    {
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id] > genotype->all_binding_sites[gene_id][i].N_hindered)?
                                      genotype->max_hindered_sites[gene_id] : genotype->all_binding_sites[gene_id][i].N_hindered;
    }
}

/*
 * compress a cluster of BS of the same TF into a single binding site  
 *
 */
void compress_binding_sites(Genotype *genotype, int gene_id)
{
//    int i,j,k;
//    int first_BS_in_cluster=0;
//    int last_BS_in_cluster=0;
//    int next_BS=1;
//    AllTFBindingSites temp[MAXELEMENTS];
//    int N_BS;
//    int current_tf;
//    int N_clusters=0;
//    int N_hindered=0;
//    int N_BS_evaluated=0;
//    int pos_left_hindrance,pos_right_hindrance;
//    int max_hindered_BS;
//    genotype->max_hindered_clusters[gene_id]=0;
//    FILE *fp;
//    char filename[32];
//    
//    for(j=0;j<NGENES;j++)
//    {
//        for(k=0;k<2;k++)
//            genotype->avg_N_BS_in_cluster[gene_id][j][k]=0;
//    }   
//    
//    while(N_BS_evaluated<genotype->binding_sites_num[gene_id]) /*loop through all the binding sites*/
//    {
//        N_BS=0;
//        next_BS=first_BS_in_cluster;
//        last_BS_in_cluster=first_BS_in_cluster;
//        current_tf=genotype->all_binding_sites[gene_id][first_BS_in_cluster].tf_id;
//        while(genotype->all_binding_sites[gene_id][next_BS].N_hindered>=(next_BS-last_BS_in_cluster) &&
//                next_BS<genotype->binding_sites_num[gene_id] && N_BS<MAX_BS_IN_CLUSTER)// considering compress at most 6 BS
//        {          
//            if(genotype->all_binding_sites[gene_id][next_BS].tf_id== current_tf) 
//            {
//                last_BS_in_cluster=next_BS;
//                temp[N_BS]=genotype->all_binding_sites[gene_id][next_BS];
//                genotype->all_binding_sites[gene_id][next_BS].compressed=1;
//                temp[N_BS].N_hindered=0;                
//                N_BS++; 
//                N_BS_evaluated++;
//            }
//            next_BS++;
//        }
//        genotype->avg_N_BS_in_cluster[gene_id][current_tf][0]++;
//        genotype->avg_N_BS_in_cluster[gene_id][current_tf][1]+=N_BS;
//        /* recalculate N_hindered for each BS in the cluster*/
//        max_hindered_BS=0;
//        for(i=0;i<N_BS;i++)
//        {
//            for(j=0;j<i;j++) 
//            {
//                if(temp[i].BS_pos-temp[j].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)           
//                    temp[i].N_hindered++;
//            }
//            max_hindered_BS=(max_hindered_BS>temp[i].N_hindered)?max_hindered_BS:temp[i].N_hindered;
//        }
//        /* update compressed_binding_sites*/
//        genotype->compressed_binding_sites[gene_id][N_clusters].tf_id=current_tf;
//        genotype->compressed_binding_sites[gene_id][N_clusters].off_flux=calc_flux( temp,
//                                                                                    N_BS, 
//                                                                                    max_hindered_BS,
//                                                                                    genotype->compressed_binding_sites[gene_id][N_clusters].coeff_on_flux);
//        genotype->compressed_binding_sites[gene_id][N_clusters].N_hindered=N_hindered;
//        genotype->compressed_binding_sites[gene_id][N_clusters].start_pos=genotype->all_binding_sites[gene_id][first_BS_in_cluster].BS_pos;  
//        genotype->compressed_binding_sites[gene_id][N_clusters].end_pos=genotype->all_binding_sites[gene_id][last_BS_in_cluster].BS_pos;
//        /*advance first_BS_in_cluster. */        
//        while(first_BS_in_cluster<genotype->binding_sites_num[gene_id]-1)
//        {
//            first_BS_in_cluster++;
//            if(genotype->all_binding_sites[gene_id][first_BS_in_cluster].tf_id!=current_tf || 
//                genotype->all_binding_sites[gene_id][first_BS_in_cluster].BS_pos>genotype->all_binding_sites[gene_id][last_BS_in_cluster].BS_pos) 
//            {
//                if(genotype->all_binding_sites[gene_id][first_BS_in_cluster].compressed==0)
//                {
//                    N_clusters++; /* this new cluster hinders the previous clusters*/
//                    break;
//                }
//            }  
//        }
//        /*calc clusters hindered by the new cluster*/
//        N_hindered=0;
//        for(i=0;i<N_clusters;i++)
//        {
//            pos_left_hindrance=genotype->compressed_binding_sites[gene_id][i].end_pos+TF_ELEMENT_LEN+2*HIND_LENGTH;
//            pos_right_hindrance=genotype->compressed_binding_sites[gene_id][i].start_pos;
//            if(genotype->all_binding_sites[gene_id][first_BS_in_cluster].BS_pos<=pos_left_hindrance &&
//                genotype->all_binding_sites[gene_id][first_BS_in_cluster].BS_pos>=pos_right_hindrance)
//                N_hindered++;
//        }
//        genotype->max_hindered_clusters[gene_id]=(genotype->max_hindered_clusters[gene_id]>N_hindered)?genotype->max_hindered_clusters[gene_id]:N_hindered;
//    }    
//    genotype->cluster_num[gene_id]=N_clusters+1;
//    for(i=0;i<genotype->binding_sites_num[gene_id];i++)
//        genotype->all_binding_sites[gene_id][i].compressed=0;
}

//double calc_flux(AllTFBindingSites *BS_info,
//                int N_BS,             
//                int max_hindered_BS,
//                double coeff_on_flux[MAX_MODE])
//{
//    double ratio_matrices[max_hindered_BS+1][MAX_MODE];   
//    double transition_matrix[MAX_MODE];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,m,n;
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_hindered_BS+1;i++)
//    {
//        for(j=0;j<MAX_MODE;j++)
//            ratio_matrices[i][j]=0.0;
//    }
//    for(j=0;j<MAX_MODE;j++)
//        transition_matrix[j]=0.0;   
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    ratio_matrices[pos_next_record][0]=BS_info[0].Koff; 
//    ratio_matrices[pos_next_record][1]=kon; 
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_hindered_BS+1);
//        product_of_freq = kon;
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//                product_of_freq*=BS_info[n].Koff;            
//        }
//       
//        if(m-BS_info[m].N_hindered!=0)
//        {
//          // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put at 0 and  
//          // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//            pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//            transition_matrix[0]=0.0;
//            for(j=1;j<MAX_MODE;j++)
//                transition_matrix[j]=ratio_matrices[pos_of_last_record][j-1]; 
//        }
//        else
//        {
//            for(i=0;i<MAX_MODE;i++)
//                transition_matrix[i]=0.0;            
//            transition_matrix[1]=1.0;
//        }
//
//        pos_of_last_record=mod(pos_next_record-1,max_hindered_BS+1);  //find last record              
//
//        for(i=0;i<MAX_MODE;i++)
//            ratio_matrices[pos_next_record][i]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i]+
//                                              product_of_freq*transition_matrix[i];
//    }
//    
//    for(i=1;i<MAX_MODE;i++)
//        coeff_on_flux[i]=ratio_matrices[pos_next_record][i];
//    coeff_on_flux[0]=0.0;
//    
//    return ratio_matrices[pos_next_record][0];   
//}

/*
 * compute the list of binding sites for the specified number of gene
 * copies
 */
void calc_all_binding_sites(Genotype *genotype)
{    
    int gene_id;
    for(gene_id=N_SIGNAL_TF;gene_id < genotype->ngenes;gene_id++)
    {        
        if(genotype->recalc_TFBS[gene_id]) /* do not calculate the binding sites if there's no mutation in the promoter or in TF binding seq*/
        {
            calc_all_binding_sites_copy(genotype,gene_id);
//            compress_binding_sites(genotype,gene_id);        
//            calc_configurations(genotype, gene_id);            
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

/* calc with approximation, but greatly enhance performance*/
/**further modified to handle collapsed BS clusters*/
//float calc_TF_dist_from_compressed_BS(CompressedBindingSites *BS_info,
//                                            int ntfgenes,
//                                            int max_N_hindered_BS,
//                                            int N_BS,                                           
//                                            int activating[NPROTEINS],
//                                            float protein_number[NGENES])
//{
//    double ratio_matrices[max_N_hindered_BS+1][MAX_MODE][MAX_MODE];   
//    double transition_matrix[MAX_MODE][MAX_MODE];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;    
//    float protein_number_power[ntfgenes][MAX_MODE];
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,k,m,n;
//    /*calc nth power of TF concentration*/
//    for(i=0;i<ntfgenes;i++)
//    {
//        protein_number_power[i][0]=1.0;
//        protein_number_power[i][1]=protein_number[i];
//        for(j=2;j<MAX_MODE;j++)
//            protein_number_power[i][j]=protein_number_power[i][j-1]*protein_number[i];
//    }
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_N_hindered_BS+1;i++)
//    {
//        for(j=0;j<MAX_MODE;j++)            
//        {
//            for(k=0;k<MAX_MODE;k++)
//                ratio_matrices[i][j][k]=0.0;
//        }
//    }
//    for(i=0;i<MAX_MODE;i++)
//    {
//        for(j=0;j<MAX_MODE;j++)
//                transition_matrix[i][j]=0.0;
//    }
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    ratio_matrices[pos_next_record][0][0]=BS_info[0].off_flux;   
//    
//    if(activating[BS_info[0].tf_id]==1) // if a activator binds to this BS
//    {
//        for(i=1;i<MAX_MODE;i++)
//            ratio_matrices[pos_next_record][0][i]=BS_info[0].coeff_on_flux[i]*protein_number_power[BS_info[0].tf_id][i];
//    }
//    else
//    {
//        for(i=1;i<MAX_MODE;i++)
//            ratio_matrices[pos_next_record][i][0]=BS_info[0].coeff_on_flux[i]*protein_number_power[BS_info[0].tf_id][i];
//    }    
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
//        product_of_freq = 1.0;
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//                product_of_freq*=BS_info[n].off_flux;            
//        }
//
//        switch(activating[BS_info[m].tf_id])
//        {
//            case 1: // a BS of activators
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
//                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//                    pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//                    // the next couple lines construct transition_matrix that has two properties:                    
//                    // 1. make elements of the transition_matrix, T_ij, from the elements of ratio_matrices[pos_of_last_record], R_ij
//                    // T_ij=Sigma_(k=1){R_i(j-k)*coeff_on_flux_of_mode_k*tf_concentration_raised_to_mode_k}
//                    // 2. moving all columns of transition_matrix 1 col to the right, and make an all-zero first column.
//                    for(i=0;i<MAX_MODE;i++) 
//                    {         
//                        transition_matrix[i][0]=0.0;
//                        for(j=1;j<MAX_MODE;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                            for(k=0;k<MAX_MODE;k++)
//                                transition_matrix[i][j]+=ratio_matrices[pos_of_last_record][i][j-k]*
//                                                         BS_info[m].coeff_on_flux[k]*
//                                                         protein_number_power[BS_info[m].tf_id][k];
//                        }
//                    }
//                }
//                else
//                {
//                    for(i=0;i<MAX_MODE;i++) 
//                    {  
//                        for(j=0;j<MAX_MODE;j++)                        
//                            transition_matrix[i][j]=0.0;                         
//                    }
//                    for(j=1;j<MAX_MODE;j++)
//                    {                       
//                        transition_matrix[0][j]=BS_info[m].coeff_on_flux[j]*protein_number_power[BS_info[m].tf_id][j];
//                    }
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              
//
//                for(i=0;i<MAX_MODE;i++)
//                {
//                    for(j=0;j<MAX_MODE;j++)                    
//                        ratio_matrices[pos_next_record][i][j]=  BS_info[m].off_flux*
//                                                                ratio_matrices[pos_of_last_record][i][j]+
//                                                                product_of_freq*transition_matrix[i][j];
//                }
//                break;
//
//            case 0: // a BS of repressors               
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                    pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); 
//                    for(j=0;j<MAX_MODE;j++)
//                        transition_matrix[0][j]=0.0;
//                    
//                    for(i=1;i<MAX_MODE;i++) 
//                    {   
//                        for(j=0;j<MAX_MODE;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                            for(k=0;k<MAX_MODE;k++)
//                                transition_matrix[i][j]+=ratio_matrices[pos_of_last_record][i-k][j]*
//                                                         BS_info[m].coeff_on_flux[k]*
//                                                         protein_number_power[BS_info[m].tf_id][k];
//                        }
//                    }
//                }
//                else
//                {
//                    for(i=0;i<MAX_MODE;i++) 
//                    {   
//                        for(j=0;j<MAX_MODE;j++)                        
//                            transition_matrix[i][j]=0.0;
//                    }
//                    for(i=1;i<MAX_MODE;i++)
//                    {                        
//                        transition_matrix[i][0]=BS_info[m].coeff_on_flux[i]*protein_number_power[BS_info[m].tf_id][i];
//                    }
//                }                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
//                
//                for(i=0;i<MAX_MODE;i++)
//                {                 
//                    for(j=0;j<MAX_MODE;j++)                    
//                        ratio_matrices[pos_next_record][i][j]=  BS_info[m].off_flux*
//                                                                ratio_matrices[pos_of_last_record][i][j]+
//                                                                product_of_freq*transition_matrix[i][j]; 
//                }
//                break;
//        }
//    }
//    sum=0.0;
//    for(i=0;i<MAX_MODE;i++)
//    {
//        for(j=0;j<MAX_MODE;j++)
//        {
//            sum+=ratio_matrices[pos_next_record][i][j];
//        }
//    }
//    for(i=0;i<MAX_MODE;i++)
//    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
//        
//        for(;j<MAX_MODE;j++)
//        {
//            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
//        }
//    } 
//    return (float)(prob_act_over_rep/sum);
//     // end of the forward algorithm   
//}
/**below is the original algorithm***/
/* calc with approximation, but greatly enhance performance*/

float calc_TF_dist_from_all_BS( AllTFBindingSites *BS_info,
                                int nproteins,
                                int max_N_hindered_BS,
                                int N_BS,
                                int N_act_BS,
                                int N_rep_BS, 
                                int activating[NPROTEINS],
                                int approximation,
                                float protein_number[NGENES],
                                int min_act_to_transcr)
{
    approximation=(approximation<(N_BS+1))?approximation:(N_BS+1);
    double ratio_matrices[max_N_hindered_BS+1][approximation][approximation];   
    double transition_matrix[approximation][approximation];
    double sum,prob_act_over_rep=0.0;    
    double product_of_freq;    
    float Kon[nproteins-1];
    int pos_of_last_record;    
    int pos_next_record;
    int i,j,k,m,n;
//    int n_col1=1;
//    int n_col2=1;
//    int n_row=1;
    
    /*calc Kon based on TF concentration*/
    for(i=0;i<nproteins-1;i++)    
        Kon[i]=lumped_kon*protein_number[i];
    
    /* initializing matrices to all zeros */
    for(i=0;i<max_N_hindered_BS+1;i++)
    {
        for(j=0;j<approximation;j++)            
        {
            for(k=0;k<approximation;k++)
                ratio_matrices[i][j][k]=0.0;
        }
    }
    for(j=0;j<approximation;j++)
    {
        for(k=0;k<approximation;k++)
            transition_matrix[j][k]=0.0;
    }
    
    /* body of the forward algorithm*/    
    pos_next_record=0; //where in the ratio_matrices to put the next record
    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
    
    if(activating[BS_info[0].tf_id]==1) // if a activator binds to this BS
    {
        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col2++;
    }
    else
    {
        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
    }    
    
    for(m=1;m<N_BS;m++)
    {
        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
        product_of_freq = Kon[BS_info[m].tf_id];

        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
        {
            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
                product_of_freq*=BS_info[n].Koff;            
        }

        switch(activating[BS_info[m].tf_id])
        {
            case 1: // a BS of activators
//                n_col2++;                
                if(m-BS_info[m].N_hindered!=0)
                {
                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
                    pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              

                    for(i=0;i<approximation;i++)
                    {
                        transition_matrix[i][0]=0.0;                      
//                      n_col1=(n_col2<N_act_bound[i])? n_col2:N_act_bound[i];
                        for(j=1;j<approximation;j++)
                            transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];                      
                    }
                }
                else
                {
                    for(i=0;i<approximation;i++)
                    {
//                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];                        
                        for(j=0;j<approximation;j++)
                            transition_matrix[i][j]=0.0;
                    }
                    transition_matrix[0][1]=1.0;
                }
                    
                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              

                for(i=0;i<approximation;i++)
                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
                    
                    for(j=0;j<approximation;j++)                    
                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
                                                      product_of_freq*transition_matrix[i][j];
                }
                break;

            case 0: // a BS of repressors
//                n_row++;//                
//                n_row=(n_row < max_N_rep_bound+1)?n_row:max_N_rep_bound+1;                
                if(m-BS_info[m].N_hindered!=0)
                {
                    pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); 
//                  n_col1=(n_col2<N_act_bound[0])?n_col2:N_act_bound[0];
                    for(j=0;j<approximation;j++)
                        transition_matrix[0][j]=0.0;
                  
                    for(i=1;i<approximation;i++)
                    {
//                      n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];                      
                        for(j=0;j<approximation;j++)
                            transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
                    }
                }
                else
                {
                    for(i=0;i<approximation;i++)
                    {
//                      n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
                        for(j=0;j<approximation;j++)                        
                            transition_matrix[i][j]=0.0;                       
                    }
                    transition_matrix[1][0]=1.0;
                }                    
                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
                
                for(i=0;i<approximation;i++)
                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];                     
                    for(j=0;j<approximation;j++)                    
                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
                                                      product_of_freq*transition_matrix[i][j]; 
                }
                break;
        }
    }

    sum=0.0;

    for(i=0;i<approximation;i++)
    {
        for(j=0;j<approximation;j++)
        {
            sum+=ratio_matrices[pos_next_record][i][j];
        }
    }   

    for(i=0;i<approximation;i++)
    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
        j=(min_act_to_transcr>round(fabs((i-0.31)/0.33)))?min_act_to_transcr:round(fabs((i-0.31)/0.33));
        for(;j<approximation;j++)
        {
            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
        }
    } 

    return (float)(prob_act_over_rep/sum);
     // end of the forward algorithm   
}

float calc_TF_dist_from_all_BS_simple(  AllTFBindingSites *BS_info,
                                        int nproteins,                                      
                                        int N_BS,
                                        int activating[NPROTEINS],                                        
                                        float protein_number[NGENES],
                                        int min_act_to_transcr)
{
    int i;
    float Kon[nproteins-1];
    float t_unbound_rep=1.0;
    float t_unbound_act=1.0;
    float t_unbound_special_tf=1.0;
    
    /*calc Kon based on TF concentration*/
    for(i=0;i<nproteins-1;i++)    
        Kon[i]=lumped_kon*protein_number[i];
    
    if(min_act_to_transcr==1 || activating[nproteins-1]==0)
    {
        for(i=0;i<N_BS;i++)
        {
            if(activating[BS_info[i].tf_id]==0) // a repressor        
                t_unbound_rep*=BS_info[i].Koff/(BS_info[i].Koff+Kon[BS_info[i].tf_id]);
            else
                t_unbound_act*=BS_info[i].Koff/(BS_info[i].Koff+Kon[BS_info[i].tf_id]);
        }
        return t_unbound_rep*(1.0-t_unbound_act);
    }
    else
    {
        for(i=0;i<N_BS;i++)
        {
            if(activating[BS_info[i].tf_id]==0) // a repressor        
                t_unbound_rep*=BS_info[i].Koff/(BS_info[i].Koff+Kon[BS_info[i].tf_id]);
            else
            {   
                if(BS_info[i].tf_id==nproteins-1) // this is the signal protein                
                    t_unbound_special_tf*=BS_info[i].Koff/(BS_info[i].Koff+Kon[BS_info[i].tf_id]);                
                else
                    t_unbound_act*=BS_info[i].Koff/(BS_info[i].Koff+Kon[BS_info[i].tf_id]);
            }
        }
        return t_unbound_rep*(1.0-t_unbound_special_tf)*(1.0-t_unbound_act);
    }
}

/**/
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
    newtime->gene_id = i;  
    newtime->time = t;
    pos = sls_store(newtime, start, last);
    return pos;
}

/**returns 0 if new fixed event won't happen concurrently with any exisiting event*/
int check_concurrence(  float t, 
                        FixedEvent *translation, 
                        FixedEvent *transcription, 
                        FixedEvent *envA, 
                        FixedEvent *envB,
                        FixedEvent *burn_in)
{   
    while(translation!=NULL)
    {
        if(t<=translation->time)            
            return 1;
        translation=translation->next;
    }
    while(transcription!=NULL)
    {
        if(t<=transcription->time)        
            return 1;
        transcription=transcription->next;
    }
    while(envA!=NULL)
    {
        if(t==(envA->time))        
            return 1; 
        envA=envA->next;
    }
    while(envB!=NULL)
    {
        if(t==envB->time)            
            return 1;
        envB=envB->next;
    }
    while(burn_in!=NULL)
    {
        if(t==burn_in->time)        
            return 1; 
        burn_in=burn_in->next;
    }
    return 0;
}

void add_time_point(float time,
                    float conc,
                    TimeCourse **start,
                    TimeCourse **last)
{
    TimeCourse *newtime;
    FILE *fp;

    newtime = malloc(sizeof*newtime);
    if (!newtime) 
    {
        fp=fopen("output.txt","a+");
        fprintf(fp,"error in add_time_point\n");
        fclose(fp);        
        exit(1);
    }
    newtime->time = time;
    newtime->concentration = conc;
    sls_store_end2(newtime, start, last);
}

//void add_fixed_event_end(int gene_id,                         
//                         float t,
//                         FixedEvent **start,
//                         FixedEvent **last)
//{
//  FixedEvent *newtime;
//  FILE *fp;
//  newtime = malloc(sizeof*newtime);
//  if (!newtime) {
//    fp=fopen("output.txt","a+");
//    fprintf(fp,"error in add_fixed_event_end\n");
//    fclose(fp);
//      
//    printf("error in add_fixed_event_end\n");
//    exit(1);
//  }
//  newtime->gene_id = gene_id; 
//  newtime->time = t;
////  LOG_VERBOSE_NOCELLID("adding event end at time=%f for gene=%d (copy=%d)\n", t, gene_id, gene_copy);
//  sls_store_end(newtime, start, last);
//}

/*delete linked table from anywhere*/
void delete_fixed_event(int gene_id,                      
                        int gene_specific_event_id,
                        FixedEvent **start,
                        FixedEvent **last)
{
    FixedEvent *info, *lastinfo=NULL;
    int j, done;
    FILE *fp;
  
    j = -1;
    done = 0;
    info = *start;
    while (info) 
    {
        if (info->gene_id==gene_id) 
        {
            j++;
            if (j == gene_specific_event_id) 
            {
                if (info == *start) 
                {
                    *start = info->next;
                    if (info == *last) 
                        *last = NULL;                       
                } 
                else 
                {
                    lastinfo->next = info->next;
                    if (info == *last) 
                        *last = lastinfo;
                }                
                done = 1;
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
    if (done == 0) 
    {
        fp=fopen("output.txt","a+");
        fprintf(fp,"error in delete_fixed_event\n");
        fclose(fp);        
        exit(1);
    }
    free(info);
}

void delete_fixed_event_start(FixedEvent **start,
                              FixedEvent **last)
{
    FixedEvent *info;

    info = *start;
    *start = info->next;
    if (*last == info) 
        *last = NULL;
    free(info);
}

/*
 * initialize the cell state with the specified initial protein
 * concentration, mean mRNA number and mRNA decay and whether to do
 * burn-in of high kon rate or not
 */
void initialize_cell(   CellState *state,
                        int ngenes,
                        int nproteins,
                        int protein_pool[NPROTEINS][2][NGENES], 
                        float mRNAdecay[NGENES],
                        int init_mRNA_number[NGENES],
                        float init_protein_number[NPROTEINS],                       
                        int plotting,
                        RngStream RS)
{
    int i, j;
    float t;
    int N_data_points; 
    
#if RANDOM_INIT_CELL
    int totalmRNA;
#endif

    /* start cell size at 1.0 */
    state->cell_size = 1.0;     
    state->cell_size_after_burn_in = 1.0;

    /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
    state->growth_rate = 0.0;    

    state->mRNA_transcr_time_end = NULL;
    state->mRNA_transcr_time_end_last = NULL;
    state->mRNA_transl_init_time_end = NULL;
    state->mRNA_transl_init_time_end_last = NULL;
    state->signalB_starts_end = NULL;
    state->signalB_starts_end_last = NULL;
    state->signalA_starts_end = NULL;
    state->signalA_starts_end_last = NULL;
    state->burn_in_growth_rate =NULL;
    state->burn_in_growth_rate_last=NULL;
    state->sampling_point_end=NULL;
    state->sampling_point_end_last=NULL;
    
    for (i=N_SIGNAL_TF; i < ngenes; i++) 
    {
        state->active[i]= NUC_NO_PIC;       
        state->mRNA_cyto_num[i]=init_mRNA_number[i];
        state->mRNA_transl_cyto_num[i]=0;
        state->mRNA_transcr_num[i]=0;
        state->last_Pact[i]=0.0;
    }
       
    /* initiate protein concentration*/
    for (i=N_SIGNAL_TF; i < ngenes; i++) 
        state->gene_specific_protein_number[i] = init_protein_number[i];
    
    for(i=N_SIGNAL_TF;i<nproteins;i++)
    {
        state->protein_number[i]=0.0;
        
        for(j=0;j<protein_pool[i][0][0];j++)
            state->protein_number[i]+=state->gene_specific_protein_number[protein_pool[i][1][j]];
//        state->protein_conc[i]=state->protein_number[i]/NA/vol_nuc;
    }
//    state->protein_conc[nproteins-2]=state->protein_number[nproteins-2]/NA/(vol_cell-vol_nuc);
//    state->protein_conc[nproteins-1]=state->protein_number[nproteins-1]/NA/(vol_cell-vol_nuc);
    
    /* deal with the signal--the fake tf*/
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        state->mRNA_cyto_num[i]=0;        
        state->mRNA_transcr_num[i]=0;
        state->mRNA_transl_cyto_num[i]=0;
        state->gene_specific_protein_number[i]=0.0;
    }    
    
    /*mark when to start calculating average growth rate*/
    add_fixed_event(-1,duration_of_burn_in_growth_rate,&(state->burn_in_growth_rate),&(state->burn_in_growth_rate_last)); 
    
    /*plot protein concentration and fitness vs time ?*/
    if(plotting==1)
    {        
        t=1.0;
        N_data_points=(int)tdevelopment;
        for(i=0;i<N_data_points;i++)
        {
            add_fixed_event(-1,t,&(state->sampling_point_end),&(state->sampling_point_end_last));
            t+=1.0;
        } 
    }    
}

void set_env(CellState *state, char env, float duration_env1, float duration_env2)
{
    float t=0.0;   

    if(env=='A')
        state->protein_number[1]=signal_strength;
    else            
        state->protein_number[1]=basal_signal_strength;
    state->protein_number[0]=background_signal_strength;
       
    while(t<tdevelopment)
    {
        if(env=='A')
        {
            add_fixed_event(-1,t+duration_env1,&(state->signalB_starts_end),&(state->signalB_starts_end_last));
            env='B';
            t=t+duration_env1;                    
        }    
        else
        {
            add_fixed_event(-1,t+duration_env2,&(state->signalA_starts_end),&(state->signalA_starts_end_last));
            env='A';
            t=t+duration_env2;
        }
    }    
}
/* 
 * change in the number of mRNAs in the cytoplasm affects the kon
 * rates, note we must have already updated protein_number first
 */
void change_mRNA_cytoplasm(int i,
                           Genotype *genotype,
                           CellState *state,
                           GillespieRates *rates)
{
  float salphc;  
  /* number of mRNAs in cytoplasm affects kon rates */
  salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / state->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  state->konvalues[i][KON_SALPHAC_INDEX] = salphc;
}

void calc_all_rates(Genotype *genotype,
                    CellState *state,
                    GillespieRates *rates,                     
                    int UPDATE_WHAT)
{
    int i,cluster_id,too_much_error;
    float salphc; 
    float protein_decay;
    
    if(UPDATE_WHAT==UPDATE_PACT || UPDATE_WHAT==INITIALIZATION) /*Since these values are also updated in update_protein_number_cell_size, this is only used in initilization*/
    {
        for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
        {
            /* if protein decay is otherwise going to fall below cut-off, use aging term */
            protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;
            salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);
            state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = protein_decay;
            state->konvalues[i][KON_SALPHAC_INDEX] = salphc;
        }
        for (i=0;i<N_SIGNAL_TF;i++)
        {
            state->konvalues[i][KON_PROTEIN_DECAY_INDEX]=0.0;
            state->konvalues[i][KON_SALPHAC_INDEX]=0.0;
        }
    }
    
    /* reset rates and operations */
//    rates->transport=0.0;
    rates->mRNAdecay=0.0;
    rates->pic_disassembly=0.0;
    rates->acetylation=0.0;
    rates->deacetylation=0.0;
    rates->transcript_init=0;
    rates->pic_assembly=0.0;
    rates->subtotal=0.0;
    
    for(i=0;i<genotype->ngenes;i++)
    {
        rates->acetylation_rate[i]=0.0;
        rates->deacetylation_rate[i]=0.0;
        rates->pic_assembly_rate[i]=0.0;
        rates->pic_disassembly_rate[i]=0.0;
        rates->mRNAdecay_rate[i]=0.0;
//        rates->transport_rate[i]=0.0; 
        rates->transcript_init_rate[i]=0;
        state->Pact[i]=0.0;
    } 

    /* update rates*/
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
         /* update RNA transport rate*/
//        rates->transport_rate[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
//        rates->transport += rates->transport_rate[i];
        rates->mRNAdecay_rate[i] = genotype->mRNAdecay[i] * ((float) state->mRNA_cyto_num[i] + (float) state->mRNA_transl_cyto_num[i]);
        rates->mRNAdecay += rates->mRNAdecay_rate[i];        
    }
    
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {    
        cluster_id=genotype->which_cluster[i];        
        if(genotype->cisreg_cluster[cluster_id][0]!=i)  /*if this gene is not the first of its cluster*/ 
            state->Pact[i]=state->Pact[genotype->cisreg_cluster[cluster_id][0]]; /* copy TF distribution from the first*/

        else /* otherwise, we need to calc the ratio*/
        {
            if(genotype->binding_sites_num[i]!=0)
#if !IGNORE_BS_OVERLAPPING
                state->Pact[i]=calc_TF_dist_from_all_BS(    genotype->all_binding_sites[i],
                                                            genotype->nproteins,
                                                            genotype->max_hindered_sites[i],
                                                            genotype->binding_sites_num[i],
                                                            genotype->N_act_BS[i],
                                                            genotype->N_rep_BS[i],
                                                            genotype->activating, 
                                                            MAX_MODE,
                                                            state->protein_number,
                                                            genotype->min_act_to_transc[i]);        
#else
                state->Pact[i]=calc_TF_dist_from_all_BS_simple( genotype->all_binding_sites[i],
                                                                genotype->nproteins,                                                                
                                                                genotype->binding_sites_num[i],                                                                
                                                                genotype->activating,                                                                 
                                                                state->protein_number,
                                                                genotype->min_act_to_transc[i]);
#endif
            else
                state->Pact[i]=0.0;
        }
        
        /* calc other rates*/
        switch (state->active[i])
        {
            case NUC_NO_PIC:
                rates->acetylation_rate[i]=state->Pact[i]*ACETYLATE*Pact_scaling;             
                rates->acetylation+=rates->acetylation_rate[i];
                rates->deacetylation_rate[i]=0.0;
                rates->pic_assembly_rate[i]=0.0;
                rates->pic_disassembly_rate[i]=0.0;
                break;
                
            case NO_NUC_NO_PIC:
                rates->deacetylation_rate[i]=DEACETYLATE*(1.0-state->Pact[i]);
                rates->deacetylation+=rates->deacetylation_rate[i];
                rates->pic_assembly_rate[i]=PICASSEMBLY*state->Pact[i];
//                rates->pic_assembly_rate[i]=PICASSEMBLY;
                rates->pic_assembly+=rates->pic_assembly_rate[i]; 
                rates->pic_disassembly_rate[i]=0.0;
                rates->acetylation_rate[i]=0.0;
                break;
                
            case PIC_NO_NUC: 
                rates->pic_disassembly_rate[i]=genotype->pic_disassembly[i];
                rates->pic_disassembly+=rates->pic_disassembly_rate[i]; 
                rates->transcript_init_rate[i]= 1;
                rates->transcript_init+=1;
                rates->deacetylation_rate[i]=0.0; 
                rates->acetylation_rate[i]=0.0;
                rates->pic_assembly_rate[i]=0.0;
                break;
        }
    }
    rates->subtotal+=rates->deacetylation;
    rates->subtotal+=rates->pic_assembly;
    rates->subtotal+=rates->acetylation;
//    rates->subtotal+=rates->transport;
    rates->subtotal+=rates->mRNAdecay;
    rates->subtotal+=rates->pic_disassembly;
    rates->subtotal+=(float)rates->transcript_init*TRANSCRIPTINIT;
    
    too_much_error=0;
    if(UPDATE_WHAT==UPDATE_PACT || UPDATE_WHAT==INITIALIZATION)
        rates->rate_update_Pact=MAX_RATE_TO_UPDATE_PACT;
    else
    {
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
        {
            if(state->last_Pact[i]!=0.0 && state->Pact[i]!=0.0)
            {
                if(fabs(state->Pact[i]-state->last_Pact[i])/state->last_Pact[i]>MAX_CHANGE_IN_PACT)
                {
                    too_much_error=1;  
                    break;                        
                }
            }
            if(state->last_Pact[i]!=0.0 && state->Pact[i]==0.0)
            {
                too_much_error=1;
                break;
            }
            if(state->last_Pact[i]==0.0 && state->Pact[i]!=0.0)
            {
                too_much_error=1;
                break;
            }
        }
        if(!too_much_error)
            rates->rate_update_Pact/=2.0;  
    }
    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
        state->last_Pact[i]=state->Pact[i];  
    
    rates->subtotal+=rates->rate_update_Pact;    
}

/* 
 * check to see if a fixed event ends within the time t
 *
 * returns:
 *  0 if there is no fixed event occuring before time t
 *  1 if a transcription event happens before time t
 *  2 if a translation event happens before time t
 *  3 if a gene replication event happens before time t
 */
 
int does_fixed_event_end(FixedEvent *mRNA_transl_init_time_end,
                         FixedEvent *mRNA_transcr_time_end,
                         FixedEvent *signalB_starts_end,
                         FixedEvent *signalA_starts_end,	
			 FixedEvent *burn_in_growth_rate,
			 float t) 
{
    int retval=0;
    float t1;
    float t2;
    float t3;
    float t4;
    float t5;    

    if(mRNA_transcr_time_end == NULL && mRNA_transl_init_time_end==NULL && signalB_starts_end == NULL && signalA_starts_end == NULL && burn_in_growth_rate == NULL)
    {
        retval =0;
    }
    else
    {
        t1 = mRNA_transcr_time_end ? mRNA_transcr_time_end->time : TIME_INFINITY;
        t2 = mRNA_transl_init_time_end ? mRNA_transl_init_time_end->time : TIME_INFINITY;
        t3 = signalB_starts_end ? signalB_starts_end->time : TIME_INFINITY;
        t4 = signalA_starts_end ? signalA_starts_end->time : TIME_INFINITY;
	t5 = burn_in_growth_rate ? burn_in_growth_rate->time : TIME_INFINITY;

        if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5))
	{
            retval = 1;	
        }
        else
        {
            if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5)) 
            {
                retval = 2;
            }
            else
            {
                if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5)) 
                {
                    retval = 3;
                }
                else
                {
                    if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5)) 
                    {
                        retval = 4;
                    }
                    else
                    {
			if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4))
			{
                            retval = 5;
			}
			else                     
                            retval = 0;
                    }
                }
            }
        }
    }
    return retval;
}


/*
 * end transcription: update the mRNAs in the nucleus, cytoplasm
 * etc. accordingly and delete the event from the queue
 */
void end_transcription(float *dt,
                        float t,
                        CellState *state,
                        GillespieRates *rates,
                        Genotype *genotype,
                        int *end_state,
                        char *error,
                        int mut_step,
                        Mutation *mut_record,
                        char *env)
{
    int i;
    float concurrent;
    float endtime;
    /* recompute the delta-t based on difference between now and the
       time of transcription end */
    *dt = state->mRNA_transcr_time_end->time - t;    
    update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, error, mut_step, mut_record);
    /* get the gene which is ending transcription */
    i = state->mRNA_transcr_time_end->gene_id;    
    /* increase number of mRNAs in nucleus */
//    (state->mRNA_nuclear_num[i])++;
    /* increase number of mRNAs that are initializing translation*/
    (state->mRNA_transl_cyto_num[i])++;
    /* decrease the number of mRNAs undergoing transcription */
    (state->mRNA_transcr_num[i])--;
    /* delete the fixed even which has just occurred */
    delete_fixed_event_start(&(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));
    /* add rate KRNA to transport and Gillespie rates */ 
    
    /*add transcription initialization event*/  
    endtime=t+*dt+TTRANSLATION;
    concurrent=check_concurrence(   endtime,
                                    state->mRNA_transl_init_time_end,
                                    state->mRNA_transcr_time_end,
                                    state->signalB_starts_end,
                                    state->signalA_starts_end,
                                    state->burn_in_growth_rate);
    while(concurrent)
    {
        endtime+=0.0001;
        concurrent=check_concurrence(   endtime,
                                        state->mRNA_transl_init_time_end,
                                        state->mRNA_transcr_time_end,
                                        state->signalB_starts_end,
                                        state->signalA_starts_end,
                                        state->burn_in_growth_rate);        
    }  
    
    add_fixed_event(i, endtime, &(state->mRNA_transl_init_time_end), &(state->mRNA_transl_init_time_end_last));
}

void end_translation_init(  Genotype *genotype, 
                            CellState *state, 
                            GillespieRates *rates, 
                            float *dt,
                            float t,
                            int *end_state,
                            char *error,
                            int mut_step,
                            Mutation *mut_record,
                            char *env)
{
    int i;   
    float salphc;
    
    *dt = state->mRNA_transl_init_time_end->time - t;         /* make dt window smaller */
    update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, error, mut_step, mut_record);
    /* get identity of gene that has just finished translating */
    i=state->mRNA_transl_init_time_end->gene_id;   

    /* there is one less mRNA that is initializing translation */
    (state->mRNA_transl_cyto_num[i])--;  
    /* delete the event that just happened */
    delete_fixed_event_start(&(state->mRNA_transl_init_time_end), &(state->mRNA_transl_init_time_end_last));
    
    /* there is one more mRNA that is now in cytoplasm */
    (state->mRNA_cyto_num[i])++;
    salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / state->konvalues[i][KON_PROTEIN_DECAY_INDEX];
    state->konvalues[i][KON_SALPHAC_INDEX] = salphc;
//    change_mRNA_cytoplasm(state->mRNA_transl_init_time_end->gene_id, genotype, state, rates); 
}

/*
 * compute tprime factor used in the integration of growth rate
 */
float compute_tprime(Genotype *genotype, CellState* state, float *number_of_selection_protein, float dt) 
{
    int n_copies;
    int i;    
      
    n_copies=genotype->protein_pool[genotype->nproteins-1][0][0];
    float as[n_copies],c[n_copies];
    for(i=0;i<n_copies;i++)
    {
        c[i]=state->konvalues[genotype->protein_pool[genotype->nproteins-1][1][i]][KON_PROTEIN_DECAY_INDEX];
        as[i]=c[i]*state->konvalues[genotype->protein_pool[genotype->nproteins-1][1][i]][KON_SALPHAC_INDEX];
    }       
    return rtsafe(  &calc_fx_dfx,
                    n_copies,
                    number_of_selection_protein,
                    as,
                    c,
                    0.0,
                    dt,
                    1.0/60.0);
}

/*
 * get integral for growth rate
 */
float compute_integral(Genotype *genotype, CellState *state, float *protein, float dt, float Pp)
{
    int i,n_copies,ids[NPROTEINS];
    float integral=0.0,ect1;    
   
    n_copies=genotype->protein_pool[genotype->nproteins-1][0][0];
    for(i=0;i<n_copies;i++)
        ids[i]=genotype->protein_pool[genotype->nproteins-1][1][i];
        
    for(i=0;i<n_copies;i++)
    {
        ect1=exp(-state->konvalues[ids[i]][KON_PROTEIN_DECAY_INDEX]*dt)-1;
    
        integral+=1.0/Pp*(state->konvalues[ids[i]][KON_SALPHAC_INDEX]*ect1/state->konvalues[ids[i]][KON_PROTEIN_DECAY_INDEX]-
                protein[i]*ect1/state->konvalues[ids[i]][KON_PROTEIN_DECAY_INDEX]+ 
                state->konvalues[ids[i]][KON_SALPHAC_INDEX]*dt);
    }
    return integral;
}

/*
 * return the instantaneous growth rate given the current cell state,
 * also return the integrated growth rate as a pointer
 */
float compute_growth_rate_dimer(float *integrated_growth_rate,
                                Genotype *genotype,
                                CellState *state,
                                float* number_of_selection_protein,                               
                                float t, 
                                float dt,
				char env,
                                int *end_state,
                                char *error,
                                int mut_step,
                                Mutation *mut_record)
{
    int i;
    float instantaneous_growth_rate=0.0;  /* this is returned from the function */
    float total_alpha_s = 0.0;
    float dt_prime, dt_rest;
    FILE *fp;
    float P_next_s=state->protein_number[genotype->nproteins-1];   
    float P_s=0.0;  
  
    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++) 
       P_s+=number_of_selection_protein[i];
                      
    switch (env)
    {   
        case 'A': /* protein A is necessary!*/
            /* choose the appropriate piecewise linear integral */           
            if ( (P_s >= Pp_s) && (P_next_s >= Pp_s) )  /* P > Pp throughout */
            {   
                *integrated_growth_rate = gmax_a * dt;
            }
            else if ((P_s <= Pp_s) && (P_next_s <= Pp_s)) /* P < Pp throughout */
            {    
                *integrated_growth_rate = gmax_a* compute_integral(genotype, state, number_of_selection_protein, dt, Pp_s);
            }
            else if ((Pp_s >= P_s) && (P_next_s >= Pp_s)) /* P < Pp up until t' then P > Pp */
            {                    
                dt_prime = compute_tprime(genotype,state,number_of_selection_protein,dt);                
                dt_rest = dt - dt_prime;                
                *integrated_growth_rate = gmax_a * compute_integral(genotype, state, number_of_selection_protein, dt_prime, Pp_s);                
                *integrated_growth_rate += gmax_a * dt_rest;
            }
            else if ((P_s >= Pp_s) && (Pp_s >= P_next_s)) /* P > Pp up until t' then P < Pp */
            {   
                dt_prime = compute_tprime(genotype,state,number_of_selection_protein,dt);                 
                dt_rest = dt - dt_prime;
                *integrated_growth_rate = gmax_a * dt_prime;                
                *integrated_growth_rate += gmax_a* compute_integral(genotype, state, number_of_selection_protein, dt_rest, Pp_s);
            }
            else 
            {               
                fp=fopen(error,"a+");
                fprintf(fp,"Step %d, %c %d %d %s %d %f error in compute_growth_rate_dimer\n",
                        mut_step,
                        mut_record->mut_type,
                        mut_record->pos_g,
                        mut_record->pos_n,
                        mut_record->nuc_diff,
                        mut_record->kinetic_type,
                        mut_record->kinetic_diff);
//                fprintf(fp,"env=%c, P_a=%f, P_next_a=%f, Pp_a=%f, P_b=%f, P_next_b=%f, Pp_b=%f\n",
//                        env,
//                        P_a,
//                        P_next_a,
//                        Pp_a,
//                        P_b,
//                        P_next_b,
//                        Pp_b);
                fclose(fp);               
                *end_state=0;
                return 0;
            }
            /* compute instantaneous growth rate at t */
            if (P_next_s < Pp_s)
                instantaneous_growth_rate = gmax_a*P_next_s/Pp_s;
            else
                instantaneous_growth_rate = gmax_a;
            break;
    
        case 'B': /* protein b is necessary! */            
            *integrated_growth_rate = gmax_b*dt-penalty*compute_integral(genotype, state, number_of_selection_protein, dt, 1.0);	
            instantaneous_growth_rate = gmax_b - penalty*P_next_s;            
            break;
    }
    /* compute the total cost of translation across all genes  */
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++)        
#if NO_REGULATION_COST
    if(genotype->which_protein[i]<genotype->nproteins-2)
#endif
        total_alpha_s += genotype->translation[i] * state->mRNA_cyto_num[i];
    
    /* add constant term for integrated rate */
    *integrated_growth_rate += -h * dt * (total_alpha_s);
    /* and instantaneous integrated rate */
    instantaneous_growth_rate += -h * (total_alpha_s);
#if !RdcPdup
    *integrated_growth_rate -= h_extra_copy*genotype->ncopies_under_penalty*dt;
    instantaneous_growth_rate -= h_extra_copy*genotype->ncopies_under_penalty;
#endif
    /* make sure growth rates can't be negative */
    if (*integrated_growth_rate < 0.0)
        *integrated_growth_rate = 0.0;
    if (instantaneous_growth_rate < 0.0)
        instantaneous_growth_rate = 0.0;
    /* return the instantaneous growth rate */
    return (instantaneous_growth_rate);
}

void calc_fx_dfx(float x, int n_copies, float *p, float *as, float *c,float *fx, float *dfx)
{
    int i;
    
    *fx=0;
    *dfx=0;    
    for(i=0;i<n_copies;i++)
    {
        *fx+=(p[i]-as[i]/c[i])*exp(-c[i]*x)+as[i]/c[i];
        *dfx+=(as[i]-c[i]*p[i])*exp(-c[i]*x);
    }
    *fx-=Pp_s;    
}

/* 
 * update both the protein concentration and current cell size
 *
 * need some sort of control in case it declines to essentially zero.
 * Add in discrete, stochastic and/or zero values, but this may create
 * false attractor if time steps are small and rising tide keeps
 * getting rounded down
 */
void update_protein_number_cell_size( Genotype *genotype,
                                    CellState *state,                                   
                                    GillespieRates *rates,
                                    float dt,                                                                      
                                    float t,                                   
                                    char env,
                                    int *end_state,
                                    char *error,
                                    int mut_step,
                                    Mutation *mut_record)
{
    int i,j;
    float ct, ect, ect1;
    float number_of_selection_protein[genotype->protein_pool[genotype->nproteins-1][0][0]];
    float instantaneous_growth_rate = 0.0;
    float integrated_growth_rate = 0.0;
//    float adjusted_decay;
    FILE *fp;
    
    /* store the concentration of the selection genes before updating*/   
    for(i=0;i<genotype->protein_pool[genotype->nproteins-1][0][0];i++)
        number_of_selection_protein[i]=state->gene_specific_protein_number[genotype->protein_pool[genotype->nproteins-1][1][i]];
 
    for (i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {
//        /* update protein decay rates due to dilution caused by growth */
//        adjusted_decay = genotype->proteindecay[i] + state->growth_rate;
//        /* if this results in a very small or zero decay rate, use protein aging term */
//        /* NOTE: we need to update this for every gene that encodes this protein*/
//        if (adjusted_decay > protein_aging)
//            state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = adjusted_decay;
//        else 
//            state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = protein_aging;
//      
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
    /* convert protein number to concentration */
//    for(i=2;i<genotype->nproteins-2;i++)
//        state->protein_conc[i]=state->protein_number[i]*prop_nuc_tf/NA/vol_nuc;
//    
//    state->protein_conc[genotype->nproteins-2]=state->protein_conc[genotype->nproteins-2]/NA/(vol_cell-vol_nuc); /*ignore the changes in cell volumn*/
//    state->protein_conc[genotype->nproteins-1]=state->protein_conc[genotype->nproteins-1]/NA/(vol_cell-vol_nuc);
//    
    /* now find out the protein concentration at end of dt interval and compute growth rate */   
    instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genotype, 
                                                            state, 
                                                            number_of_selection_protein,                                                             
                                                            t, 
                                                            dt,                                                                                                                     
                                                            env,
                                                            end_state,
                                                            error,
                                                            mut_step,
                                                            mut_record);  
	if(*end_state==0)
	{
		fp=fopen(error,"a+");
                fprintf(fp,"step %d error in compute_growth_rate_dimer",mut_step);
//		fprintf(fp,"mRNA_a=%d tranl_rate_a=%f mRNA_decay_a=%f protein_decay_a=%f salphc_a=%f c_a=%f ect_a=%f dt=%f\n",
//			state->mRNA_cyto_num[8],
//			genotype->translation[8],
//			genotype->mRNAdecay[8],
//			genotype->proteindecay[8],
//			state->konvalues[8][KON_SALPHC_INDEX],
//			state->konvalues[8][KON_PROTEIN_DECAY_INDEX],
//			ect1_s,
//			dt);
		fclose(fp);
		return;
	}  

    /* use the integrated growth rate to compute the cell size in the next timestep */
    state->cell_size = (state->cell_size)*exp(integrated_growth_rate);    
    /* update the instantaneous growth rate for the beginning of the next timestep */
    state->growth_rate = instantaneous_growth_rate;      
}

/*
 * START
 * Functions that handle each possible Gillespie event 
 *
 */
void transport_event(   Genotype *genotype,
                        CellState *state,
                        GillespieRates *rates,
                        float dt,
                        float t,
//                     long int *seed,
                        RngStream RS)
{
//    int gene_id,concurrent;
//    float x;
//    float endtime = t + dt + TTRANSLATION;
//
//    while(1)
//    {
//        x=RngStream_RandU01(RS)*rates->transport; 
//#if N_SIGNAL_TF==1
//        gene_id=0;
//#else
//        gene_id=1;
//#endif
//        /* choose gene product (mRNA) that gets transported to cytoplasm
//           based on weighting in transport[] array */
//        while (gene_id < genotype->ngenes-1 && x > 0.0) 
//        {
//            gene_id++;
//            x -= rates->transport_rate[gene_id];
//        }
//        if(rates->transport_rate[gene_id]>0.0)
//            break;
//    }
//    
//    (state->mRNA_nuclear_num[gene_id])--;   /* one less mRNAs in nucleus */
//    (state->mRNA_transl_cyto_num[gene_id])++;   /* it has just arrived in cytoplasm, ready to be translated */
//
//    /* add the endtime for translation */
//    concurrent=check_concurrence(   endtime,
//                                    state->mRNA_transl_init_time_end,
//                                    state->mRNA_transcr_time_end,
//                                    state->signalB_starts_end,
//                                    state->signalA_starts_end,
//                                    state->burn_in_growth_rate);
//    while(concurrent)
//    {
//        endtime+=0.0001;
//        concurrent=check_concurrence(   endtime,
//                                        state->mRNA_transl_init_time_end,
//                                        state->mRNA_transcr_time_end,
//                                        state->signalB_starts_end,
//                                        state->signalA_starts_end,
//                                        state->burn_in_growth_rate);        
//    }  
//    
//    add_fixed_event(gene_id, endtime, &(state->mRNA_transl_init_time_end), &(state->mRNA_transl_init_time_end_last));    
}

void mRNA_decay_event(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id;
    float x;
    int mRNA_id;

    while(1)/*in case of numerical error*/    
    {           
        x=RngStream_RandU01(RS)*rates->mRNAdecay;   
#if N_SIGNAL_TF==1
        gene_id=0;
#else
        gene_id=1;
#endif
        /* loop through mRNA products, to choose the mRNA with the
        proportionally higher decay rate */
        while (gene_id < genotype->ngenes-1 && x > 0.0) 
        {
            gene_id++;
            x-= rates->mRNAdecay_rate[gene_id];
        }
        /*rarely, numerical error picks up a gene that has not mRNA at all*/
        if((state->mRNA_cyto_num[gene_id]+state->mRNA_transl_cyto_num[gene_id])>=1)
            break;
    } 
    
    /* assume mRNA cytoplasm transport events equally likely */
    x = RngStream_RandInt(RS,1,state->mRNA_cyto_num[gene_id] + state->mRNA_transl_cyto_num[gene_id]);
    
    /* decay mRNA in cytoplasm */
    if (x <= state->mRNA_cyto_num[gene_id])
    {
        /* remove the mRNA from the cytoplasm count */
        (state->mRNA_cyto_num[gene_id])--;  
        
        change_mRNA_cytoplasm(gene_id, genotype, state, rates); 

    } 
    else 
    {
        /* decay mRNA in process of translating */       
        mRNA_id = RngStream_RandInt(RS,0,state->mRNA_transl_cyto_num[gene_id]-1);
        /* delete this fixed event: this mRNA will never be translated */
        delete_fixed_event(gene_id, mRNA_id, &(state->mRNA_transl_init_time_end), &(state->mRNA_transl_init_time_end_last));
       
        /* remove the mRNA from the count */
        (state->mRNA_transl_cyto_num[gene_id])--; 
    }
}

void histone_acteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x= RngStream_RandU01(RS)*rates->acetylation;
#if N_SIGNAL_TF==1
        gene_id=0;
#else
        gene_id=1;
#endif
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->acetylation_rate[gene_id];
        }
        if(rates->acetylation_rate[gene_id]>0.0)
            break;
    }
    /* set state: eject nucleosome, but there is no PIC yet */
    state->active[gene_id] = NO_NUC_NO_PIC;
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    int gene_id; 
    float x; 
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->deacetylation;
#if N_SIGNAL_TF==1
        gene_id=0;
#else
        gene_id=1;
#endif
        /* choose a particular gene and copy to change state */
        while(gene_id<genotype->ngenes-1 && x>0.0)
        {
            gene_id++;
            x-=rates->deacetylation_rate[gene_id];
        }
        if(rates->deacetylation_rate[gene_id]>0.0)
            break;
    }
    /* set state: nucleosome returns */
    state->active[gene_id] = NUC_NO_PIC;
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genotype, RngStream RS)
{
    float x; 
    int gene_id;
    while(1)
    {    
        x= RngStream_RandU01(RS)*rates->pic_assembly;
#if N_SIGNAL_TF==1
        gene_id=0;
#else
        gene_id=1;
#endif
        /* choose a particular gene and copy to change state */
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

void disassemble_PIC_event(Genotype *genotype, CellState *state,GillespieRates *rates, RngStream RS)
{
    int gene_id;
    float x;
    while(1)
    {
        x=RngStream_RandU01(RS)*rates->pic_disassembly;   
#if N_SIGNAL_TF==1
        gene_id=0;
#else
        gene_id=1;
#endif
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

void transcription_init_event(GillespieRates *rates, CellState *state, Genotype *genotype, float dt, float t, RngStream RS)
{
    int gene_id;  
    int x;
    float candidate_t;
    int concurrent;
    
#if N_SIGNAL_TF==1
    gene_id=0;
#else
    gene_id=1;
#endif
    
    x=RngStream_RandInt(RS,1,rates->transcript_init);
    while(gene_id<genotype->ngenes-1 && x>0)
    {
        gene_id++;
        x-=rates->transcript_init_rate[gene_id];
    }
     
    /* now that transcription of gene has been initiated, 
     * we add the timepoint at which the transcription ends, 
     * which is dt+time of transcription from now */
    candidate_t=t+dt+TTRANSCRIPTION;
    concurrent=check_concurrence(   candidate_t,
                                    state->mRNA_transl_init_time_end,
                                    state->mRNA_transcr_time_end,
                                    state->signalB_starts_end,
                                    state->signalA_starts_end,
                                    state->burn_in_growth_rate);
    while(concurrent)
    {
        candidate_t+=0.0001;
        concurrent=check_concurrence(   candidate_t,
                                        state->mRNA_transl_init_time_end,
                                        state->mRNA_transcr_time_end,
                                        state->signalB_starts_end,
                                        state->signalA_starts_end,
                                        state->burn_in_growth_rate);        
    }    
    add_fixed_event(gene_id, candidate_t, 
                        &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));

    /* increase the number mRNAs being transcribed */
    (state->mRNA_transcr_num[gene_id])++;                      
}
/*
 * END
 * Functions that handle each possible Gillespie event 
 */

/*copy genotype from the acestor to offsprings*/
void clone_cell_forward(Genotype *genotype_orig,                
                        Genotype *genotype_clone,
                        int clone_type)
{
    int i, j;
           
    if(clone_type!=MUT_KCONST) /* not a mutation in kinetic constants*/
    {
        for(i=0;i<NGENES;i++)
            genotype_clone->which_cluster[i]=-1;
        if(clone_type!=COPY_ALL)
        {
            for(i=0; i< genotype_orig->ngenes;i++)
            {
                genotype_clone->which_cluster[i]=genotype_orig->which_cluster[i];
                genotype_clone->recalc_TFBS[i]=0;
                if(genotype_clone->clone_info[i])  /* only copy to places that mutated*/                               
                {               
                    if(clone_type!= MUT_TFSEQ) /* if the mutation was not in binding sequence */
                        memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_orig->cisreg_seq[i][0],CISREG_LEN*sizeof(char));
                    
                    genotype_clone->recalc_TFBS[i]=1;
                }
            }
        }
        else
        {
            for(i=0; i< genotype_orig->ngenes;i++)
            {
                genotype_clone->which_cluster[i]=genotype_orig->which_cluster[i];
                genotype_clone->recalc_TFBS[i]=0;
                memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_orig->cisreg_seq[i][0],CISREG_LEN*sizeof(char));                    
                genotype_clone->recalc_TFBS[i]=1;                
            }
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
        /*then copy from orig*/
        i=0;
        while(genotype_orig->cisreg_cluster[i][0]!=-1)
        {
            j=0;
            while(genotype_orig->cisreg_cluster[i][j]!=-1)
            {
                genotype_clone->cisreg_cluster[i][j]=genotype_orig->cisreg_cluster[i][j];
                j++;
            }
            i++;
        }
    }
    
    if(clone_type!=MUT_CISSEQ) /* no change in gene numbers*/
    {
        for(i=0;i<genotype_clone->ngenes;i++)
        {
            genotype_clone->which_protein[i]=-1;
            genotype_clone->min_act_to_transc[i]=MAX_MODE+1;
        }        
        for(i=0; i < genotype_orig->ngenes; i++) 
        {            
            genotype_clone->mRNAdecay[i]=genotype_orig->mRNAdecay[i];
            genotype_clone->proteindecay[i]=genotype_orig->proteindecay[i];
            genotype_clone->translation[i]=genotype_orig->translation[i];            
            genotype_clone->pic_disassembly[i]=genotype_orig->pic_disassembly[i];
            genotype_clone->which_protein[i]=genotype_orig->which_protein[i];
            genotype_clone->min_act_to_transc[i]=genotype_orig->min_act_to_transc[i];
        }        
        for(i=0;i<genotype_clone->nproteins;i++)
        {            
            for(j=0;j<genotype_clone->protein_pool[i][0][0];j++)
                genotype_clone->protein_pool[i][1][j]=-1;
            genotype_clone->protein_pool[i][0][0]=0;            
        }
        for(i=0;i<genotype_orig->nproteins;i++)
        {            
            genotype_clone->protein_pool[i][0][0]=genotype_orig->protein_pool[i][0][0];            
            for(j=0;j<genotype_orig->protein_pool[i][0][0];j++)
                genotype_clone->protein_pool[i][1][j]=genotype_orig->protein_pool[i][1][j];                     
        }
    }
    
    /* since there is no tag to mark which binding seq is mutated, we copy all*/  
    for(i=0; i < genotype_orig->ntfgenes; i++) 
    {          
        for(j=0;j<TF_ELEMENT_LEN;j++)
        {    
            genotype_clone->tf_seq[i][j]=genotype_orig->tf_seq[i][j];
            genotype_clone->tf_seq_rc[i][j]=genotype_orig->tf_seq_rc[i][j];
        }
    }
    for(i=0;i<NPROTEINS;i++)
    {
        genotype_clone->activating[i]=genotype_orig->activating[i];
        genotype_clone->koff[i]=genotype_orig->koff[i];
    }
    
    /* these are easy update. Just do the update everytime*/
    genotype_clone->fitness=genotype_orig->fitness;
//    genotype_clone->var_fitness=genotype_orig->var_fitness;
    genotype_clone->ngenes=genotype_orig->ngenes;
    genotype_clone->ntfgenes=genotype_orig->ntfgenes;
    genotype_clone->nproteins=genotype_orig->nproteins;
    genotype_clone->N_act=genotype_orig->N_act;
    genotype_clone->N_rep=genotype_orig->N_rep; 
    genotype_clone->var_GR1=genotype_orig->var_GR1;
    genotype_clone->var_GR2=genotype_orig->var_GR2;
    genotype_clone->avg_GR1=genotype_orig->avg_GR1;
    genotype_clone->avg_GR2=genotype_orig->avg_GR2;
    genotype_clone->ncopies_under_penalty=genotype_orig->ncopies_under_penalty;    
    
    /* do not copy info for this gene, unless mutation changes it*/
    for(i=0;i<genotype_orig->ngenes;i++)    
        genotype_clone->clone_info[i]=0;
    
#if SET_BS_MANUALLY
    for(i=2;i<genotype_orig->ngenes;i++)
    {
        genotype_clone->binding_sites_num[i]=genotype_orig->binding_sites_num[i];
        genotype_clone->max_hindered_sites[i]=genotype_orig->max_hindered_sites[i];
        for(j=0;j<genotype_orig->binding_sites_num[i];j++)
        {
            genotype_clone->all_binding_sites[i][j].tf_id=genotype_orig->all_binding_sites[i][j].tf_id;
            genotype_clone->all_binding_sites[i][j].N_hindered=genotype_orig->all_binding_sites[i][j].N_hindered;
            genotype_clone->all_binding_sites[i][j].Koff=genotype_orig->all_binding_sites[i][j].Koff; 
            genotype_clone->all_binding_sites[i][j].mis_match=genotype_orig->all_binding_sites[i][j].mis_match;
        }
    }
#endif
}

/*
 * diagnostic function to dump a copy of the GillespieRates and some
 * of the CellState to the current error file
 */
void log_snapshot(Genotype *genotype, 
                  CellState *state,
                  GillespieRates *rates,                 
                  float x,
                  float t)
{
//  int i, p, nkon = 0;
//
//  LOG("snapshot at time=%g:\n x=%g, koff=%g = %d (tf_bound_num) * %g (koff/tf_bound_num)\n transport=%g\n decay=%g\n",
//      t, x, rates->koff, state->tf_bound_num, rates->koff/(float)state->tf_bound_num, 
//      rates->transport, rates->mRNAdecay);
//  LOG_NOFUNC(" rates->salphc=%g\n rates->max_salphc=%g rates->min_salphc=%g\n", rates->salphc, rates->max_salphc, rates->min_salphc);
//  LOG_NOFUNC(" konrate=%g\n", konrate);
//  LOG_NOFUNC(" pic_disassembly=%g\n kon=%g = %d * %g\n",
//             rates->pic_disassembly, rates->salphc+(konrate), kon_states->nkon, (rates->salphc+(konrate))/(float)kon_states->nkon);
//  
//  for (p=0; p < MAX_COPIES; p++) {
//    LOG_NOFUNC(" acetylation=%g (copy %d)\n deacetylation=%g (copy %d)\n PIC assembly=%g (copy %d)\n transcriptinit=%g (copy %d)\n",
//               (float)rates->acetylation_num[p]*ACETYLATE, p, (float)rates->deacetylation_num[p]*DEACETYLATE, p, 
//               (float)rates->pic_assembly_num[p]*PICASSEMBLY, p, (float)rates->transcript_init_num[p]*TRANSCRIPTINIT, p);
//  }
//  LOG_NOFUNC(" total rates=%g=%g+%g\n", rates->subtotal + (konrate), rates->subtotal, konrate);
//  LOG_NOFUNC(" total free=%d + total bound=%d + total hindered=%d = total sites=%d\n", 
//             kon_states->nkon, state->tf_bound_num, state->tf_hindered_num, genotype->binding_sites_num);
//  
//  LOG_NOFUNC("\n");
}

/*
 * run the model for a specified cell for a single timestep:
 */
void do_single_timestep(Genotype *genotype, 
                        CellState *state,                         
                        GillespieRates *rates, 
                        float *t,          
                        int maxbound2,
                        int maxbound3,                                      
                        char *env, 
                        float duration_env1,
                        float duration_env2,
                        int signalA_as_noise,
                        int signalA_mismatches,
                        RngStream RS,
			char *output,
                        char *error,
                        int mut_step,
                        Mutation *mut_record,
                        int *end_state,
                        int *did_burn_in) 
{    
    int event, UPDATE_WHAT;     
    float fixed_time; 
    float dt;
    float x;
    FILE *fp;

    x = expdev(RS);        /* draw random number */
    dt = x/rates->subtotal;
    if (dt < 0.0) 
    {	
        fp=fopen(error,"a+");
        fprintf(fp,"Step %d, %c %d %d %s %d %f error dt <0\n",
                mut_step,
                mut_record->mut_type,
                mut_record->pos_g,mut_record->pos_n,
                mut_record->nuc_diff,
                mut_record->kinetic_type,
                mut_record->kinetic_diff);
        fclose(fp);
        *end_state=0; /*use 0 to indicate abnormal behavior of the program.
                       *I am expecting rounding error to raise this flag.
                       *In case this flag is raised, quit the current replicate
                       *of growth and rerun a replicate.*/
        return;    
    }
    /* first check to see if a fixed event occurs in current t->dt window, or in tdevelopment if running for a fixed development time */
    fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;
    event = does_fixed_event_end(state->mRNA_transl_init_time_end,
                                 state->mRNA_transcr_time_end,
                                 state->signalB_starts_end,
                                 state->signalA_starts_end,
				 state->burn_in_growth_rate,
                                 fixed_time);
    while(event!=0)
    {           
        UPDATE_WHAT=do_fixed_event( genotype, state, rates, &dt, *t, event, 
                        duration_env1, duration_env2, env, 
                        signalA_as_noise, end_state, error, mut_step, mut_record);
        if(*end_state==0) 
            return; 
        *t += dt;  /* advance time by the dt */
#if CAUTIOUS
   float dt_copy,x_copy,rate_subtotal_copy,t_copy;  
   int event2;
        if(event==5)
        {
            *did_burn_in=1;
            state->cell_size_after_burn_in=state->cell_size;
        }          
        fixed_time=*t+dt;        
       
        event2=does_fixed_event_end(state->mRNA_transl_init_time_end,
                                 state->mRNA_transcr_time_end,
                                 state->signalB_starts_end,
                                 state->signalA_starts_end,
				 state->burn_in_growth_rate,
                                 fixed_time); 
        if(event2!=0)
        {
            fp=fopen(error,"a+");
            fprintf(fp,"step %d event %d overlaps with event %d\n",mut_step,event,event2);
            fclose(fp);
        }                       
	t_copy=*t; 
	dt_copy=dt;
	x_copy=x; 
	rate_subtotal_copy=rates->subtotal; 
#endif        
        x -= dt*rates->subtotal;  /* we've been running with rates->subtotal for dt, so substract it from x*/        
        calc_all_rates(genotype, state, rates, UPDATE_WHAT);  /* update rates->subtotal and re-compute a new dt */      
        dt = x/rates->subtotal;
	/*deal with rounding error*/
	if(dt<0.0)
	{  	
#if CAUTIOUS
		//fp=fopen(error,"a+");
		//fprintf(fp,"step %d, rounding error in dt\n",mut_step);
        	//fprintf(fp,"Step %d, %c %d %d %s %d %f, t=%f, rounding error in dt\n",
		//	mut_step,
		//	mut_record->mut_type,
		//	mut_record->pos_g,
		//	mut_record->pos_n,
		//	mut_record->nuc_diff,
		//	mut_record->kinetic_type,
		//	mut_record->kinetic_diff,
		//	*t);
		//fprintf(fp,"original dt=%f, original x=%f, original Rsubtotal=%f, new x=%f, new dt=%f\n",
		//	dt_copy,
		//	x_copy,
		//	rate_subtotal_copy,
		//	x,
		//	dt);	
        	//fclose(fp);
#endif 
		dt=0.0001; 
	}
        fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment; 
        /* check to see there aren't more fixed events to do */
        event = does_fixed_event_end(state->mRNA_transl_init_time_end, 
                                     state->mRNA_transcr_time_end,
                                     state->signalB_starts_end,
                                     state->signalA_starts_end,
				     state->burn_in_growth_rate,
                                     fixed_time);                                    
    } 
  /* no remaining fixed events to do in dt, now do stochastic events */  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (*t+dt < tdevelopment)
    { 
      /* if the total rates falls below zero, we do an emergency recalibration of cell */
        if (!(rates->subtotal > 0.0)) 
        {
            calc_all_rates(genotype, state, rates, INITIALIZATION);
            /* if this still results in either zero or negative total rates,
               this most likely due the cell being "dead" no TFs bound, no
               activity etc.  We mark cell as "dead" in this case, and
               remove from queue. */
            if (!(rates->subtotal > 0.0)) 
            {  
                fp=fopen(error,"a+");
		fprintf(fp,"Step %d, %c %d %d %s %d %f unresponding cell\n",
                mut_step,
                mut_record->mut_type,
                mut_record->pos_g,mut_record->pos_n,
                mut_record->nuc_diff,
                mut_record->kinetic_type,
                mut_record->kinetic_diff);
                fclose(fp);
                *end_state=-1;
                return;        
            }
        } 
        update_protein_number_cell_size(genotype, state, rates, dt, *t, *env, end_state, error, mut_step, mut_record);        
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, *t, RS, end_state, error, mut_step, mut_record);
        if(*end_state==0)
            return;               
        calc_all_rates(genotype,state,rates,UPDATE_WHAT);
        /* Gillespie step: advance time to next event at dt */
        *t += dt;
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - *t;
        /* final update of protein concentration */
        update_protein_number_cell_size(genotype, state, rates, dt, *t, *env, end_state, error, mut_step, mut_record); 
        if(*end_state==0)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        *t = tdevelopment;
    }
}
/* while there are either transcription or translation events
     occuring in current t->dt window */
int do_fixed_event(Genotype *genotype, 
                    CellState *state, 
                    GillespieRates *rates, 
                    float *dt,
                    float t,        
                    int event, 
                    float duration_env1,
                    float duration_env2,
                    char *env,
                    int signalA_as_noise,                
                    int *end_state,
                    char *error,
                    int mut_step,
                    Mutation *mut_record)
{     
    int return_value;
    return_value=0;
    switch (event) 
    {
        case 1:     /* if a transcription event ends */
            end_transcription(dt, t, state, rates, genotype,end_state,error,mut_step,mut_record,env);                 
            break;
        case 2:     /* if a translation initialization event ends */ 
            end_translation_init(genotype, state, rates, dt, t,end_state,error,mut_step,mut_record,env);  
            return_value=UPDATE_PACT;
            break;
        case 3:     /* environment B starts*/ 
            *dt = state->signalB_starts_end->time - t;     
            update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, error, mut_step, mut_record); 
            delete_fixed_event_start(&(state->signalB_starts_end),&(state->signalB_starts_end_last));
            *env='B';
            state->protein_number[1]=basal_signal_strength;
            break;
        case 4:     /*environment A starts*/
            *dt = state->signalA_starts_end->time - t;   
            update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, error, mut_step, mut_record);  
            delete_fixed_event_start(&(state->signalA_starts_end),&(state->signalA_starts_end_last));
            if(!signalA_as_noise)
            {
                state->protein_number[1]=signal_strength;                
                *env='A';
            }
            else                 
                state->protein_number[1]=signal_strength; 
            break;	
	case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-t;     
            update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, error, mut_step, mut_record);
            state->cell_size_after_burn_in=state->cell_size;           
            delete_fixed_event_start(&(state->burn_in_growth_rate),&(state->burn_in_growth_rate_last));
            break;
#if CAUTIOUS
        FILE *fp;
        default:
            fp=fopen(error,"a+");
            fprintf(fp,"Step %d, %c %d %d %s %d %f error in do_fixed_event\n",mut_step,mut_record->mut_type,mut_record->pos_g,mut_record->pos_n,mut_record->nuc_diff,mut_record->kinetic_type,mut_record->kinetic_diff);
            fclose(fp);
            *end_state=0;      
            return;
#endif
    }  
    return return_value;
}

int do_Gillespie_event(Genotype *genotype,
                        CellState *state,
                        GillespieRates *rates,
                        float dt,
                        float t,
                        RngStream RS,
                        int *end_state,
                        char *error,
                        int mut_step,
                        Mutation *mut_record)
{
    float x;
    int return_value;
    return_value=0;
//    FILE *fp;  
    
    while(1)
    {    
        x=RngStream_RandU01(RS)*(rates->subtotal);
//        if (x <= rates->transport)   /* transportation event */ 
//        {  
//            transport_event(genotype, state, rates, dt, t, RS);            
//            break;
//        } 
//        else 
//        {   
//            x -= rates->transport;           
            if (x <= rates->mRNAdecay)  /*STOCHASTIC EVENT: an mRNA decay event */
            {    
                mRNA_decay_event(rates, state, genotype, RS);
                return_value=UPDATE_PACT;
                break;
            } 
            else 
            {               
                x -= rates->mRNAdecay;             
                if (x <= rates->pic_disassembly) /* STOCHASTIC EVENT: PIC disassembly*/
                {                   
                    disassemble_PIC_event(genotype, state, rates, RS);
                    break;
                } 
                else 
                {  
                    x -= rates->pic_disassembly;  
                    if (x <= rates->acetylation)  /* acetylation*/
                    {    
                        histone_acteylation_event(rates, state, genotype, RS);
                        break;
                    } 
                    else 
                    { 
                        x-= rates->acetylation;                  
                        if (x <= rates->deacetylation)/* STOCHASTIC EVENT: histone deacetylation */ 
                        {   
                            histone_deacteylation_event(rates, state, genotype, RS);
                            break;
                        } 
                        else 
                        {
                            x -= rates->deacetylation;                        
                            if (x <= rates->pic_assembly)/* STOCHASTIC EVENT: PIC assembly*/
                            {  
                                assemble_PIC_event(rates, state, genotype, RS); 
                                break;
                            } 
                            else 
                            {
                                x -= rates->pic_assembly;                            
                                if (x <= (float)rates->transcript_init * TRANSCRIPTINIT) /* STOCHASTIC EVENT: transcription initiation */
                                {   
                                    transcription_init_event(rates, state, genotype, dt, t, RS);
                                    break;
                                }
                                else // the rest of the rate is to mandatorily update Pact
                                    break;
//                                else 
//                                {
//                                    /*
//                                     * FALLBACK: shouldn't get here, previous
//                                     * events should be exhaustive
//                                     */
//                                   // fp=fopen(error,"a+");
//                                   // fprintf(fp,"Step %d, %c %d %d %s %d %f error in do_Gillespie_event\n",
//                                   //         mut_step,
//                                   //         mut_record->mut_type,
//                                   //         mut_record->pos_g,
//                                   //         mut_record->pos_n,
//                                   //         mut_record->nuc_diff,
//                                   //         mut_record->kinetic_type,
//                                    //        mut_record->kinetic_diff);
//                                    //fprintf(fp,"x=%f, x_copy=%f, subtotal=%f, transport=%f, mRNAdecay=%f, pic_disassembly=%f, acetylation=%f, deacetylation=%f, pic_assembly=%f, trancrip_init=%0f\n",
//                                    //        x,
//                                    //	x_copy,
//                                    //	rates->subtotal,
//                                    //        rates->transport,
//                                    //        rates->mRNAdecay,
//                                    //        rates->pic_disassembly,
//                                    //        rates->acetylation,
//                                    //        rates->deacetylation,
//                                    //        rates->pic_assembly,
//                                    //        (float)rates->transcript_init*TRANSCRIPTINIT);
//                                    //fclose(fp);
//                                    *end_state=0;
//                                    return;                                
//                                }
                            }
                        }
                    }
                }       
            }
//        }
    }
    return return_value;
}

void free_fixedevent(CellState *state)
{
    FixedEvent *temp1, *temp2;

    temp1=state->signalA_starts_end;

    while(temp1){		
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signalA_starts_end=NULL;
    state->signalA_starts_end_last=NULL;
    
    temp1=state->signalB_starts_end;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signalB_starts_end=NULL;
    state->signalB_starts_end_last=NULL;
    
    temp1=state->mRNA_transcr_time_end;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;           
            free(temp2);	
    }
    state->mRNA_transcr_time_end=NULL;
    state->mRNA_transcr_time_end_last=NULL;
    
    temp1=state->mRNA_transl_init_time_end;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->mRNA_transl_init_time_end=NULL;
    state->mRNA_transl_init_time_end_last=NULL;

    temp1=state->burn_in_growth_rate;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->burn_in_growth_rate=NULL;
    state->burn_in_growth_rate_last=NULL;
}

void calc_avg_growth_rate(  Genotype *genotype,
                            int init_mRNA[NGENES],
                            float init_protein_number[NPROTEINS],
                            float maxbound2,
                            float maxbound3,
                            RngStream RS_parallel[N_THREADS],
                            char *output,
                            char *error,
                            int mut_step,
                            Mutation *mut_record)       
{       
    float GR1[N_replicates],GR2[N_replicates]; 
#if !RdcPdup
    int k;
    genotype->ncopies_under_penalty=0;
    for(k=N_SIGNAL_TF;k<genotype->nproteins;k++)    
        genotype->ncopies_under_penalty+=(genotype->protein_pool[k][0][0]<=MAX_COPIES)?0:(genotype->protein_pool[k][0][0]-MAX_COPIES);   
#endif    
        
    omp_set_num_threads(N_THREADS);
    #pragma omp parallel
    {
        int ID=omp_get_thread_num();
//        int ID=0;
        int i,j,k;
        int N_replicates_per_thread=N_replicates/N_THREADS;
        int end_state;
        char env;
        Genotype genotype_offspring;
        CellState state_offspring;
        GillespieRates rate_offspring;
        int init_mRNA_offspring[NGENES]; 
        float init_protein_number_offspring[NGENES];
        float gr1[N_replicates_per_thread],gr2[N_replicates_per_thread];        
        float t;
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];       
        int did_burn_in;
	
        initialize_cache(&genotype_offspring);  
        
        for(i=0;i<N_replicates_per_thread;i++)
        {
            gr1[i]=0.0;
            gr2[i]=0.0;
        }
        #pragma omp critical
        {            
            genotype_offspring.ngenes=genotype->ngenes;
            genotype_offspring.ntfgenes=genotype->ntfgenes;
            genotype_offspring.nproteins=genotype->nproteins;
            clone_cell_forward(genotype, &genotype_offspring, COPY_ALL);
            for(j=0; j < NGENES; j++) 
            {  
                init_mRNA_offspring[j] = init_mRNA[j];
                init_protein_number_offspring[j] = init_protein_number[j];
            }
#if SET_BS_MANUALLY
            for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            {
                genotype_offspring.binding_sites_num[i]=genotype->binding_sites_num[i];
                genotype_offspring.max_hindered_sites[i]=genotype->max_hindered_sites[i];
                for(j=0;j<genotype->binding_sites_num[i];j++)
                {
                    genotype_offspring.all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i][j].tf_id;
                    genotype_offspring.all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i][j].N_hindered;
                    genotype_offspring.all_binding_sites[i][j].Koff=genotype->all_binding_sites[i][j].Koff;                    
                }
            }   
#endif
        }      
        calc_all_binding_sites(&genotype_offspring); 
        
        /* now calc growth rate under two environments*/
        for(i=0;i<N_replicates_per_thread;i++) /* env 1, usually a constant signal that matches env*/
        {	 
            env=init_env1;
            end_state=1;
            did_burn_in=0;
            /*the first gene is a only a signal that takes the form of a tf*/
            for(j=N_SIGNAL_TF; j < genotype_offspring.ngenes; j++)             
                mRNA[j] = init_mRNA_offspring[j];                       
            for(j=N_SIGNAL_TF; j<genotype_offspring.nproteins;j++)
            {
                for(k=0;k<genotype_offspring.protein_pool[j][0][0];k++)
                    protein[genotype_offspring.protein_pool[j][1][k]]=(float)init_protein_number_offspring[j]/genotype_offspring.protein_pool[j][0][0];                
            }
            initialize_cell(&state_offspring, genotype_offspring.ngenes, genotype_offspring.nproteins,
                            genotype_offspring.protein_pool,genotype_offspring.mRNAdecay, mRNA, protein, 
                            0, RS_parallel[ID]);
            set_env(&state_offspring,env,env1_t_signalA,env1_t_signalB);
            calc_all_rates(&genotype_offspring, &state_offspring, &rate_offspring, INITIALIZATION);
            t = 0.0;

            while(t<tdevelopment && end_state==1)
            {
               do_single_timestep(&genotype_offspring, &state_offspring, &rate_offspring, 
                                    &t,                                        
                                    maxbound2,
                                    maxbound3,                                                                                         
                                    &env,
                                    env1_t_signalA,
                                    env1_t_signalB,
                                    env1_signalA_as_noise,
                                    env1_signalA_mismatches,
                                    RS_parallel[ID],
				    output,
                                    error,
                                    mut_step,
                                    mut_record,
                                    &end_state,
                                    &did_burn_in);
            } 
            if(end_state==1)
	    {	
                gr1[i]=log(state_offspring.cell_size/state_offspring.cell_size_after_burn_in)/(tdevelopment-duration_of_burn_in_growth_rate);
#if CAUTIOUS
                FILE *fp;
		if(gr1[i]<0.0)
		{
			fp=fopen(error,"a+");
			fprintf(fp,"negative growth rate at %d\n",mut_step);
			fprintf(fp,"cell size after burn-in=%.20f cell size=%.20f\n",state_offspring.cell_size_after_burn_in,state_offspring.cell_size);
			fprintf(fp,"tdevelopment=%f duration of burn-in=%f\n",tdevelopment,duration_of_burn_in_growth_rate);
			fclose(fp);
		}
                if(did_burn_in==0)
                {
                    fp=fopen(error,"a+");
                    fprintf(fp,"error at step %d, thread %d, replication %d, env %c, burn-in did not engage\n",mut_step,ID,i,env);
                    fclose(fp);
                }
#endif
	    }
            else
	    { 
		if(end_state==-1) // unresponding cell. There is no dead cell in this model. 
                    gr1[i]=0.0;
            	else   /* if do_single_timestep throws out an error, ignore it and re-run*/
                    i--;
            }		
            
            free_fixedevent(&state_offspring);           
        }
        
        /*******************env2*********************/
        for(i=0;i<N_replicates_per_thread;i++) 
        {	 
            env=init_env2;  
            end_state=1;
            did_burn_in=0;
            /*the first gene is a only a signal that takes the form of a tf*/
            for(j=N_SIGNAL_TF; j < genotype_offspring.ngenes; j++)             
                mRNA[j] = init_mRNA_offspring[j];                       
            for(j=N_SIGNAL_TF; j<genotype_offspring.nproteins;j++)
            {
                for(k=0;k<genotype_offspring.protein_pool[j][0][0];k++)
                    protein[genotype_offspring.protein_pool[j][1][k]]=(float)init_protein_number_offspring[j]/genotype_offspring.protein_pool[j][0][0];                
            }
            initialize_cell(&state_offspring, genotype_offspring.ngenes, genotype_offspring.nproteins,
                            genotype_offspring.protein_pool,genotype_offspring.mRNAdecay, mRNA, protein, 
                            0, RS_parallel[ID]);
            set_env(&state_offspring,env,env2_t_signalA,env2_t_signalB);
            calc_all_rates(&genotype_offspring, &state_offspring, &rate_offspring, INITIALIZATION);	
            t = 0.0;

            while(t<tdevelopment && end_state==1)
            {
                do_single_timestep(&genotype_offspring, &state_offspring, &rate_offspring, 
                                    &t,                                        
                                    maxbound2,
                                    maxbound3,                                                                                         
                                    &env,
                                    env2_t_signalA,
                                    env2_t_signalB,
                                    env2_signalA_as_noise,
                                    env2_signalA_mismatches,                                  
                                    RS_parallel[ID],
				    output,
                                    error,
                                    mut_step,
                                    mut_record,
                                    &end_state,
                                    &did_burn_in);
            } 
            if(end_state==1)
            {
		gr2[i]=log(state_offspring.cell_size/state_offspring.cell_size_after_burn_in)/(tdevelopment-duration_of_burn_in_growth_rate);
#if CAUTIOUS
                FILE *fp;
                if(gr2[i]<0.0)
		{
                    fp=fopen(error,"a+");
                    fprintf(fp,"negative growth rate at %d\n",mut_step);
                    fprintf(fp,"cell size after burn-in=%.20f cell size=%.20f\n",state_offspring.cell_size_after_burn_in,state_offspring.cell_size);		
                    fclose(fp);
		}
                if(did_burn_in==0)
                {
                    fp=fopen(error,"a+");
                    fprintf(fp,"error at step %d, thread %d, replication %d, env %c, burn-in did not engage\n",mut_step,ID,i,env);
                    fclose(fp);
                }
#endif
	    }
            else
	    {	
                if(end_state==-1)
                    gr2[i]=0.0;
            	else
                    i--;
	    }
            
            free_fixedevent(&state_offspring);            
        }    
        
        for(j=0;j<NGENES;j++)
            free(genotype_offspring.all_binding_sites[j]);
        #pragma omp critical
        {
            j=0;
            for(i=ID*N_replicates_per_thread;i<(ID+1)*N_replicates_per_thread;i++)
            {
                GR1[i]=gr1[j];
                GR2[i]=gr2[j];
                j++;
            }
        }        
    }
    
    int i;
    float avg_GR1,avg_GR2,var_GR1,var_GR2;  
  
    avg_GR1=0.0;
    avg_GR2=0.0;
    var_GR1=0.0;
    var_GR2=0.0;   
    
    for(i=0;i<N_replicates;i++)
    {
        avg_GR1+=GR1[i];
        avg_GR2+=GR2[i];        
    }    
    avg_GR1=avg_GR1/(float)N_replicates;
    avg_GR2=avg_GR2/(float)N_replicates;
    
    for(i=0;i<N_replicates;i++)
    {
        var_GR1+=pow(GR1[i]-avg_GR1,2.0);
        var_GR2+=pow(GR2[i]-avg_GR2,2.0);
    }
    var_GR1=var_GR1/(float)(N_replicates*(N_replicates-1));
    var_GR2=var_GR2/(float)(N_replicates*(N_replicates-1));     
    genotype->avg_GR1=avg_GR1;
    genotype->avg_GR2=avg_GR2;
    genotype->var_GR1=var_GR1;
    genotype->var_GR2=var_GR2;
    genotype->fitness=(env1_occurence*avg_GR1+env2_occurence*avg_GR2)/(env1_occurence+env2_occurence);
}

#if PLOTTING
int does_fixed_event_end_plotting(  FixedEvent *mRNA_transl_init_time_end,
                                    FixedEvent *mRNA_transcr_time_end,
                                    FixedEvent *signalB_starts_end,
                                    FixedEvent *signalA_starts_end,	
                                    FixedEvent *burn_in_growth_rate,
                                    FixedEvent *sampling_point_end,
                                    float t) 
{
    int retval=0;
    float t1;
    float t2;
    float t3;
    float t4;
    float t5;    
    float t6;

    if(mRNA_transcr_time_end == NULL && mRNA_transl_init_time_end==NULL && signalB_starts_end == NULL && signalA_starts_end == NULL && burn_in_growth_rate == NULL && sampling_point_end == NULL)
    {
        retval =0;
    }
    else
    {
        t1 = mRNA_transcr_time_end ? mRNA_transcr_time_end->time : TIME_INFINITY;
        t2 = mRNA_transl_init_time_end ? mRNA_transl_init_time_end->time : TIME_INFINITY;
        t3 = signalB_starts_end ? signalB_starts_end->time : TIME_INFINITY;
        t4 = signalA_starts_end ? signalA_starts_end->time : TIME_INFINITY;
	t5 = burn_in_growth_rate ? burn_in_growth_rate->time : TIME_INFINITY;
        t6 = sampling_point_end ? sampling_point_end->time : TIME_INFINITY;

        if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4) && (t1<=t5) && (t1<=t6))
	{
            retval = 1;	
        }
        else
        {
            if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4) && (t2<=t5) && (t2<=t6)) 
            {
                retval = 2;
            }
            else
            {
                if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4) && (t3<=t5) && (t3<=t6)) 
                {
                    retval = 3;
                }
                else
                {
                    if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3) && (t4<=t5) && (t4<=t6)) 
                    {
                        retval = 4;
                    }
                    else
                    {
			if((t5 <= t1) && (t5 <= t) && (t5 <= t2) && (t5 <= t3) && (t5<=t4) && (t5<=t6))
			{
                            retval = 5;
			}
			else
                        {
                            if((t6 <= t1) && (t6 <= t) && (t6 <= t2) && (t6 <= t3) && (t6<=t4) && (t6<=t5))
                                retval=6;
                            else
                                retval = 0;
                        }				
                    }
                }
            }
        }
    }
    return retval;
}

int do_fixed_event_plotting(    Genotype *genotype, 
                                CellState *state, 
                                GillespieRates *rates, 
                                float *dt,
                                float t,        
                                int event, 
                                float duration_env1,
                                float duration_env2,
                                char *env,
                                int signalA_as_noise,                              
                                int *end_state)
{  
    FILE *fp;
    int return_value;
    return_value=0;
    switch (event) 
    {        
        case 1:     /* if a transcription event ends */
            end_transcription(dt, t, state, rates, genotype,end_state,NULL,0,NULL,env);                 
            break;
        case 2:     /* if a translation event ends */ 
            end_translation_init(genotype, state, rates, dt, t,end_state,NULL,0,NULL,env);   
            return_value=UPDATE_PACT;
            break;
        case 3:     /* environment B starts*/ 
            *dt = state->signalB_starts_end->time - t;     
            update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, NULL, 0, NULL);  
            delete_fixed_event_start(&(state->signalB_starts_end),&(state->signalB_starts_end_last)); 
            state->protein_number[1]=basal_signal_strength;
            *env='B';
            break;
        case 4:     /*environment A starts. See "Environments" for detail*/
            *dt = state->signalA_starts_end->time - t;   
            update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, NULL, 0, NULL);  
            delete_fixed_event_start(&(state->signalA_starts_end),&(state->signalA_starts_end_last));
            if(!signalA_as_noise)
            {
                state->protein_number[1]=signal_strength;                
                *env='A';
            }
            else                 
                state->protein_number[1]=signal_strength;
            break;	
	case 5: /* finishing burn-in growth rate*/
            *dt=duration_of_burn_in_growth_rate-t;     
            update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, NULL, 0, NULL);
            delete_fixed_event_start(&(state->burn_in_growth_rate),&(state->burn_in_growth_rate_last));
            break;
        case 6:
            *dt=state->sampling_point_end->time-t;
            update_protein_number_cell_size(genotype, state, rates, *dt, t, *env, end_state, NULL, 0, NULL);
            delete_fixed_event_start(&(state->sampling_point_end),&(state->sampling_point_end_last));
            break;
        default:          
            *end_state=0;      
            break;
    }    
    return return_value;
}

void do_single_timestep_plotting(   Genotype *genotype, 
                                    CellState *state,                         
                                    GillespieRates *rates, 
                                    float *t,          
                                    int maxbound2,
                                    int maxbound3,                                      
                                    char *env, 
                                    float duration_env1,
                                    float duration_env2,
                                    int signalA_as_noise,                                                                       
                                    float (*phenotype)[149],                                    
                                    float fitness[149],                                    
                                    RngStream RS,                                    
                                    int *timepoint,                                   
                                    int *end_state) 
{    
    int event,i,mut_step,UPDATE_WHAT;     /* boolean to keep track of whether FixedEvent has ended */   
    float fixed_time; 
    float dt;
    float x;    
    mut_step=0;

    x = expdev(RS);        /* draw random number */
    dt = x/rates->subtotal;
    /* first check to see if a fixed event occurs in current t->dt window, or in tdevelopment if running for a fixed development time */
    fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;
    event = does_fixed_event_end_plotting(  state->mRNA_transl_init_time_end,
                                            state->mRNA_transcr_time_end,
                                            state->signalB_starts_end,
                                            state->signalA_starts_end,
                                            state->burn_in_growth_rate,
                                            state->sampling_point_end,
                                            fixed_time);
    while(event!=0)
    {           
        UPDATE_WHAT=do_fixed_event_plotting(genotype, state, rates, &dt, *t, 
                                            event, duration_env1, duration_env2, env,
                                            signalA_as_noise, end_state);
        if(*end_state==0)
            return;        
        if(event==6)
        {    
            for(i=0;i<genotype->nproteins;i++)
                phenotype[i][*timepoint]=state->protein_number[i];
            
            fitness[*timepoint]=state->growth_rate;          
            (*timepoint)++;    
        }              
        fixed_time=*t+dt; 
        *t += dt;                  /* advance time by the dt */	    
        x -= dt*rates->subtotal;  /* we've been running with rates->subtotal for dt, so substract it from x*/        
        calc_all_rates(genotype, state, rates, UPDATE_WHAT);  /* update rates->subtotal and re-compute a new dt */ 
        dt = x/rates->subtotal;
	/*deal with rounding error*/
	if(dt<0.0)	
            dt=0.0001;         
        /* check to see there aren't more fixed events to do */
        fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;         
        event = does_fixed_event_end_plotting(  state->mRNA_transl_init_time_end, 
                                                state->mRNA_transcr_time_end,
                                                state->signalB_starts_end,
                                                state->signalA_starts_end,
                                                state->burn_in_growth_rate,
                                                state->sampling_point_end,
                                                fixed_time);                                    
    } 
  /* no remaining fixed events to do in dt, now do stochastic events */  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */          
    if (*t+dt < tdevelopment)
    { 
        update_protein_number_cell_size(genotype, state, rates, dt, *t, *env, end_state, NULL, 0, NULL);  
        if(*end_state==0)
            return; 
        UPDATE_WHAT=do_Gillespie_event(genotype, state, rates, dt, *t, RS, end_state, NULL, 0, NULL);
        if(*end_state==0)
            return;        
        calc_all_rates(genotype,state,rates,UPDATE_WHAT);       
        /* Gillespie step: advance time to next event at dt */
        *t += dt;
    } 
    else 
    { 
        /* do remaining dt */
        dt = tdevelopment - *t;
        /* final update of protein concentration */
        update_protein_number_cell_size(genotype, state, rates, dt, *t, *env, end_state, NULL, 0, NULL); 
       
        if(*end_state==0)
            return;         
        /* advance to end of development (this exits the outer while loop) */
        *t = tdevelopment;
    }
}

void calc_avg_growth_rate_plotting( Genotype *genotype,
                                    int init_mRNA[NGENES],
                                    float init_protein_number[NPROTEINS],
                                    float maxbound2,
                                    float maxbound3,
                                    RngStream RS[N_THREADS])
{     
    float phenotypeA[N_replicates][genotype->ngenes][149];
    float phenotypeB[N_replicates][genotype->ngenes][149];
    float fitnessA[N_replicates][149];
    float fitnessB[N_replicates][149];
    int l,m,n;
    char filename1[32],filename2[32];
    FILE *fp1,*fp2;
    int N_update=0;
    
    for(l=0;l<N_replicates;l++)
    {
        for(m=0;m<genotype->ngenes;m++)
        {
            for(n=0;n<149;n++)
            {
                phenotypeA[l][m][n]=0.0;
                phenotypeB[l][m][n]=0.0;               
            }
        }
    }    
    for(l=0;l<N_replicates;l++)
    {   
        for(n=0;n<149;n++)
            {
                fitnessA[l][n]=0.0;
                fitnessB[l][n]=0.0;               
            }        
    }
    
    omp_set_num_threads(N_THREADS);
   
    #pragma omp parallel
    {
         printf("%d",omp_get_num_threads());
    
        int ID=omp_get_thread_num();
//        int ID=0;
        int i,j,timepoint,k;    
        int end_state;  
        int N_replicates_per_thread=N_replicates/N_THREADS;
//        int N_replicates_per_thread=100;
        char env;
        Genotype genotype_offspring;
        CellState state_offspring;
        GillespieRates rate_offspring;
        float init_mRNA_offspring[NGENES]; 
        float init_protein_number_offspring[NGENES];
        float t;
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];       

        initialize_cache(&genotype_offspring); 
        clone_cell_forward(genotype, &genotype_offspring, COPY_ALL);
        
        #pragma omp critical
        {            
            genotype_offspring.ngenes=genotype->ngenes;
            genotype_offspring.ntfgenes=genotype->ntfgenes;
            genotype_offspring.nproteins=genotype->nproteins;
            clone_cell_forward(genotype, &genotype_offspring, COPY_ALL);
            for(j=0; j < NGENES; j++) 
            {  
                init_mRNA_offspring[j] = init_mRNA[j];
                init_protein_number_offspring[j] = init_protein_number[j];
            }
        #if SET_BS_MANUALLY
            for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            {
                genotype_offspring.binding_sites_num[i]=genotype->binding_sites_num[i];
                genotype_offspring.max_hindered_sites[i]=genotype->max_hindered_sites[i];
                for(j=0;j<genotype->binding_sites_num[i];j++)
                {
                    genotype_offspring.all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i][j].tf_id;
                    genotype_offspring.all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i][j].N_hindered;
                    genotype_offspring.all_binding_sites[i][j].Koff=genotype->all_binding_sites[i][j].Koff;                    
                }
            }   
        #endif
        }
        
    #if !SET_BS_MANUALLY        
        calc_all_binding_sites(&genotype_offspring); 
    #endif      
      
        for(i=0;i<N_replicates_per_thread;i++)        
        {             
            env=init_env1 ; 
                
            /*the first one or two genes are a only a signal that takes the form of a tf*/
            for(j=N_SIGNAL_TF; j < genotype_offspring.ngenes; j++)             
                mRNA[j] = init_mRNA_offspring[j];                       
            for(j=N_SIGNAL_TF; j<genotype_offspring.nproteins;j++)
            {
                for(k=0;k<genotype_offspring.protein_pool[j][0][0];k++)
                    protein[genotype_offspring.protein_pool[j][1][k]]=(float)init_protein_number_offspring[j]/genotype_offspring.protein_pool[j][0][0];                
            }
            initialize_cell(&state_offspring, genotype_offspring.ngenes, genotype_offspring.nproteins,
                            genotype_offspring.protein_pool,genotype_offspring.mRNAdecay, mRNA, protein, 
                            PLOTTING, RS[ID]);
            set_env(&state_offspring,env,env1_t_signalA,env1_t_signalB);
            calc_all_rates(&genotype_offspring, &state_offspring, &rate_offspring, INITIALIZATION);
            
            t = 0.0;
            end_state=1;          
            timepoint=0;
            while(t<tdevelopment && end_state==1)
            {    
                do_single_timestep_plotting(&genotype_offspring, &state_offspring, &rate_offspring, 
                                            &t,                                        
                                            maxbound2,
                                            maxbound3,                                                                                         
                                            &env,
                                            env1_t_signalA,
                                            env1_t_signalB,
                                            env1_signalA_as_noise,                                                                     
                                            phenotypeA[ID*N_replicates_per_thread+i],                                           
                                            fitnessA[ID*N_replicates_per_thread+i],  
                                            RS[ID],                                           
                                            &timepoint,                                           
                                            &end_state);        
            }
            if(end_state==0)
                i--;            
            free_fixedevent(&state_offspring);
        }

          
        for(i=0;i<N_replicates_per_thread;i++)        
        {  
            
            env=init_env2 ;        
            /*the first gene is a only a signal that takes the form of a tf*/
            for(j=N_SIGNAL_TF; j < genotype_offspring.ngenes; j++)             
                mRNA[j] = init_mRNA_offspring[j];                       
            for(j=N_SIGNAL_TF; j<genotype_offspring.nproteins;j++)
            {
                for(k=0;k<genotype_offspring.protein_pool[j][0][0];k++)
                    protein[genotype_offspring.protein_pool[j][1][k]]=(float)init_protein_number_offspring[j]/genotype_offspring.protein_pool[j][0][0];                
            }
            initialize_cell(&state_offspring, genotype_offspring.ngenes, genotype_offspring.nproteins,
                            genotype_offspring.protein_pool,genotype_offspring.mRNAdecay, mRNA, protein, 
                            PLOTTING, RS[ID]);
            set_env(&state_offspring,env,env2_t_signalA,env2_t_signalB);
            calc_all_rates(&genotype_offspring, &state_offspring, &rate_offspring, INITIALIZATION);
            t = 0.0;
            end_state=1;           
            timepoint=0;
            while(t<tdevelopment && end_state==1)
            {
                do_single_timestep_plotting(&genotype_offspring, &state_offspring, &rate_offspring, 
                                            &t,                                        
                                            maxbound2,
                                            maxbound3,                                                                                         
                                            &env,
                                            env2_t_signalA,
                                            env2_t_signalB,
                                            env2_signalA_as_noise,                                            
                                            phenotypeB[ID*N_replicates_per_thread+i],
                                            fitnessB[ID*N_replicates_per_thread+i],                                            
                                            RS[ID],                                          
                                            &timepoint,                                            
                                            &end_state);
            }
            if(end_state==0)
                i--;            
            free_fixedevent(&state_offspring);
        }        
    }//end of parallel section   
    
    for(l=0;l<genotype->nproteins;l++)
    {
        snprintf(filename1,sizeof(char)*32,"phenotypeA_%i",l);
        snprintf(filename2,sizeof(char)*32,"phenotypeB_%i",l);
        fp1=fopen(filename1,"w");
        fp2=fopen(filename2,"w");
        for(m=0;m<149;m++)
        {
            for(n=0;n<N_replicates;n++)
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
    fp1=fopen("fitnessA","w");
    fp2=fopen("fitnessB","w");    
    for(m=0;m<149;m++)
    {
        for(n=0;n<N_replicates;n++)
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
#endif

/************* begin of mutation functions ***************/

void mut_substitution(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int pos_n, pos_g;
    char *Genome,nucleotide;
    float random;
    Genome= &genotype->cisreg_seq[0][0];
    #if !SIMPLE_SUBSTITUTION    
    while(1)
    {
        pos_n=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
        random=RngStream_RandU01(RS);
        if (Genome[pos_n]=='a'||Genome[pos_n]=='t')
        {
             if(random<=0.2)break;
	}	            
        else
	{ 
             if(random>0.2)break;
	}
    }        
    random=RngStream_RandU01(RS);
    switch (Genome[pos_n])
    {
        case 'a':
            if(random<=0.2)
                Genome[pos_n]='c';
            else
            {
                if((random-0.2)<=0.3)
                    Genome[pos_n]='t';
                else
                    Genome[pos_n]='g';
            }
            break;
        case 't':
            if(random<=0.2)
                Genome[pos_n]='g';
            else
            {
                if((random-0.2)<=0.3)
                    Genome[pos_n]='a';
                else
                    Genome[pos_n]='c';
            }
            break;           
        case 'c':
            if(random<=0.2)
                Genome[pos_n]='g';
            else
            {
                if((random-0.2)<=0.3)
                    Genome[pos_n]='t';
                else
                    Genome[pos_n]='a';
            }
            break;           
        case 'g':  
            if(random<=0.2)
                Genome[pos_n]='c';
            else
            {
                if((random-0.2)<=0.3)
                    Genome[pos_n]='a';
                else
                    Genome[pos_n]='t';
            }
            break;
    }   
    #else
    pos_n=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
    nucleotide=set_base_pair(RngStream_RandU01(RS));
    while(nucleotide==Genome[pos_n])
        nucleotide=set_base_pair(RngStream_RandU01(RS));
    Genome[pos_n]=nucleotide;
    #endif
    
    pos_g=pos_n/CISREG_LEN;    
    genotype->recalc_TFBS[pos_g]=1;  /*recalc TFBS*/
    genotype->clone_info[pos_g]=1;   /*restore info in clone_cell*/
    update_cisreg_cluster(genotype,pos_g,'s');    
    
    /*record mutation info*/
    mut_record->pos_n=pos_n;
    mut_record->pos_g=pos_g;
    mut_record->nuc_diff[0]=Genome[pos_n];
}

void reproduce_substitution(Genotype *genotype, Mutation *mut_record)
{    
    char *Genome;
    Genome= &genotype->cisreg_seq[0][0];
    Genome[mut_record->pos_n]=mut_record->nuc_diff[0];   
        
    genotype->recalc_TFBS[mut_record->pos_g]=1;  /*recalc TFBS*/
    genotype->clone_info[mut_record->pos_g]=1;   /*restore info in clone_cell*/
    update_cisreg_cluster(genotype,mut_record->pos_g,'s');    
}

/* MAXCOPIES has to be 1*/
void mut_insertion(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int inset_size=0;
    int pos_g,pos_n,i;
    char n; 

    inset_size=RngStream_RandInt(RS,1,max_inset);
    pos_g=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    pos_n=RngStream_RandInt(RS,0,CISREG_LEN-1);                /* at which new seq will be inserted*/
    
    if (pos_n+inset_size>CISREG_LEN) inset_size=CISREG_LEN-pos_n;

    for(i=1;i<=CISREG_LEN-inset_size-pos_n;i++)    
        genotype->cisreg_seq[pos_g][CISREG_LEN-i]=genotype->cisreg_seq[pos_g][CISREG_LEN-inset_size-i];
    for(i=pos_n;i<pos_n+inset_size;i++)
    {        					
        n=set_base_pair(RngStream_RandU01(RS));	
        genotype->cisreg_seq[pos_g][i]=n;
        /*record mutation info*/
        mut_record->nuc_diff[i-pos_n]=n;
    }
    
    genotype->recalc_TFBS[pos_g]=1;  /*recalc TFBS*/
    genotype->clone_info[pos_g]=1;   /*restore info in clone_cell*/
    update_cisreg_cluster(genotype,pos_g,'i');    
    
    /*record mutation info*/
    mut_record->pos_n=pos_n;
    mut_record->pos_g=pos_g;    
}

void reproduce_insertion(Genotype *genotype, Mutation *mut_record)
{
    int inset_size;    
    int pos_g,pos_n,i;    
    
    inset_size=(int)strlen(mut_record->nuc_diff);
    pos_g=mut_record->pos_g;
    pos_n=mut_record->pos_n;               /* at which new seq will be inserted*/
    
    if (pos_n+inset_size>CISREG_LEN) inset_size=CISREG_LEN-pos_n;

    for(i=1;i<=CISREG_LEN-inset_size-pos_n;i++)    
        genotype->cisreg_seq[pos_g][CISREG_LEN-i]=genotype->cisreg_seq[pos_g][CISREG_LEN-inset_size-i];
    for(i=pos_n;i<pos_n+inset_size;i++)      	
        genotype->cisreg_seq[pos_g][i]=mut_record->nuc_diff[i-pos_n];    
    genotype->recalc_TFBS[pos_g]=1;  /*recalc TFBS*/
    genotype->clone_info[pos_g]=1;   /*restore info in clone_cell*/
    update_cisreg_cluster(genotype,pos_g,'i');   
}

void mut_partial_deletion(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int delet_size=0;
    int pos_g, pos_n, i;
    char n;    

    delet_size = RngStream_RandInt(RS,1,max_delet);              
    pos_g=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);/* from which a seq will be deleted */
    pos_n=RngStream_RandInt(RS,0,CISREG_LEN-1);
    
    /*record mutation info*/
    mut_record->pos_g=pos_g;
    mut_record->pos_n=pos_n;

    if (pos_n+delet_size>CISREG_LEN)	/* if only the tail is deleted*/
    {
        delet_size=CISREG_LEN-pos_n;
        
        for(i=pos_n;i<pos_n+delet_size;i++)
        {					
            n=set_base_pair(RngStream_RandU01(RS));
            genotype->cisreg_seq[pos_g][i]=n;
             /*record mutation info*/
            mut_record->nuc_diff[i-pos_n]=n;
        }				
    }
    else /* otherwise, join the two fragments aside the deletion */
    {
        for(i=pos_n;i<CISREG_LEN-delet_size;i++)
        {
            genotype->cisreg_seq[pos_g][i]=genotype->cisreg_seq[pos_g][i+delet_size];
        }                				
        for(i=CISREG_LEN-delet_size;i<CISREG_LEN;i++) /* and fill the gab by generating new seq */
        {				
            n= set_base_pair(RngStream_RandU01(RS));
            genotype->cisreg_seq[pos_g][i]=n;
             /*record mutation info*/
            mut_record->nuc_diff[i-(CISREG_LEN-delet_size)]=n;
        }						
    }	
    
    genotype->recalc_TFBS[pos_g]=1;  /*recalc TFBS*/
    genotype->clone_info[pos_g]=1;   /*restore info in clone_cell*/
    update_cisreg_cluster(genotype,pos_g,'p');  
}

void reproduce_partial_deletion(Genotype *genotype, Mutation *mut_record)
{
    int delet_size;   
    int pos_g, pos_n, i;    
    
    delet_size=(int)strlen(mut_record->nuc_diff);		
    pos_g=mut_record->pos_g;    		
    pos_n=mut_record->pos_n;                /* from which a seq will be deleted */
   
    if (pos_n+delet_size>CISREG_LEN)	/* if only the tail is deleted*/
    {
        delet_size=CISREG_LEN-pos_n;
        
        for(i=pos_n;i<pos_n+delet_size;i++)
        {
            genotype->cisreg_seq[pos_g][i]=mut_record->nuc_diff[i-pos_n];            
        }				
    }
    else /* else, join the two fragments aside the deletion */
    {
        for(i=pos_n;i<CISREG_LEN-delet_size;i++)
        {
            genotype->cisreg_seq[pos_g][i]=genotype->cisreg_seq[pos_g][i+delet_size];
        }                				
        for(i++;i<CISREG_LEN;i++) /* and fill the gab by generating new seq */
        {          
            genotype->cisreg_seq[pos_g][i]=mut_record->nuc_diff[i-pos_n];           
        }						
    }	
    
    genotype->recalc_TFBS[pos_g]=1;  /*recalc TFBS*/
    genotype->clone_info[pos_g]=1;   /*restore info in clone_cell*/
    update_cisreg_cluster(genotype,pos_g,'p'); 
}

void mut_whole_gene_deletion(Genotype *genotype, Mutation *mut_record, RngStream RS) // any gene can be deleted
{
    int pos_g, offset, i,j,cluster_id,cluster_id_copy;
    char *temp1;
    int protein_id, n_gene_of_selection_protein, puppet_gene,pos_puppet_gene;    
   
    /* check first whether the fitness proteins have extra copies of genes*/   
    n_gene_of_selection_protein=genotype->protein_pool[genotype->nproteins-1][0][0];
    
    if(n_gene_of_selection_protein==1)
    {
        protein_id=genotype->nproteins-1;
        while(protein_id==genotype->nproteins-1)
        {   
            pos_g=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
            protein_id=genotype->which_protein[pos_g];
        }
    }
    else 
        pos_g=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    
    /*record mutation info*/
    mut_record->pos_g=pos_g;
    
    /*if the original selection genes are to be deleted, we find one of their duplicated genes and copy it to 
     the default selection genes' position. Then we deleted the duplicated gene. For example, assuming ngenes=10,
     gene 8 is selection gene A and is marked for deletion. We know gene 5 was duplicated from gene 8, so we 
     overwrite gene 8 with all the info of gene 5, then we delete gene 5 (the puppet). We go through these troubles to keep 
     the structure of genotype information.*/
    if(pos_g==genotype->ngenes-1)
    {
        cluster_id=genotype->which_cluster[pos_g];
        if(genotype->cisreg_cluster[cluster_id][1]!=-1)
        {            
            puppet_gene=genotype->cisreg_cluster[cluster_id][0];
        }
        else /*this means pos_g has a unique cisreg, which is a trouble*/
        {
            protein_id=genotype->which_protein[pos_g];
            pos_puppet_gene=0;
            if(genotype->protein_pool[protein_id][1][pos_puppet_gene]==pos_g)pos_puppet_gene++;
            puppet_gene=genotype->protein_pool[protein_id][1][pos_puppet_gene];
            /*dirty work before update_cisreg_cluster*/
            /*first add pos_g to the cis_reg cluster of i*/
            cluster_id_copy=genotype->which_cluster[puppet_gene];            
            j=0;
            while(genotype->cisreg_cluster[cluster_id_copy][j]!=-1)j++;
            genotype->cisreg_cluster[cluster_id_copy][j]=pos_g;            
            genotype->which_cluster[pos_g]=cluster_id_copy;
            /*then delete pos_g's cisreg cluster*/
            cluster_id_copy=cluster_id;            
            while(genotype->cisreg_cluster[cluster_id_copy][0]!=-1)
            {
                /*reset cluster=cluster_id*/
                j=0;
                while(genotype->cisreg_cluster[cluster_id_copy][j]!=-1)
                {
                    genotype->cisreg_cluster[cluster_id_copy][j]=-1;
                    j++;
                }
                /*then copy from cluster_id+1, and update which_cluster*/                    
                j=0;
                while(genotype->cisreg_cluster[cluster_id_copy+1][j]!=-1)
                {                        
                    genotype->cisreg_cluster[cluster_id_copy][j]=genotype->cisreg_cluster[cluster_id_copy+1][j];
                    genotype->which_cluster[genotype->cisreg_cluster[cluster_id_copy+1][j]]--;
                    j++;
                }
                cluster_id_copy++;
            }
        }       
        /*overwrite pos_g with the info of puppet_gene*/
        for(j=0;j<CISREG_LEN;j++)        
            genotype->cisreg_seq[pos_g][j]=genotype->cisreg_seq[puppet_gene][j];
        genotype->min_act_to_transc[pos_g]=genotype->min_act_to_transc[puppet_gene];  /* shift elements in the array*/
        genotype->pic_disassembly[pos_g]=genotype->pic_disassembly[puppet_gene];           
        genotype->mRNAdecay[pos_g]=genotype->mRNAdecay[puppet_gene];
        genotype->proteindecay[pos_g]=genotype->proteindecay[puppet_gene];
        genotype->translation[pos_g]=genotype->translation[puppet_gene]; 
        genotype->recalc_TFBS[pos_g]=1;
        genotype->clone_info[pos_g]=1;       
        /*now pos_g is essentially i, so delete the original gene i */
        pos_g=puppet_gene;        
    }
    
    temp1 = &genotype->cisreg_seq[pos_g][0];			
    offset=CISREG_LEN;
    /* remove pos_g from cis_seq */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-pos_g-1);i++) 
    {				
        *temp1=*(temp1+offset);         /* move sequence to fill the gap */
        temp1++;				
    }  
    
    /* if pos_g is a tf gene*/
    protein_id=genotype->which_protein[pos_g];   
    if(protein_id<genotype->nproteins-1)
    {  
        /*if this tf has only one copy of gene, then we'll delete the binding seq of this tf*/ 
        if(genotype->protein_pool[protein_id][0][0]==1) 
        {            
            /* remove from binding seq */
            temp1=&genotype->tf_seq[protein_id][0];
            offset=TF_ELEMENT_LEN;
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp1=*(temp1+offset);
                temp1++;
            }
            /* remove from rc binding seq */    
            temp1=&genotype->tf_seq_rc[protein_id][0];
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp1=*(temp1+offset);
                temp1++;
            }   
        }
    } 

    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
    for(i=pos_g;i<genotype->ngenes;i++)
    {
        genotype->min_act_to_transc[i]=genotype->min_act_to_transc[i+1];
        genotype->pic_disassembly[i]=genotype->pic_disassembly[i+1];             /* shift elements in the array*/
        genotype->mRNAdecay[i]=genotype->mRNAdecay[i+1];
        genotype->proteindecay[i]=genotype->proteindecay[i+1];
        genotype->translation[i]=genotype->translation[i+1];         
        genotype->recalc_TFBS[i]=1;
        genotype->clone_info[i]=1; 
    }

    /* if necessary, ntfgenes as well*/
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes--;    
    /* now change protein_pool and cisreg_cluster*/   
    update_protein_pool(genotype,protein_id,pos_g,'w'); 
    update_cisreg_cluster(genotype,pos_g,'w');  
    genotype->ngenes--;   
}

void reproduce_whole_gene_deletion(Genotype *genotype, Mutation *mut_record) // any gene can be deleted
{    
    int pos_g, offset, i,j, cluster_id, cluster_id_copy,puppet_gene,pos_puppet_gene;
    char *temp1;    
    int protein_id;
       
    pos_g=mut_record->pos_g;
    
    /*if the original selection genes are to be deleted, we find one of their duplicated genes and copy it to 
     the default selection genes' position. Then we deleted the duplicated gene. For example, assuming ngenes=10,
     gene 8 is selection gene A and is marked for deletion. We know gene 5 was duplicated from gene 8, so we 
     overwrite gene 8 with all the info of gene 5, then we delete gene 5 (the puppet). We go through these troubles to keep 
     the structure of genotype information.*/
    if(pos_g==genotype->ngenes-1)
    {
        cluster_id=genotype->which_cluster[pos_g];
        if(genotype->cisreg_cluster[cluster_id][1]!=-1)
        {            
            puppet_gene=genotype->cisreg_cluster[cluster_id][0];
        }
        else /*this means pos_g has a unique cisreg, which is a trouble*/
        {
            protein_id=genotype->which_protein[pos_g];
            pos_puppet_gene=0;
            if(genotype->protein_pool[protein_id][1][pos_puppet_gene]==pos_g)pos_puppet_gene++;
            puppet_gene=genotype->protein_pool[protein_id][1][pos_puppet_gene];
            /*dirty work before update_cisreg_cluster*/
            /*first add pos_g to the cis_reg cluster of i*/
            cluster_id_copy=genotype->which_cluster[puppet_gene];            
            j=0;
            while(genotype->cisreg_cluster[cluster_id_copy][j]!=-1)j++;
            genotype->cisreg_cluster[cluster_id_copy][j]=pos_g;            
            genotype->which_cluster[pos_g]=cluster_id_copy;
            /*then delete pos_g's cisreg cluster*/
            cluster_id_copy=cluster_id;            
            while(genotype->cisreg_cluster[cluster_id_copy][0]!=-1)
            {
                /*reset cluster=cluster_id*/
                j=0;
                while(genotype->cisreg_cluster[cluster_id_copy][j]!=-1)
                {
                    genotype->cisreg_cluster[cluster_id_copy][j]=-1;
                    j++;
                }
                /*then copy from cluster_id+1, and update which_cluster*/                    
                j=0;
                while(genotype->cisreg_cluster[cluster_id_copy+1][j]!=-1)
                {                        
                    genotype->cisreg_cluster[cluster_id_copy][j]=genotype->cisreg_cluster[cluster_id_copy+1][j];
                    genotype->which_cluster[genotype->cisreg_cluster[cluster_id_copy+1][j]]--;
                    j++;
                }
                cluster_id_copy++;
            }
        }       
        /*overwrite pos_g with the info of puppet_gene*/
        for(j=0;j<CISREG_LEN;j++)        
            genotype->cisreg_seq[pos_g][j]=genotype->cisreg_seq[puppet_gene][j];
        genotype->min_act_to_transc[pos_g]=genotype->min_act_to_transc[puppet_gene];  /* shift elements in the array*/
        genotype->pic_disassembly[pos_g]=genotype->pic_disassembly[puppet_gene];           
        genotype->mRNAdecay[pos_g]=genotype->mRNAdecay[puppet_gene];
        genotype->proteindecay[pos_g]=genotype->proteindecay[puppet_gene];
        genotype->translation[pos_g]=genotype->translation[puppet_gene]; 
        genotype->recalc_TFBS[pos_g]=1;
        genotype->clone_info[pos_g]=1;       
        /*now pos_g is essentially i, so delete the original gene i */
        pos_g=puppet_gene;        
    }
    
    temp1 = &genotype->cisreg_seq[pos_g][0];			
    offset=CISREG_LEN;
    /* remove pos_g from cis_seq */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-pos_g-1);i++) 
    {				
        *temp1=*(temp1+offset);         /* move sequence to fill the gap */
        temp1++;				
    }   
    
    /* if pos_g is a tf gene*/
    protein_id=genotype->which_protein[pos_g];   
    if(protein_id<genotype->nproteins-1)
    {  
        /*if this tf has only one copy of gene, then we'll delete the binding seq of this tf*/ 
        if(genotype->protein_pool[protein_id][0][0]==1) 
        {            
            /* remove from binding seq */
            temp1=&genotype->tf_seq[protein_id][0];
            offset=TF_ELEMENT_LEN;
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp1=*(temp1+offset);
                temp1++;
            }
            /* remove from rc binding seq */    
            temp1=&genotype->tf_seq_rc[protein_id][0];
            for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
            {
                *temp1=*(temp1+offset);
                temp1++;
            }   
        }
    } 
    
    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
    for(i=pos_g;i<genotype->ngenes;i++)
    {
        genotype->min_act_to_transc[i]=genotype->min_act_to_transc[i+1];
        genotype->pic_disassembly[i]=genotype->pic_disassembly[i+1];             /* shift elements in the array*/
        genotype->mRNAdecay[i]=genotype->mRNAdecay[i+1];
        genotype->proteindecay[i]=genotype->proteindecay[i+1];
        genotype->translation[i]=genotype->translation[i+1];         
        genotype->recalc_TFBS[i]=1;
        genotype->clone_info[i]=1; 
    }   
   
    /* if necessary, ntfgenes as well*/
    if(protein_id<genotype->nproteins-1)
        genotype->ntfgenes--;
    
    /* now change protein_pool and cisreg_cluster*/
    update_protein_pool(genotype,protein_id,pos_g,'w'); 
    update_cisreg_cluster(genotype,pos_g,'w');  
    genotype->ngenes--;   
}

void mut_duplication(Genotype *genotype, Mutation *mut_record, RngStream RS) //any gene can be duplicated
{
    int pos_g, pos_g_copy, i, protein_id;
    char *temp1, *temp2;  
    
#if RdcPdup     /* the method of reducing probability*/ 
    float random1;
    int random2,tag; 
    tag=0;
    while(tag==0)
    {
        random1=RngStream_RandU01(RS);
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(random1<=genotype->P_dup_per_protein[i])
            {    
                random2=RngStream_RandInt(RS,0,genotype->protein_pool[i][0][0]-1);
                pos_g=genotype->protein_pool[i][1][random2];
                tag=1;
                break;
            }
            else
                random1-=genotype->P_dup_per_protein[i];
        }
    }
#else       /* the method of reducing growth rate */
    pos_g=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);     
#endif     
    
    if(pos_g==genotype->ngenes-1)
        pos_g_copy=pos_g+1; /* note that if the selection genes are to be duplicated, shifting sequences and info will cause problem*/
    else
        pos_g_copy=pos_g;
    
    /*record mutation info*/
    mut_record->pos_g=pos_g;    
    /* copy the promoter*/
    temp1=&genotype->cisreg_seq[pos_g_copy][0]; /* points to the gene to be duplicated*/
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][CISREG_LEN-1]; /* points to the end of the selection gene */
    
    for(i=0;i<CISREG_LEN;i++) 
    {
        *(temp2+CISREG_LEN)=*temp2; /* shift the sequences of the selection genes CISREG_LEN bp */
        temp2--;
    }    
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0]; /* point temp2 to the start of the original position of the selection gene */    
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;  /* put the duplicated gene at the original place of the selection gene */ 
    
    /* add it to PIC_assembly, mRNAdecay, proteinDecay, translation, and which_protein*/
    /* shift the selection genes to make a slot*/
    genotype->min_act_to_transc[genotype->ngenes]=genotype->min_act_to_transc[genotype->ngenes-1];
    genotype->pic_disassembly[genotype->ngenes]=genotype->pic_disassembly[genotype->ngenes-1]; 
    genotype->mRNAdecay[genotype->ngenes]=genotype->mRNAdecay[genotype->ngenes-1];
    genotype->proteindecay[genotype->ngenes]=genotype->proteindecay[genotype->ngenes-1];
    genotype->translation[genotype->ngenes]=genotype->translation[genotype->ngenes-1];   
    genotype->clone_info[genotype->ngenes]=1;
    genotype->recalc_TFBS[genotype->ngenes]=1; 
        
    /* copy and paste info to the slot*/
    genotype->min_act_to_transc[genotype->ngenes-1]=genotype->min_act_to_transc[pos_g_copy];
    genotype->pic_disassembly[genotype->ngenes-1]=genotype->pic_disassembly[pos_g_copy];
    genotype->mRNAdecay[genotype->ngenes-1]=genotype->mRNAdecay[pos_g_copy];
    genotype->proteindecay[genotype->ngenes-1]=genotype->proteindecay[pos_g_copy];
    genotype->translation[genotype->ngenes-1]=genotype->translation[pos_g_copy];   
    genotype->clone_info[genotype->ngenes-1]=1;
    genotype->recalc_TFBS[genotype->ngenes-1]=1;   
    
    /* update protein_pool*/
    protein_id=genotype->which_protein[pos_g];    
    update_protein_pool(genotype,protein_id,pos_g,'d'); 
    /* update cisreg_cluster*/    
    update_cisreg_cluster(genotype,pos_g,'d');
    /* update gene numbers*/  
    if(protein_id<genotype->nproteins-1)/* note duplication do not change nproteins*/            
        genotype->ntfgenes++;     
    genotype->ngenes++;   
}

void reproduce_gene_duplication(Genotype *genotype, Mutation *mut_record) //any gene can be duplicated
{   
    int pos_g,pos_g_copy, i, protein_id;
    char *temp1, *temp2;
    
    pos_g=mut_record->pos_g;
    
    if(pos_g==genotype->ngenes-1)
        pos_g_copy=pos_g+1; /* note that if the selection genes are to be duplicated, shifting sequences and info will cause problem*/
    else
        pos_g_copy=pos_g;
    
    /*record mutation info*/
    mut_record->pos_g=pos_g;    
    /* copy the promoter*/
    temp1=&genotype->cisreg_seq[pos_g_copy][0]; /* points to the gene to be duplicated*/
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][CISREG_LEN-1]; /* points to the end of the selection gene */
    
    for(i=0;i<CISREG_LEN;i++) 
    {
        *(temp2+CISREG_LEN)=*temp2; /* shift the sequences of the selection genes CISREG_LEN bp */
        temp2--;
    }    
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0]; /* point temp2 to the start of the original position of the selection gene */    
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;  /* put the duplicated gene at the original place of the selection gene */ 
    
    /* add it to PIC_assembly, mRNAdecay, proteinDecay, translation, and which_protein*/
    /* shift the selection genes to make a slot*/
    genotype->min_act_to_transc[genotype->ngenes]=genotype->min_act_to_transc[genotype->ngenes-1];
    genotype->pic_disassembly[genotype->ngenes]=genotype->pic_disassembly[genotype->ngenes-1]; 
    genotype->mRNAdecay[genotype->ngenes]=genotype->mRNAdecay[genotype->ngenes-1];
    genotype->proteindecay[genotype->ngenes]=genotype->proteindecay[genotype->ngenes-1];
    genotype->translation[genotype->ngenes]=genotype->translation[genotype->ngenes-1];   
    genotype->clone_info[genotype->ngenes]=1;
    genotype->recalc_TFBS[genotype->ngenes]=1;
    
    /* copy and paste info to the slot*/
    genotype->min_act_to_transc[genotype->ngenes-1]=genotype->min_act_to_transc[pos_g_copy];
    genotype->pic_disassembly[genotype->ngenes-1]=genotype->pic_disassembly[pos_g_copy];
    genotype->mRNAdecay[genotype->ngenes-1]=genotype->mRNAdecay[pos_g_copy];
    genotype->proteindecay[genotype->ngenes-1]=genotype->proteindecay[pos_g_copy];
    genotype->translation[genotype->ngenes-1]=genotype->translation[pos_g_copy];   
    genotype->clone_info[genotype->ngenes-1]=1;
    genotype->recalc_TFBS[genotype->ngenes-1]=1;
    
    /* update protein_pool*/
    protein_id=genotype->which_protein[pos_g];
    update_protein_pool(genotype,protein_id,pos_g,'d'); 
    /* update cisreg_cluster*/    
    update_cisreg_cluster(genotype,pos_g,'d');
    /* update gene numbers*/  
    if(protein_id<genotype->nproteins-1)/* note duplication do not change nproteins*/ 
        genotype->ntfgenes++;      
    genotype->ngenes++;     
}

void mut_binding_sequence(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int pos_g, pos_n, protein_id, i;
    char n;      
    char *tf_seq, *tf_seq_rc,*temp1,*temp2; 
    
    pos_g=RngStream_RandInt(RS,0,genotype->ngenes-2); // ngenes-1 is the selection gene, so tf genes must be 0 to ngenes-2
    while(genotype->which_protein[pos_g]==genotype->nproteins-1)
        pos_g=RngStream_RandInt(RS,0,genotype->ngenes-2);
    
    protein_id=genotype->which_protein[pos_g];
    /*if this TF has more than one copies of gene, then the mutation creates a new tf
     *which requires a new slot in tf_seq and tf_seq_rc to store the new binding seq*/
    if(genotype->protein_pool[protein_id][0][0]>1)
    {
        /*first, copy the binding seq to the new slot*/
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }
        /*point tf_seq and tf_seq_rc to the new slot*/
        tf_seq=&genotype->tf_seq[genotype->nproteins-1][0];
        tf_seq_rc=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        /*and update protein pool*/
        update_protein_pool(genotype,protein_id,pos_g,'c');  
    }    
    else /*if this tf has only one copy of gene, no new slot is required*/
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
    }
    
    pos_n=RngStream_RandInt(RS,0,TF_ELEMENT_LEN-1);
    n=set_base_pair(RngStream_RandU01(RS));
        
    while (n == tf_seq[pos_n])
        n=set_base_pair(RngStream_RandU01(RS));
    tf_seq[pos_n]=n;    
    /*record mutation info*/
    mut_record->pos_g=pos_g;
    mut_record->pos_n=pos_n;
    mut_record->nuc_diff[0]=n;    
    /* update the complement sequence*/
    switch (n)
    {
        case 'g':
            tf_seq_rc[pos_n]='c'; break;
        case 'c':
            tf_seq_rc[pos_n]='g'; break;
        case 'a':
            tf_seq_rc[pos_n]='t'; break;
        case 't':
            tf_seq_rc[pos_n]='a'; break;
    } 
        
    for(i=0;i<genotype->ngenes;i++) 
    {
        genotype->recalc_TFBS[i]=1;   /* recalculate binding sites on every promoter */
        genotype->clone_info[i]=1;   /* copy info back to every gene in clone_cell */
    }
}

void reproduce_mut_binding_sequence(Genotype *genotype, Mutation *mut_record)
{    
    int pos_g, pos_n, protein_id, i;
    
    char *tf_seq, *tf_seq_rc, *temp1,*temp2; 
    
    pos_g=mut_record->pos_g;
    protein_id=genotype->which_protein[pos_g];
    /*if this TF has more than one copies of gene, then the mutation creates a new tf
     *which requires a new slot in tf_seq and tf_seq_rc to store the new binding seq*/
    if(genotype->protein_pool[protein_id][0][0]>1)
    {
        /*first, copy the binding seq to the new slot*/
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
        temp1=&genotype->tf_seq[genotype->nproteins-1][0];
        temp2=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_seq++;
            *temp2++=*tf_seq_rc++;
        }
        /*point tf_seq and tf_seq_rc to the new slot*/
        tf_seq=&genotype->tf_seq[genotype->nproteins-1][0];
        tf_seq_rc=&genotype->tf_seq_rc[genotype->nproteins-1][0];
        /*and update protein pool*/
        update_protein_pool(genotype,protein_id,pos_g,'c');  
    }    
    else /*if this tf has only one copy of gene, no new slot is required*/
    {    
        tf_seq=&genotype->tf_seq[protein_id][0];
        tf_seq_rc=&genotype->tf_seq_rc[protein_id][0];
    }
        			
    pos_n=mut_record->pos_n;
        
    tf_seq[pos_n]=mut_record->nuc_diff[0];
    
    /* update the complement sequence*/
    switch (tf_seq[pos_n])
    {
        case 'g':
            tf_seq_rc[pos_n]='c'; break;
        case 'c':
            tf_seq_rc[pos_n]='g'; break;
        case 'a':
            tf_seq_rc[pos_n]='t'; break;
        case 't':
            tf_seq_rc[pos_n]='a'; break;
    } 
    
    for(i=0;i<genotype->ngenes;i++) 
    {
        genotype->recalc_TFBS[i]=1;   /* recalculate binding sites on every promoter */
        genotype->clone_info[i]=1;   /* copy info back to every gene in clone_cell */
    }
}

/* For the moment, only mRNA_decay, translation, protein_decay, and pic_disassembly 
 * will be mutated. We assume a mutation attacks the four constants with equal 
 * probability. 
 */
void mut_kinetic_constant(Genotype *genotype, Mutation *mut_record, float kdis[NUM_K_DISASSEMBLY],RngStream RS)
{
    float random1, random2;
    int pos_kdis, pos_g;
      
    random1=RngStream_RandU01(RS)-proportion_mut_identity-proportion_mut_koff-proportion_mut_binding_seq;
    pos_g=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    
    /*record mutation info*/
    mut_record->pos_g=pos_g;
    
    if(random1<=proportion_mut_kdis) /* mut kdis */
    {   
        pos_kdis=RngStream_RandInt(RS,0,NUM_K_DISASSEMBLY-1);      
    #if UNLIMITED_MUTATION
        genotype->pic_disassembly[pos_g]=kdis[RngStream_RandInt(RS,0,NUM_K_DISASSEMBLY-1)];
    #else
        genotype->pic_disassembly[pos_g]=genotype->pic_disassembly[pos_g]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);  
        while(genotype->pic_disassembly[pos_g]<MIN_PICDISASSEMBLE || genotype->pic_disassembly[pos_g]>MAX_PICDISASSEMBLE)
            genotype->pic_disassembly[pos_g]=genotype->pic_disassembly[pos_g]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);
    #endif        
        /*record mutation info*/
        mut_record->kinetic_type=0;
        mut_record->kinetic_diff=kdis[pos_kdis];
    }
    else 
    {
        random1-=proportion_mut_kdis;
        if(random1<=proportion_mut_mRNA_decay) /* mut mRNAdecay */
        {
        #if UNLIMITED_MUTATION
            random2 = exp(0.4909*gasdev(RS)-3.20304);
            while(random2<MIN_mRNA_decay)
                random2 = exp(0.4909*gasdev(RS)-3.20304);
            genotype->mRNAdecay[pos_g]=random2;
        #else
            genotype->mRNAdecay[pos_g]=genotype->mRNAdecay[pos_g]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut); 
            while(genotype->mRNAdecay[pos_g]>MAX_mRNA_decay || genotype->mRNAdecay[pos_g]<MIN_mRNA_decay)
                genotype->mRNAdecay[pos_g]=genotype->mRNAdecay[pos_g]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);
        #endif
            /*record mutation info*/
            mut_record->kinetic_type=1;
            mut_record->kinetic_diff=genotype->mRNAdecay[pos_g];
        }
        else 
        {
            random1-=proportion_mut_mRNA_decay;
            if(random1<=proportion_mut_translation_rate) /* mut translation */
            {
        #if UNLIMITED_MUTATION
                random2= exp(0.7406*gasdev(RS)+4.56);
                while(random2>MAX_TRANSLATION)
                    random2=exp(0.7406*gasdev(RS)+4.56);
                genotype->translation[pos_g]=random2;
        #else
                genotype->translation[pos_g]=genotype->translation[pos_g]*exp(sd_mut_effect*gasdev(RS)-bias_in_mut);
                while(genotype->translation[pos_g]>MAX_TRANSLATION || genotype->translation[pos_g]< MIN_TRANSLATION)
                    genotype->translation[pos_g]=genotype->translation[pos_g]*exp(sd_mut_effect*gasdev(RS)-bias_in_mut);
        #endif
                mut_record->kinetic_type=2;
                mut_record->kinetic_diff=genotype->translation[pos_g];
            }
            else /* mut protein decay */
            {  
            #if UNLIMITED_MUTATION
                random2=genotype->proteindecay[pos_g];         
                while (random2 < 0.0 ) 
                {
                    if (RngStream_RandU01(RS) < 0.08421)
                        random2 = (float)EPSILON;
                    else random2 = exp(0.7874*gasdev(RS)-3.7665);
                }
                genotype->proteindecay[pos_g]=random2;
            #else
                genotype->proteindecay[pos_g]=genotype->proteindecay[pos_g]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut);  
                while(genotype->proteindecay[pos_g]>MAX_protein_decay || genotype->proteindecay[pos_g]<MIN_protein_decay)
                    genotype->proteindecay[pos_g]=genotype->proteindecay[pos_g]*exp(sd_mut_effect*gasdev(RS)+bias_in_mut); 
            #endif
                mut_record->kinetic_type=3;
                mut_record->kinetic_diff=genotype->proteindecay[pos_g];
            } 
        }
    }
}

void reproduce_mut_kinetic_constant(Genotype *genotype, Mutation *mut_record)
{    
    int pos_g;    
        
    pos_g=mut_record->pos_g; /* which gene */   
    
    switch (mut_record->kinetic_type)
    {
        case 0: /* mut kdis */
            genotype->pic_disassembly[pos_g]=mut_record->kinetic_diff;
            break;
        case 1: /* mut mRNAdecay */ 
            genotype->mRNAdecay[pos_g]=mut_record->kinetic_diff;
            break;
        case 2: /* mut translation */                      
            genotype->translation[pos_g]=mut_record->kinetic_diff;           
            break;        
        case 3: /* mut protein decay */
            genotype->proteindecay[pos_g]=mut_record->kinetic_diff;            
            break;        
    }
}

void mut_identity(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int tf_id,protein_id,i;   
    char *tf_seq,*tf_seq_rc,*temp1,*temp2;
 
    tf_id = RngStream_RandInt(RS,2,genotype->ngenes-2); // the first signal tf must be an activator, therefore is not subject to mutation

    protein_id=genotype->which_protein[tf_id];    
    while(protein_id==genotype->nproteins-1)
    {
        tf_id = RngStream_RandInt(RS,2,genotype->ngenes-2); // the signal tf must be an activator, therefore is not subject to mutation
        protein_id=genotype->which_protein[tf_id];
    }    
    
    mut_record->pos_g=tf_id;
       
    /* if this tf gene has more than one copies, the mutation increases nproteins*/
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
        genotype->activating[genotype->nproteins-2]=(genotype->activating[genotype->nproteins-2]==1)?0:1;  /* update_protein_pool put the new protein at nproteins-1 and then increases nproteins by 1, 
                                                                                                            * so the new protein is at nproteins-2 now.
                                                                                                            * Note that N_act and N_rep is updated in update_protein_pool*/        
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
    {
        genotype->recalc_TFBS[i]=1;   /* recalculate binding sites on every promoter */
        genotype->clone_info[i]=1;   /* copy info back to every gene in clone_cell */
    }    
}

void reproduce_mut_identity(Genotype *genotype, Mutation *mut_record)
{
    int tf_id, protein_id,i;
    char *tf_seq,*tf_seq_rc,*temp1,*temp2;
    
    tf_id = mut_record->pos_g;
    protein_id=genotype->which_protein[tf_id]; 
       
     /* if this tf gene has more than one copies, the mutation inreases nproteins*/
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
        genotype->activating[genotype->nproteins-2]=(genotype->activating[genotype->nproteins-2]==1)?0:1;  /* update_protein_pool put the new protein at nproteins-1 and then increases nproteins by 1, 
                                                                                                            * so the new protein is at nproteins-2 now.
                                                                                                            * Note that N_act and N_rep is updated in update_protein_pool*/        
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
    {
        genotype->recalc_TFBS[i]=1;   /* recalculate binding sites on every promoter */
        genotype->clone_info[i]=1;   /* copy info back to every gene in clone_cell */
    }    
}

void mut_koff(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int i;
    float new_koff;
    int tf_id,protein_id;
    char *tf_seq,*tf_seq_rc,*temp1,*temp2;

    tf_id=RngStream_RandInt(RS,0,genotype->ngenes-2);
    protein_id=genotype->which_protein[tf_id];
    while(protein_id==genotype->nproteins-1)
    {
        tf_id=RngStream_RandInt(RS,0,genotype->ngenes-2);
        protein_id=genotype->which_protein[tf_id];
    }       

    /*generate a new koff */ 
#if UNLIMITED_MUTATION
    new_koff=RngStream_RandU01(RS);
    while(new_koff<MIN_koff || new_koff>MAX_koff)
 	new_koff=RngStream_RandU01(RS);
#else
    float old_koff;
    old_koff=genotype->koff[protein_id];
    new_koff=old_koff*exp(gasdev(RS)*sd_mut_effect+bias_in_mut);
    while(new_koff<MIN_koff || new_koff>MAX_koff)
        new_koff=old_koff*exp(gasdev(RS)*sd_mut_effect+bias_in_mut);
#endif
    
    /*update protein_pool and koff*/    
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
        genotype->koff[genotype->nproteins-2]=new_koff;  /* update_protein_pool put the new protein at nproteins-1 and then increases nproteins by 1, 
                                                                                                            * so the new protein is at nproteins-2 now.*/
    }                                                                                                           
    else
        genotype->koff[protein_id]=new_koff;       

    
#if SET_BS_MANUALLY 
    int j;
    for(i=1;i<genotype->ngenes-1;i++)
    {
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {
            if(genotype->all_binding_sites[i][j].tf_id==tf_id)
                genotype->all_binding_sites[i][j].Koff=genotype->koff[protein_id]*
                        pow(KR,(float)genotype->all_binding_sites[i][j].mis_match/3.0);
        }
    }
#endif    

    mut_record->pos_g=tf_id;
    mut_record->kinetic_diff=genotype->koff[tf_id]; 
    
    for(i=0;i<genotype->ngenes;i++) 
    {
        genotype->recalc_TFBS[i]=1;   /* recalculate binding sites on every promoter */
        genotype->clone_info[i]=1;   /* copy info back to every gene in clone_cell */
    }
}

void reproduce_mut_koff(Genotype *genotype, Mutation *mut_record)
{   
    int tf_id,protein_id,i;
    char *tf_seq, *tf_seq_rc, *temp1, *temp2;
    
    tf_id=mut_record->pos_g;
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
        genotype->koff[genotype->nproteins-2]=mut_record->kinetic_diff;  /* update_protein_pool increases nprotein, 
                                                          * so the new protein will always be at nproteins-2.*/
    }    
    else
        genotype->koff[protein_id]=mut_record->kinetic_diff;    
    
    
#if SET_BS_MANUALLY 
    int i,j;
    for(i=1;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {
            if(genotype->all_binding_sites[i][j].tf_id==mut_record->pos_g)
                genotype->all_binding_sites[i][j].Koff=genotype->koff[genotype->which_protein[mut_record->pos_g]]*
                        pow(KR,(float)genotype->all_binding_sites[i][j].mis_match/3.0);
        }
    }
#endif
    for(i=0;i<genotype->ngenes;i++) 
    {
        genotype->recalc_TFBS[i]=1;   /* recalculate binding sites on every promoter */
        genotype->clone_info[i]=1;   /* copy info back to every gene in clone_cell */
    }
}

int mutate(Genotype *genotype, float kdis[NUM_K_DISASSEMBLY],RngStream RS, Mutation *mut_record)
{
    int i;
    int return_val=0;
    
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
    
    for(i=0;i<3;i++)
	mut_record->nuc_diff[i]='\0';
    switch (mut_record->mut_type)
    {
        case 's': //substitution in cis-reg       		
            mut_substitution(genotype,mut_record,RS);			
            return_val=MUT_CISSEQ;
            break;
        		
        case 'i': // insertion in cis-reg       
            mut_insertion(genotype,mut_record,RS);	
            return_val=MUT_CISSEQ;
            break;	
        		
        case 'p': // small deletion in cis-reg       
            mut_partial_deletion(genotype,mut_record,RS);
            return_val=MUT_CISSEQ;
            break;			
        		
        case 'w': // whole gene deletion.          
            mut_whole_gene_deletion(genotype,mut_record,RS);
            return_val=MUT_N_GENE;
            break;
        		
        case 'd': // Whole gene duplication also has two versions                  
            mut_duplication(genotype,mut_record,RS);
            return_val=MUT_N_GENE;
            break;
        
        case 'c': //binding sequence 
            mut_binding_sequence(genotype,mut_record,RS);
            return_val=MUT_TFSEQ;
            break;  
            
        case 'k': //mutations to kinetic constants        
            mut_kinetic_constant(genotype, mut_record,kdis,RS);
            return_val=MUT_KCONST;
            break;
            
        case 'e': //activator to repressor or the reverse
            mut_identity(genotype, mut_record, RS);
            return_val=MUT_TFSEQ;
            break; 
       
        case 'f': //mutations to the koff of a tf
            mut_koff(genotype,mut_record,RS);
            return_val=MUT_TFSEQ;
            break;
           
    }
    return return_val;
}

/* this function perform mutation indicated by input. Used to reproduce the genotype following a serial of mutations*/
int reproduce_mutate(Genotype *genotype, Mutation *mut_record)
{    
    int return_val=0;
    switch (mut_record->mut_type)
    {
        case 's': //substitution        		
            reproduce_substitution(genotype,mut_record);			
            return_val=MUT_CISSEQ;
            break;
        		
        case 'i': // insertion        
            reproduce_insertion(genotype,mut_record);	
            return_val=MUT_CISSEQ;
            break;	
        		
        case 'p': // partial deletion        
            reproduce_partial_deletion(genotype,mut_record);
            return_val=MUT_CISSEQ;
            break;			
        		
        case 'w': // whole gene deletion.         
            reproduce_whole_gene_deletion(genotype,mut_record);
            return_val=MUT_N_GENE;
            break;
        		
        case 'd': // Whole gene duplication also has two versions                  
            reproduce_gene_duplication(genotype,mut_record);
            return_val=MUT_N_GENE;
            break;
        
        case 'c': //binding sequence        
            reproduce_mut_binding_sequence(genotype,mut_record);            
            return_val=MUT_TFSEQ;
            break;  
            
        case 'k': //mutations in kinetic constants        
            reproduce_mut_kinetic_constant(genotype, mut_record);
            return_val=MUT_KCONST;
            break;    
            
        case 'e': //changing the identity of a TF
            reproduce_mut_identity(genotype, mut_record);
            return_val=MUT_TFSEQ;
            break; 
            
        case 'f':
            reproduce_mut_koff(genotype,mut_record);
            return_val=MUT_TFSEQ;
            break;  
    }
    return return_val;
}

/* this function calculates the probability of different mutations based on
 * the current genotype. It then modify the value of mut_type
 */
void draw_mutation(Genotype *genotype, char *mut_type, RngStream RS)
{
    float random;
    float tot_mut_rate=0.0;
    float tot_subs_rate, tot_indel_rate, tot_dup_rate, tot_sil_rate, tot_mut_kin_rate, tot_mut_identity_rate, tot_mut_binding_seq_rate, tot_mut_koff_rate; 
    int i;    
    
    /* duplication rate*/
    /* To restrict the duplication of selection genes, we allow each 
       selection gene to have at most 2 copies.*/  
    if(genotype->ngenes<NGENES)
    {
    #if RdcPdup    /* the method of reducing probability*/    
        tot_dup_rate=0.0;
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(genotype->protein_pool[i][0][0]<=MAX_COPIES)
                genotype->P_dup_per_protein[i]=genotype->protein_pool[i][0][0]*DUPLICATION;
            else
                genotype->P_dup_per_protein[i]=genotype->protein_pool[i][0][0]*DUPLICATION/reduction_in_P_dup;
            tot_dup_rate+=genotype->P_dup_per_protein[i];        
        }
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
            genotype->P_dup_per_protein[i]/=tot_dup_rate;
        tot_mut_rate+=tot_dup_rate;    
    #else    /* the method of increaing penalty*/    
        tot_dup_rate=(genotype->ngenes-N_SIGNAL_TF)*DUPLICATION;
        tot_mut_rate+=tot_dup_rate;
    #endif
    }
    else
    {
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
            genotype->P_dup_per_protein[i]=0.0;
        tot_dup_rate=0.0;
    }
    
    /* silencing rate*/
    tot_sil_rate=(float)(genotype->ngenes-1-N_SIGNAL_TF)*SILENCING; /* the signaling gene and one copy of the selection gene cannot be deleted */
    tot_mut_rate+=tot_sil_rate;
    
    /* calc total susbtitution rate*/
    tot_subs_rate=(float)(genotype->ngenes-N_SIGNAL_TF)*CISREG_LEN*SUBSTITUTION; /* NA to the signal gene*/
    tot_mut_rate+=tot_subs_rate;
    
    /* indel rate*/
    tot_indel_rate=(float)(genotype->ngenes-N_SIGNAL_TF)*CISREG_LEN*INDEL;    /* NA to the signal gene*/
    tot_mut_rate+=tot_indel_rate;
    
    /* mut in kinetic constants */
    tot_mut_kin_rate=(float)(genotype->ngenes-N_SIGNAL_TF)*MUTKINETIC*(1.0-proportion_mut_binding_seq-proportion_mut_identity-proportion_mut_koff); /* NA to the signal gene*/
    tot_mut_rate+=tot_mut_kin_rate;
    
    /* mut in binding seq*/
    tot_mut_binding_seq_rate=(float)(genotype->ntfgenes)*MUTKINETIC*proportion_mut_binding_seq; /* NA to the selection genes*/
    tot_mut_rate+=tot_mut_binding_seq_rate;
    
    /* mut in identity*/
    tot_mut_identity_rate=(float)(genotype->ntfgenes-2)*MUTKINETIC*proportion_mut_identity; /* if there's only one signal tf, then this tf must be an activator*/

    tot_mut_rate+=tot_mut_identity_rate;
    
    /* mut in koff*/
    tot_mut_koff_rate=(float)(genotype->ntfgenes)*MUTKINETIC*proportion_mut_koff;  
    tot_mut_rate+=tot_mut_koff_rate;
    
    
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

void update_protein_pool(Genotype *genotype, int protein_id, int gene_id, char mut_type)
{
    int i, j, protein_id_copy, gene_id_copy;
    
    switch (mut_type)
    {
        case 'w':/*a gene deletion*/        
            if(genotype->protein_pool[protein_id][0][0]==1) /* if this is the only copy */
            {   /* then we need to remove protein from protein_pool*/ 
                protein_id_copy=protein_id;
                for(i=0;i<genotype->nproteins-protein_id;i++)  
                {         
                    for(j=0;j<genotype->protein_pool[protein_id_copy][0][0];j++)
                        genotype->protein_pool[protein_id_copy][1][j]=-1;
                    genotype->protein_pool[protein_id_copy][0][0]=genotype->protein_pool[protein_id_copy+1][0][0];
                    for(j=0;j<genotype->protein_pool[protein_id_copy][0][0];j++)
                    {
                        gene_id_copy=genotype->protein_pool[protein_id_copy+1][1][j];
                        genotype->protein_pool[protein_id_copy][1][j]=(gene_id_copy>gene_id)?gene_id_copy-1:gene_id_copy;/*note that deletion changes the ids of genes!!!*/
                    }            
                    protein_id_copy++;
                }

                if(protein_id<genotype->nproteins-1) /*if a tf PROTEIN is deleted*/
                {   /* reduce the number of activator if necessary */ 
                    if(genotype->activating[protein_id] && (genotype->N_act!=0)) 
                        genotype->N_act--;
                    if(!genotype->activating[protein_id] && (genotype->N_rep!=0))
                        genotype->N_rep--;
                    /* also remove it from activating and koff */
                    protein_id_copy=protein_id;
                    for(i=0;i<genotype->nproteins-protein_id-1;i++)
                    {
                        genotype->activating[protein_id_copy]=genotype->activating[protein_id_copy+1];
                        genotype->koff[protein_id_copy]=genotype->koff[protein_id_copy+1];
                        protein_id_copy++;
                    }                   
                    /* in the case, all genes need to recalc binding sites*/
                    for(i=N_SIGNAL_TF;i<gene_id;i++)
                    {
                        genotype->recalc_TFBS[i]=1; /* recalc BS */
                        genotype->clone_info[i]=1; /* copy back the original BS in clone_cell*/
                    }                    
                }                
                /* update protein_id for gene<gene_id in which_protein*/
                for(i=N_SIGNAL_TF;i<gene_id;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i]<protein_id)?genotype->which_protein[i]:genotype->which_protein[i]-1;
                /* shift and update protein_id for gene>=gene_id in which_protein*/                
                for(i=gene_id;i<genotype->ngenes;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i+1]>protein_id)?genotype->which_protein[i+1]-1:genotype->which_protein[i+1]; /*the deletion also changes the ids of proteins*/
                  
                genotype->nproteins--;
            }  
            else /*if the protein has more than one genes*/
            {
                i=0;
                while(genotype->protein_pool[protein_id][1][i]!=gene_id) i++; /* find where is this gene_id*/
                j=i;
                for(;i<genotype->protein_pool[protein_id][0][0];i++)
                    genotype->protein_pool[protein_id][1][i]=(genotype->protein_pool[protein_id][1][i+1]>gene_id)?genotype->protein_pool[protein_id][1][i+1]-1:genotype->protein_pool[protein_id][1][i+1]; /*note that deletion changes the ids of genes!!!*/
                for(i=0;i<j;i++)
                    genotype->protein_pool[protein_id][1][i]=(genotype->protein_pool[protein_id][1][i]>gene_id)?genotype->protein_pool[protein_id][1][i]-1:genotype->protein_pool[protein_id][1][i];
                            
                genotype->protein_pool[protein_id][0][0]--;                
                /*take care of protein>protein_id*/                
                for(i=protein_id+1;i<genotype->nproteins;i++)  
                {   
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                        genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]>gene_id)?genotype->protein_pool[i][1][j]-1:genotype->protein_pool[i][1][j];/*note that deletion changes the ids of genes!!!*/
                }
                /*shift protein id for gene>=gene_id in which_protein*/
                for(i=gene_id;i<genotype->ngenes;i++)
                    genotype->which_protein[i]=genotype->which_protein[i+1];                
            }
            /*take care of protein<protein_id*/
            for(i=N_SIGNAL_TF;i<protein_id;i++)
            {
                for(j=0;j<genotype->protein_pool[i][0][0];j++)
                    genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]<gene_id)?genotype->protein_pool[i][1][j]:genotype->protein_pool[i][1][j]-1;
            }            
            break;        

        case 'd': /*a gene duplication*/
            /* add it to protein_pool, but do not change nproteins*/    
            genotype->protein_pool[protein_id][1][genotype->protein_pool[protein_id][0][0]]=genotype->ngenes-1; /*the newly duplicated gene takes the original place of the selection gene*/
            genotype->protein_pool[protein_id][0][0]++; 

            /*update the id of the original selection gene stored in protein_pool*/           
            j=0;
            while(genotype->protein_pool[genotype->nproteins-1][1][j]!=genotype->ngenes-1)j++;
            genotype->protein_pool[genotype->nproteins-1][1][j]++;                
            
            /*update which_protein*/          
            genotype->which_protein[genotype->ngenes]=genotype->which_protein[genotype->ngenes-1];            
            genotype->which_protein[genotype->ngenes-1]=protein_id;            
            break;

        case 'c': /*mutation in tf binding seq*/
            /* remove this copy of gene from the original protein*/
            i=0;
            while(genotype->protein_pool[protein_id][1][i]!=gene_id) i++;
            for(;i<genotype->protein_pool[protein_id][0][0];i++) 
            {
                genotype->protein_pool[protein_id][1][i]= genotype->protein_pool[protein_id][1][i+1]; /* rearrange data array */
            }
            genotype->protein_pool[protein_id][0][0]--;
            /* increase the protein id of selection protein*/ 
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j];
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1;
            }                       
            /* create a new protein and link it to this gene*/
            genotype->which_protein[gene_id]=genotype->nproteins-1; /*put the new protein to the pos of the first selection gene*/
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=gene_id;
            /* update activating*/
            if(genotype->activating[protein_id]) /* increase the number of activator */
                genotype->N_act++;
            else
                genotype->N_rep++;

            genotype->activating[genotype->nproteins-1]=genotype->activating[protein_id];
            /* update koff*/
            genotype->koff[genotype->nproteins-1]=genotype->koff[protein_id];
            /* finally, update protein numbers*/
            genotype->nproteins++;
            /* NOTE: this mutation does not change the number of genes*/
            break;
            
        case 'e': /*mutation in the identity of a TF*/
            /* remove this copy of gene from the original protein*/
            i=0;
            while(genotype->protein_pool[protein_id][1][i]!=gene_id) i++;
            for(;i<genotype->protein_pool[protein_id][0][0];i++) 
            {
                genotype->protein_pool[protein_id][1][i]= genotype->protein_pool[protein_id][1][i+1]; /* rearrange data array */
            }
            genotype->protein_pool[protein_id][0][0]--;
            /* increase the protein id of selection proteins*/               
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j];
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1;
            }            
            /* create a new protein and link it to this gene*/
            genotype->which_protein[gene_id]=genotype->nproteins-1; /*put the new protein to the pos of the first selection gene*/
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=gene_id;
            /* update activating*/
            if(genotype->activating[protein_id]) 
                genotype->N_rep++;  /* an activator turns into a repressor */
            else
                genotype->N_act++;

            genotype->activating[genotype->nproteins-1]=genotype->activating[protein_id];
            /* update koff*/
            genotype->koff[genotype->nproteins-1]=genotype->koff[protein_id];
            /* finally, update protein numbers*/
            genotype->nproteins++;
            /* NOTE: this mutation does not change the number of genes*/
            break;
            
        case 'f': /*mutation in tf koff*/
            /* remove this copy of gene from the original protein pool*/
            i=0;
            while(genotype->protein_pool[protein_id][1][i]!=gene_id) i++;
            for(;i<genotype->protein_pool[protein_id][0][0];i++) 
            {
                genotype->protein_pool[protein_id][1][i]= genotype->protein_pool[protein_id][1][i+1]; /* rearrange data array */
            }
            genotype->protein_pool[protein_id][0][0]--;
            /* increase the protein id of selection proteins*/               
            genotype->protein_pool[genotype->nproteins][0][0]=genotype->protein_pool[genotype->nproteins-1][0][0];
            for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)                
            {
                genotype->protein_pool[genotype->nproteins][1][j]=genotype->protein_pool[genotype->nproteins-1][1][j];
                genotype->which_protein[genotype->protein_pool[genotype->nproteins-1][1][j]]++;
                genotype->protein_pool[genotype->nproteins-1][1][j]=-1;
            }            
            /* create a new protein and link it to this gene*/
            genotype->which_protein[gene_id]=genotype->nproteins-1; /*put the new protein to the pos of the first selection gene*/
            genotype->protein_pool[genotype->nproteins-1][0][0]=1;
            genotype->protein_pool[genotype->nproteins-1][1][0]=gene_id;            
            /* update activating*/
            if(genotype->activating[protein_id]) 
                genotype->N_act++;  
            else
                genotype->N_rep++;

            genotype->activating[genotype->nproteins-1]=genotype->activating[protein_id];  
            /* finally, update protein numbers*/
            genotype->nproteins++;
            /* NOTE: this mutation does not change the number of genes*/
            break;
    }
}

void update_cisreg_cluster(Genotype *genotype, int gene_id, char mut_type)
{
    int cluster_id, cluster_id_copy, i, j;   
    
    cluster_id=genotype->which_cluster[gene_id];   
    
    if(mut_type!='d')/*not a duplication*/
    {   
        if(mut_type!='w') /*not a gene deletion*/
        {
            if(genotype->cisreg_cluster[cluster_id][1]!=-1) /*if it is not the only one of its kind*/
            {   /*then remove pos_g from the original cluster*/        
                i=0;
                while(genotype->cisreg_cluster[cluster_id][i]!=gene_id) i++;
                while(genotype->cisreg_cluster[cluster_id][i]!=-1)
                {
                    genotype->cisreg_cluster[cluster_id][i]=genotype->cisreg_cluster[cluster_id][i+1];
                    i++;
                }                     
                /* and create a new cluster*/
                i=N_SIGNAL_TF; // start from the first non-signal-tf gene
                while(genotype->cisreg_cluster[i][0]!=-1) i++; 
                genotype->cisreg_cluster[i][0]=gene_id;
                genotype->which_cluster[gene_id]=i;
            }
        }
        else /*is gene deletion*/
        {               
            if(genotype->cisreg_cluster[cluster_id][1]!=-1) /*if it is not the only one of its kind*/
            {   /*then remove pos_g from the original cluster*/        
                i=0;
                while(genotype->cisreg_cluster[cluster_id][i]!=gene_id) i++;
                while(genotype->cisreg_cluster[cluster_id][i]!=-1)
                {                    
                    genotype->cisreg_cluster[cluster_id][i]=((genotype->cisreg_cluster[cluster_id][i+1]-1)<-1)?-1:(genotype->cisreg_cluster[cluster_id][i+1]-1); /*gene id changed due to deletion*/
                    i++;
                }
                /*take care of clusters>cluster_id*/
                cluster_id_copy=cluster_id+1;
                while(genotype->cisreg_cluster[cluster_id_copy][0]!=-1)
                {
                    i=0;
                    while(genotype->cisreg_cluster[cluster_id_copy][i]!=-1)
                    {
                        genotype->cisreg_cluster[cluster_id_copy][i]=(genotype->cisreg_cluster[cluster_id_copy][i]<gene_id)?genotype->cisreg_cluster[cluster_id_copy][i]:genotype->cisreg_cluster[cluster_id_copy][i]-1;
                        i++;
                    }
                    cluster_id_copy++;
                }
                /*update which_cluster*/
                for(i=gene_id;i<genotype->ngenes;i++)                              
                    genotype->which_cluster[i]=genotype->which_cluster[i+1];
            }
            else
            {   /*need to shift cisreg_cluster*/ 
                cluster_id_copy=cluster_id;
                while(genotype->cisreg_cluster[cluster_id_copy][0]!=-1)
                {
                    /*reset cluster=cluster_id*/
                    i=0;
                    while(genotype->cisreg_cluster[cluster_id_copy][i]!=-1)
                    {
                        genotype->cisreg_cluster[cluster_id_copy][i]=-1;
                        i++;
                    }
                    /*then copy from cluster=cluster_id+1*/                    
                    i=0;
                    while(genotype->cisreg_cluster[cluster_id_copy+1][i]!=-1)
                    {                        
                        genotype->cisreg_cluster[cluster_id_copy][i]=(genotype->cisreg_cluster[cluster_id_copy+1][i]<gene_id)?genotype->cisreg_cluster[cluster_id_copy+1][i]:genotype->cisreg_cluster[cluster_id_copy+1][i]-1; /*note deletion changes gene id*/
                        i++;
                    }
                    cluster_id_copy++;
                }
                /*update which_cluster for gene id>=gene_id*/                
                for(i=gene_id;i<genotype->ngenes-1;i++)                                
                    genotype->which_cluster[i]=(genotype->which_cluster[i+1]<cluster_id)?genotype->which_cluster[i+1]:genotype->which_cluster[i+1]-1;
                /*take care of which_cluster for gene id<gene_id*/
                for(i=N_SIGNAL_TF;i<gene_id;i++)
                    genotype->which_cluster[i]=(genotype->which_cluster[i]>cluster_id)?genotype->which_cluster[i]-1:genotype->which_cluster[i];
            }
            /*take care of clusters<cluster_id*/
            for(i=N_SIGNAL_TF;i<cluster_id;i++)
            {
                j=0;
                while(genotype->cisreg_cluster[i][j]!=-1)
                {
                    genotype->cisreg_cluster[i][j]=(genotype->cisreg_cluster[i][j]<gene_id)?genotype->cisreg_cluster[i][j]:genotype->cisreg_cluster[i][j]-1;
                    j++;
                }
            }                                            
        }
    }
    else /*is gene duplication*/
    {    
        i=0;
        while(genotype->cisreg_cluster[cluster_id][i]!=-1)i++;        
        if(genotype->cisreg_cluster[cluster_id][i-1]!=genotype->ngenes-1) /*not the special case that a selection gene is duplicated*/
        {    
            i=0;
            while(genotype->cisreg_cluster[cluster_id][i]!=-1)i++;
            genotype->cisreg_cluster[cluster_id][i]=genotype->ngenes-1;
            
            /* now increase the id of the selection gene in cisreg_cluster*/           
            cluster_id_copy=genotype->which_cluster[genotype->ngenes-1];/*Note that ngenes has not been updated when update_cisreg_cluster is called*/
            i=0;
            while(genotype->cisreg_cluster[cluster_id_copy][i]!=genotype->ngenes-1)i++;
            genotype->cisreg_cluster[cluster_id_copy][i]=genotype->ngenes;
        }
        else
        {                     
            i=0;
            while(genotype->cisreg_cluster[cluster_id][i]!=-1)i++;
            genotype->cisreg_cluster[cluster_id][i]=genotype->ngenes;
            genotype->cisreg_cluster[cluster_id][i-1]=genotype->ngenes-1;                      
        }
        /*update which_cluster*/
        genotype->which_cluster[genotype->ngenes]=genotype->which_cluster[genotype->ngenes-1];
        genotype->which_cluster[genotype->ngenes-1]=cluster_id;
    }   
}

/****************** end of mutation functions *********************************/

void initialize_cache(Genotype *genotype)
{
    int j,k;
    FILE *fp;    
    
    for(j=0;j<NGENES;j++)
    {
        genotype->which_protein[j]=-1;
        genotype->clone_info[j]=1;          
        genotype->recalc_TFBS[j]=1;
        genotype->which_cluster[j]=-1; 
        genotype->min_act_to_transc[j]=MAX_MODE+1; /* by default a gene cannot be turned on*/
        genotype->koff[j]=-1.0;
        for(k=0;k<NGENES;k++)        
            genotype->cisreg_cluster[j][k]=-1;
    }    
    /* alloc space for protein_pool */
    for(j=0;j<NPROTEINS;j++)
    {
        genotype->protein_pool[j][0][0]=0;
        for(k=0;k<NGENES;k++)
        {
            genotype->protein_pool[j][1][k]=-1;
        }
        genotype->activating[j]=-1;
    }
    /* alloc space for binding sites*/
    for(j=0;j<NGENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(MAXELEMENTS*sizeof(AllTFBindingSites)); 
//        genotype->compressed_binding_sites[j]=malloc(MAXELEMENTS*sizeof(CompressedBindingSites));
//        genotype->N_configurations[j] = malloc(MAXELEMENTS*sizeof(int));
        
        if (!(genotype->all_binding_sites[j])) 
        {
            fp=fopen("output.txt","a+");
            fprintf(fp,"error in initialize_cache\n");
            fclose(fp);
            printf("error in initialize_cache\n");
            exit(1);
        }
    }
    for(j=N_SIGNAL_TF;j<NGENES;j++)
    {
        genotype->N_act_BS[j]=0;
        genotype->N_rep_BS[j]=0;
        genotype->binding_sites_num[j]=0;
    }    
}

void try_fixation(Genotype *genotype_ori, Genotype *genotype_offspring, int *fixation, float *pfix, RngStream RS)
{    
    float v1, v2;
    v1=(genotype_offspring->var_GR1*env1_occurence+genotype_offspring->var_GR2*env2_occurence)/(env1_occurence+env2_occurence);
    v2=(genotype_ori->var_GR1*env1_occurence+genotype_ori->var_GR2*env2_occurence)/(env1_occurence+env2_occurence);
//    v1=genotype_offspring->var_GR1*genotype_offspring->var_GR2+genotype_offspring->var_GR1*genotype_offspring->avg_GR2*genotype_offspring->avg_GR2+genotype_offspring->var_GR2*genotype_offspring->avg_GR1*genotype_offspring->avg_GR1;
//    v2=genotype_ori->var_GR1*genotype_ori->var_GR2+genotype_ori->var_GR1*genotype_ori->avg_GR2*genotype_ori->avg_GR2+genotype_ori->var_GR2*genotype_ori->avg_GR1*genotype_ori->avg_GR1;
        
    *pfix=qz(genotype_offspring->fitness, v1, genotype_ori->fitness, v2);
    
    if(*pfix>RngStream_RandU01(RS))
        *fixation=1;
    else
        *fixation=0;
}

void replay_mutations(  Genotype *genotype_ori,
                        Genotype *genotype_ori_copy,
                        FILE *file_mutation,
                        Mutation *mut_record,
                        int replay_N_steps)
{
    int i;
    
    for(i=0;i<replay_N_steps;i++)
    {
        clone_cell_forward(genotype_ori,genotype_ori_copy,COPY_ALL);
        fscanf(file_mutation,"%c %d %d %s %d %f\n",
                &(mut_record->mut_type),
                &(mut_record->pos_g),
                &(mut_record->pos_n),
                mut_record->nuc_diff,
                &(mut_record->kinetic_type),
                &(mut_record->kinetic_diff));
        reproduce_mutate(genotype_ori_copy,mut_record);            
        clone_cell_forward(genotype_ori_copy,genotype_ori,COPY_ALL);
        printf("%d, %c, %f\n",i,mut_record->mut_type,genotype_ori->fitness);
    }
}

#if NEUTRAL
void evolve_neutrally(  Genotype *genotype_ori,
                        Genotype *genotype_ori_copy,
                        char *file_mutations,
                        Mutation *mut_record,
                        RngStream RS_main)
{
    int i=0;
    FILE *fp;
    while(genotype_ori.ntfgenes>2)
    {
        while(1)
        {
            clone_cell_forward(genotype_ori,genotype_ori_copy,COPY_ALL);
            mutate(genotype_ori_copy,kdis,RS_main,mut_record);
            if(RngStream_RandU01(RS_main)<=0.5)
            {
                fp=fopen(file_mutations,"a+");
                fprintf(fp,"%c %d %d %s %d %f\n",
                        mut_record->mut_type,    
                        mut_record->pos_g,
                        mut_record->pos_n,
                        mut_record->nuc_diff,
                        mut_record->kinetic_type,
                        mut_record->kinetic_diff);
                fclose(fp);
                break;
            }
        }
        clone_cell_forward(genotype_ori_copy,genotype_ori,COPY_ALL);   
        i++;
        if(i%100==0 || genotype_ori.nproteins==4)
        {
            calc_all_binding_sites(genotype_ori);
            summarize_binding_sites(genotype_ori,i);
        }
    } 
}
#endif

int evolve_N_steps( Genotype *genotype_ori, 
                    Genotype *genotype_ori_copy,
                    int init_step, 
                    int steps, 
                    int *N_tot_trials,
                    int recalc_new_fitness, 
                    int maxbound2, 
                    int maxbound3,
                    char *file_genotype_summary,
                    char *file_mutations,
                    char *file_errors,
                    char *file_act_BS, 
                    char *file_rep_BS, 
                    char *file_all_BS,   
                    float init_protein_number[NPROTEINS],
                    int init_mRNA[NGENES], 
                    float kdis[NUM_K_DISASSEMBLY],
                    Mutation *mut_record, 
                    RngStream RS_main,
                    RngStream RS_parallel[N_THREADS],
                    int burn_in)
{
    int i,j;    
    int ready_to_evolve;
    int fixation;
    int N_trials;
    int next_step=0;
    float avg_GR1_copy, avg_GR2_copy, var_GR1_copy, var_GR2_copy;
    float pfix;  
    FILE *fp;
    
    if(burn_in)
        ready_to_evolve=0;
    
    for(i=init_step;i<steps;i++)
    {             
        fixation=0;      
        N_trials=0;
       
        if(genotype_ori->ntfgenes==2)
        { 
            fp=fopen(file_genotype_summary,"a+");                              
            fprintf(fp,"ntfgenes critical!\n");
            fclose(fp); 
            break;
        } 

        while(1)
        {	
            N_trials++;
            (*N_tot_trials)++;
            clone_cell_forward(genotype_ori,genotype_ori_copy,COPY_ALL);	
            mutate(genotype_ori_copy,kdis,RS_main,mut_record);
            fp=fopen("MUT_Detail.txt","a+");
            fprintf(fp,"%d %d %c %d %d %s %d %f\n",
                    i,
                    *N_tot_trials,
                    mut_record->mut_type,
                    mut_record->pos_g,
                    mut_record->pos_n,
                    mut_record->nuc_diff,
                    mut_record->kinetic_type,
                    mut_record->kinetic_diff);
            fclose(fp);
            
        #if !SET_BS_MANUALLY		
            calc_all_binding_sites(genotype_ori_copy);
        #endif

            calc_avg_growth_rate(   genotype_ori_copy,
                                    init_mRNA,
                                    init_protein_number,
                                    maxbound2,  
                                    maxbound3,
                                    RS_parallel,
                                    file_genotype_summary,
                                    file_errors,
                                    i,
                                    mut_record);
            try_fixation(genotype_ori,genotype_ori_copy,&fixation,&pfix,RS_main); 
            if(fixation==1)
            {                    
                fp=fopen(file_genotype_summary,"a+");
                fprintf(fp,"step %d %d %d %c %.10f ",i, *N_tot_trials, N_trials,mut_record->mut_type,pfix);
                fclose(fp);
                fp=fopen(file_mutations,"a+");
                fprintf(fp,"%c %d %d %s %d %f\n",
                        mut_record->mut_type,    
                        mut_record->pos_g,
                        mut_record->pos_n,
                        mut_record->nuc_diff,
                        mut_record->kinetic_type,
                        mut_record->kinetic_diff);
                fclose(fp);
                break;
            }
        }            
        clone_cell_forward(genotype_ori_copy,genotype_ori,COPY_ALL);
        
        #if !SET_BS_MANUALLY
            calc_all_binding_sites(genotype_ori); 
        #endif  
            
        /*increase the accuracy of the fitness of the new genotype*/                    
        for(j=1;j<=recalc_new_fitness;j++)
        {
            avg_GR1_copy=genotype_ori->avg_GR1;
            avg_GR2_copy=genotype_ori->avg_GR2;
            var_GR1_copy=genotype_ori->var_GR1;
            var_GR2_copy=genotype_ori->var_GR2;
            calc_avg_growth_rate(   genotype_ori, 
                                    init_mRNA,
                                    init_protein_number,
                                    maxbound2,
                                    maxbound3,
                                    RS_parallel,
                                    file_genotype_summary,
                                    file_errors,
                                    0,
                                    mut_record); 
            genotype_ori->avg_GR1=(float)(genotype_ori->avg_GR1+j*avg_GR1_copy)/(j+1);
            genotype_ori->avg_GR2=(float)(genotype_ori->avg_GR2+j*avg_GR2_copy)/(j+1);
            genotype_ori->var_GR1=(float)(genotype_ori->var_GR1+j*var_GR1_copy)/(j+1);
            genotype_ori->var_GR1=(float)(genotype_ori->var_GR2+j*var_GR2_copy)/(j+1);
            genotype_ori->fitness=(genotype_ori->avg_GR1*env1_occurence+genotype_ori->avg_GR2*env2_occurence)/(env1_occurence+env2_occurence);
        }            

        calc_all_binding_sites(genotype_ori); 
        
        if(i%100==0 && i!=0) /*output genotype every 100 steps*/
            summarize_binding_sites(genotype_ori,i);

        output_genotype(file_genotype_summary, file_act_BS, file_rep_BS, file_all_BS, genotype_ori, i);
        
        if(burn_in)
        {
            if(genotype_ori->avg_GR1>0.5*gpeak && genotype_ori->avg_GR2>0.5*gpeak)
                ready_to_evolve++;
            else
                ready_to_evolve=0;
            if(ready_to_evolve>=10)
            {
                next_step=i+1;
                break;
            }
        }
    }     
    return next_step;
}

int init_run_pop(   float kdis[NUM_K_DISASSEMBLY], 
                    char *RuntimeSumm, 
                    char *file_genotype_summary, 
                    char *file_mutations, 
                    char *file_errors, 
                    char *file_act_BS,                         
                    char *file_rep_BS, 
                    char *file_all_BS, 
                    unsigned long int seeds[6])
{  
    int i;
    Genotype genotype_ori;
    Genotype genotype_ori_copy;    
    int init_mRNA[NGENES]; 
    float init_protein_number[NGENES];
    int N_tot_trials=0; 
    int init_step=0;
    int maxbound2 = MAXBOUND;
    int maxbound3 = 10*MAXBOUND; 
    Mutation mut_record;
    FILE *MUT,*fp; 
    RngStream RS_main,RS_parallel[N_THREADS];      	
	
    /* initialize random number seeds*/
    RngStream_SetPackageSeed(seeds);    
    RS_main=RngStream_CreateStream("Main");
    for(i=0; i < N_THREADS; i++)
        RS_parallel[i]=RngStream_CreateStream("");
    
    /* set initial numbers of mRNA and protein*/
    for(i=N_SIGNAL_TF; i < NGENES; i++) /* loop through tf genes*/
    {
        init_mRNA[i]=0;
        init_protein_number[i]=10.0;
    }

    /* set koff */
    for(i=0;i<TF_ELEMENT_LEN-NMIN+1;i++)          
        Koff[i]=koff*pow(KR,(float)i/3.0);
    
    /*manually input binding sites info below*/
    #if SET_BS_MANUALLY 
    // this is a FFL
        int N_BS_per_gene[3]={2,2,2};
        int tf_id_per_site[3][3]={{0,0,-1},{0,2,2},{1,1,-1}};
        int hindrance_per_site[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        float koff_per_site[3][3]={{Koff[0],Koff[0],Koff[0]},{Koff[0],Koff[0],Koff[0]},{Koff[0],Koff[0],Koff[0]}};    

    // this is a BRA    
    //    int N_BS_per_gene[3]={2,3,2};
    //    int tf_id_per_site[3][3]={{1,1,-1},{0,2,0},{1,1,-1}};
    //    int hindrance_per_site[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    //    float koff_per_site[3][3]={{Koff[0],Koff[0],Koff[0]},{Koff[0],Koff[0],Koff[0]},{Koff[0],Koff[0],Koff[0]}};    
    #endif    
    
    /* initialize genotype */
    initialize_cache(&genotype_ori);
    initialize_cache(&genotype_ori_copy);    
    initialize_genotype(&genotype_ori, kdis, RS_main);    
#if SET_BS_MANUALLY
    set_binding_sites(&genotype_ori,N_BS_per_gene,tf_id_per_site,hindrance_per_site,koff_per_site);
#endif    
    genotype_ori_copy.ngenes=genotype_ori.ngenes;
    genotype_ori_copy.ntfgenes=genotype_ori.ntfgenes;
    genotype_ori_copy.nproteins=genotype_ori.nproteins; 
    
    /* record the initial network topology*/
    summarize_binding_sites(&genotype_ori,init_step); /*snapshot of the initial (0) distribution binding sites */   
    
    /* initialize mut_record */
    mut_record.kinetic_diff=0.0;
    mut_record.kinetic_type=-1;
    mut_record.mut_type='\0';
    mut_record.nuc_diff[0]='\0';
    mut_record.nuc_diff[1]='\0';
    mut_record.nuc_diff[2]='\0';
    mut_record.pos_g=-1;
    mut_record.pos_n=-1;
    
    /* plot the phenotype and fitness of a genotype */
    #if PLOTTING        
        fp=fopen("MUT_1.txt","r");    
        if(MUT!=NULL)
        {
            int replay_N_steps;
            printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
            replay_N_steps=21;
            replay_mutations(&genotype_ori, &genotype_ori_copy, fp, &mut_record, replay_N_steps);       
            fclose(fp);   
            calc_all_binding_sites(&genotype_ori);        
            summarize_binding_sites(&genotype_ori,1);   
            
            /* conditions under which the phenotype and fitness is measured */
            init_env1='A';
            init_env2='A'; 
            tdevelopment = 149.9;
            duration_of_burn_in_growth_rate = 5.0; 
            env1_t_signalA=180.0;     
            env1_t_signalB=0.0;
            env2_t_signalA=10.0;
            env2_t_signalB=40.0;
            env1_signalA_as_noise=0;    
            env2_signalA_as_noise=1;       
            N_replicates=400;
            env1_occurence=0.5;
            env2_occurence=0.5;
            
            calc_avg_growth_rate_plotting(  &genotype_ori, 
                                            init_mRNA, 
                                            init_protein_number, 
                                            maxbound2, 
                                            maxbound3, 
                                            RS_parallel);
        }    
    #else
        /* summarize the initial genotype */ 
        tdevelopment=149.9;
        duration_of_burn_in_growth_rate=5.0;
        init_env1='A';
        init_env2='B'; 
        env1_t_signalA=180.0;     
        env1_t_signalB=0.0;
        env2_t_signalA=0.0;
        env2_t_signalB=180.0;
        env1_signalA_as_noise=0;    
        env2_signalA_as_noise=1;     
        N_replicates=120;
        recalc_new_fitness=1;
        env1_occurence=0.5;
        env2_occurence=0.5;
        calc_avg_growth_rate(   &genotype_ori, 
                                init_mRNA,
                                init_protein_number,
                                maxbound2,
                                maxbound3,
                                RS_parallel,
                                file_genotype_summary,
                                file_errors,
                                0,
                                &mut_record);
        
	fp=fopen(file_genotype_summary,"a+");
        fprintf(fp,"intial: %.10f %.10f %.10f %.10f %d %d %d %d \n",            
                genotype_ori.avg_GR1,
                genotype_ori.avg_GR2,
                sqrt(genotype_ori.var_GR1)/genotype_ori.avg_GR1,
                sqrt(genotype_ori.var_GR2)/genotype_ori.avg_GR2,
                genotype_ori.ngenes,
                genotype_ori.nproteins,
                genotype_ori.N_act,
                genotype_ori.N_rep);
	fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step fixed_mutation Pfix avg_GR1 avg_GR2 cv_GR1 cv_GR2 N_genes N_proteins N_activator N_repressor \n");	
	fclose(fp);	
    #endif        

    /* Run simulation burn_in */
    #if BURN_IN && !NEUTRAL && !SUDO_REPLICATES
        /* burn_in conditions*/
        tdevelopment = 149.9;
        duration_of_burn_in_growth_rate = 5.0; 
        init_env1='A';
        init_env2='B';
        env1_t_signalA=180.0;    
        env1_t_signalB=0.0;     
        env2_t_signalA=10.0;
        env2_t_signalB=40.0;
        env1_signalA_as_noise=0;    
        env2_signalA_as_noise=1;  
        env1_signalA_mismatches=0;    
        env2_signalA_mismatches=1; 
        N_replicates=1200;
        recalc_new_fitness=0;
        env1_occurence=0.5;
        env2_occurence=0.5;	
        init_step=0; 
        
        fp=fopen(RuntimeSumm,"a+");
        fprintf(fp,"**********Burn-in conditions**********\n");
        fprintf(fp,"BURN_IN=%d\n",BURN_IN);                
        fprintf(fp,"N_replicates=%d\n",N_replicates);        
        fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
        fprintf(fp,"T-development=%f\n",tdevelopment);
        fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);        
        fprintf(fp,"Test 1: initial signal=%c, Env1 duration=%f min, Env2 duration=%f min, signalA as noise=%d, occurrence=%f\n",init_env1,env1_t_signalA, env1_t_signalB=0.0, env1_signalA_as_noise, env1_occurence);
        fprintf(fp,"Test 2: initial signal=%c, Env1 duration=%f min, Env2 duration=%f min, signalA as noise=%d, occurrence=%f\n",init_env2,env2_t_signalA, env2_t_signalB=0.0, env2_signalA_as_noise, env2_occurence);
        fclose(fp); 
        
        /* run burn_in */        
        init_step=evolve_N_steps(   &genotype_ori, 
                                    &genotype_ori_copy,
                                    init_step, 
                                    BURN_IN, 
                                    &N_tot_trials,
                                    recalc_new_fitness, 
                                    maxbound2,
                                    maxbound3,
                                    file_genotype_summary,
                                    file_mutations,
                                    file_errors,
                                    file_act_BS, 
                                    file_rep_BS, 
                                    file_all_BS   
                                    init_protein_number,
                                    init_mRNA, 
                                    kdis,
                                    &mut_record, 
                                    RS_main,
                                    RS_parallel,
                                    1); // 1 means we are running burn_in
        
        fp=fopen(RuntimeSumm,"a");
        fprintf(fp,"Burn_in finished after %dth step\n",i);
        fclose(fp);        
    #endif
    /************************* End of Burn-in ******************************/
        

    /* Start with the genotype at the end of the burn_in in previous simulations*/
    #if SUDO_REPLICATES        
        MUT=fopen("previous_record.txt","r"); 
        int replay_N_steps;
        if(MUT!=NULL)
        {
            fp=fopen(RuntimeSumm,"a+");
            fprintf(fp,"LOAD MUTATION RECORD SUCCESSFUL!\n");            
            replay_N_steps=100;
            fprintf(fp,"Replay %d mutations\n",replay_N_steps);
            fclose(fp);
            replay_mutations(&genotype_ori, &genotype_ori_copy, MUT, &mut_record, replay_N_steps);
            fclose(MUT);
        }
        else
        {
            fp=fopen(RuntimeSumm,"a+");
            fprintf(fp,"Cannot find previous Burn_in file...\n");
            return 0;
        }         
        calc_all_binding_sites(&genotype_ori);
        calc_avg_growth_rate(   &genotype_ori, 
                                init_mRNA,
                                init_protein_number,
                                maxbound2,
                                maxbound3,
                                RS_parallel,
                                file_genotype_summary,
                                file_errors,
                                replay_N_steps,
                                &mut_record);
        fp=fopen(file_genotype_summary,"a+");
        fprintf(fp,"after burn_in: ");
        fclose(fp);
        output_genotype(file_genotype_summary, file_act_BS, file_rep_BS, file_all_BS, &genotype_ori, replay_N_steps); 
    #endif
       
    /* run full simulation */
    DUPLICATION=1.5e-6*0.33;                 /* per gene per cell division (using 120min), excluding chromosome duplication. Lynch 2008*/
    SILENCING = 1.3e-6*0.33;          /* per gene per cell division (120min), excluding chromosome deletion.Lynch 2008*/
    tdevelopment=149.9;
    duration_of_burn_in_growth_rate=5.0;
    init_env1='A';
    init_env2='B'; 
    env1_t_signalA=180.0;     
    env1_t_signalB=0.0;
    env2_t_signalA=0.0;
    env2_t_signalB=180.0;
    env1_signalA_as_noise=0;    
    env2_signalA_as_noise=1;     
    N_replicates=120;
    recalc_new_fitness=1;
    env1_occurence=0.5;
    env2_occurence=0.5;
        
    #if !NEUTRAL
        fp=fopen(RuntimeSumm,"a");
        fprintf(fp,"**********second phase conditions**********\n");
        fprintf(fp,"second phase steps=%d\n",MAX_MUT_STEP-BURN_IN);                
        fprintf(fp,"N_replicates=%d\n",N_replicates);        
        fprintf(fp,"N_recalc_fitness=%d\n",recalc_new_fitness);
        fprintf(fp,"T-development=%f\n",tdevelopment);
        fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);         
        fprintf(fp,"Test 1: initial signal=%c, Env1 duration=%f min, Env2 duration=%f min, signalA as noise=%d, occurrence=%f\n",init_env1,env1_t_signalA, env1_t_signalB, env1_signalA_as_noise, env1_occurence);
        fprintf(fp,"Test 2: initial signal=%c, Env1 duration=%f min, Env2 duration=%f min, signalA as noise=%d, occurrence=%f\n",init_env2,env2_t_signalA, env2_t_signalB, env2_signalA_as_noise, env2_occurence);        
        fclose(fp); 
        
        init_step=evolve_N_steps(   &genotype_ori, 
                                    &genotype_ori_copy,
                                    init_step, 
                                    MAX_MUT_STEP, 
                                    &N_tot_trials,
                                    recalc_new_fitness, 
                                    maxbound2,
                                    maxbound3,
                                    file_genotype_summary,
                                    file_mutations,
                                    file_errors,
                                    file_act_BS, 
                                    file_rep_BS, 
                                    file_all_BS,   
                                    init_protein_number,
                                    init_mRNA, 
                                    kdis,
                                    &mut_record, 
                                    RS_main,
                                    RS_parallel,
                                    0); // 0 means we are not running burn_in
    #else
        evolve_neutrally(   &genotype_ori,
                            &genotype_ori_copy,
                            file_mutations,
                            &mut_record,
                            RngStream RS_main);        
    #endif
    
    calc_all_binding_sites(&genotype_ori);
    summarize_binding_sites(&genotype_ori,init_step); /*snapshot of the final distribution binding sites */
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

void output_genotype(char *file_overview, char *file_act_BS, char *file_rep_BS, char *file_all_BS, Genotype *genotype, int step_i)
{
    FILE *OUTPUT;
    char filename[32];
    int j,k;
    
    OUTPUT=fopen(file_overview,"a+");
    fprintf(OUTPUT,"%.10f %.10f %.10f %.10f %d %d %d %d \n",            
            genotype->avg_GR1,
            genotype->avg_GR2,
            sqrt(genotype->var_GR1)/genotype->avg_GR1,
            sqrt(genotype->var_GR2)/genotype->avg_GR2,
            genotype->ngenes,
            genotype->nproteins,
            genotype->N_act,
            genotype->N_rep);
    fclose(OUTPUT);    
    OUTPUT=fopen(file_act_BS,"a+");	     
    fprintf(OUTPUT,"%d %d %d %d %d %d %d %d %d\n",
            step_i,
            genotype->N_act_BS[2],
            genotype->N_act_BS[3],
            genotype->N_act_BS[4],
            genotype->N_act_BS[5],
            genotype->N_act_BS[6],
            genotype->N_act_BS[7],
            genotype->N_act_BS[8],
            genotype->N_act_BS[9]);
    fclose(OUTPUT);
    OUTPUT=fopen(file_rep_BS,"a+");
    fprintf(OUTPUT,"%d %d %d %d %d %d %d %d %d\n",
            step_i,
            genotype->N_rep_BS[2],
            genotype->N_rep_BS[3],
            genotype->N_rep_BS[4],
            genotype->N_rep_BS[5],
            genotype->N_rep_BS[6],
            genotype->N_rep_BS[7],
            genotype->N_rep_BS[8],
            genotype->N_rep_BS[9]);
    fclose(OUTPUT);
    OUTPUT=fopen(file_all_BS,"a+");
    fprintf(OUTPUT,"%d %d %d %d %d %d %d %d %d\n",
            step_i,
            genotype->binding_sites_num[2],
            genotype->binding_sites_num[3],
            genotype->binding_sites_num[4],
            genotype->binding_sites_num[5],
            genotype->binding_sites_num[6],
            genotype->binding_sites_num[7],
            genotype->binding_sites_num[8],
            genotype->binding_sites_num[9]);
    fclose(OUTPUT);
    for(k=N_SIGNAL_TF;k<genotype->ngenes;k++)
    {
            snprintf(filename,sizeof(char)*32,"CIS_%i.txt",k);			
            OUTPUT=fopen(filename,"a+");
            fprintf(OUTPUT,"%d",step_i);
            for(j=0;j<CISREG_LEN;j++)
                    fprintf(OUTPUT,"%c",genotype->cisreg_seq[k][j]);
            fprintf(OUTPUT,"\n");
            fclose(OUTPUT);
    }
    for(k=0;k<genotype->ntfgenes;k++)
    {
            snprintf(filename,sizeof(char)*32,"TF_%i.txt",k);
            OUTPUT=fopen(filename,"a+");
            fprintf(OUTPUT,"%d",step_i);
            for(j=0;j<TF_ELEMENT_LEN;j++)
                    fprintf(OUTPUT,"%c",genotype->tf_seq[k][j]);
            fprintf(OUTPUT,"\n");
            fclose(OUTPUT);
            snprintf(filename,sizeof(char)*32,"TF_r_%i.txt",k);
            OUTPUT=fopen(filename,"a+");
            fprintf(OUTPUT,"%d",step_i);
            for(j=0;j<TF_ELEMENT_LEN;j++)
                    fprintf(OUTPUT,"%c",genotype->tf_seq_rc[k][j]);
            fprintf(OUTPUT,"\n");
            fclose(OUTPUT);
    }
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
    int table_reduced_sites[NGENES][NGENES];
        
    for(i=0;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->ngenes;j++)
        {
            table[i][j]=0; 
            table_reduced_sites[i][j]=0;
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
    OUTPUT1=fopen("summary_BS","a+");
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
    fprintf(OUTPUT1,"\n");
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
        fprintf(OUTPUT1,"\n");
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

void print_binding_sites_distribution(Genotype *genotype, int gene_id, int init_or_end)
{
    FILE *fp;
    int i;
    char filename[32];
    snprintf(filename,sizeof(char)*32,"distribution_on_%i_%i.txt",gene_id,init_or_end);
    fp=fopen(filename,"w");
    fprintf(fp,"y tf_id pos AR\n");   
    for(i=0;i<genotype->binding_sites_num[gene_id];i++)
    {
        fprintf(fp,"0 %d %d %d\n",genotype->all_binding_sites[gene_id][i].tf_id,genotype->all_binding_sites[gene_id][i].BS_pos+8,genotype->activating[genotype->all_binding_sites[gene_id][i].tf_id]);
    }
    fclose(fp);   
}

void resolve_overlapping_sites(Genotype *genotype, int gene_id, int N_non_overlapping_sites[NGENES])
{
    int i,head,tail;
    int temp[NGENES][100]; // assuming there are at most 100 BS of the same TF on a promoter.
    int temp2[NGENES]; // counter of the number of 
    
    for(i=0;i<genotype->ngenes-1;i++)
        temp2[i]=0;
    
    for(i=0;i<genotype->binding_sites_num[gene_id];i++)
    {
        temp[genotype->all_binding_sites[gene_id][i].tf_id][temp2[genotype->all_binding_sites[gene_id][i].tf_id]]=i;
        temp2[genotype->all_binding_sites[gene_id][i].tf_id]++;
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
                if(genotype->all_binding_sites[gene_id][temp[i][tail]].BS_pos-genotype->all_binding_sites[gene_id][temp[i][head]].BS_pos-1<TF_ELEMENT_LEN+2*HIND_LENGTH)
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

//void whole_gene_deletion(Genotype *genotype,long int *seed) // only tf genes get deleted. Haven't been debugged
//{
//    float random;
//    int pos_g, offset, i,j;
//    char *temp1;
//    float *mRNAdecay, *proteindecay, *translation, *pic_dis;
//    int protein_id;
//    
//    random=ran1(seed)*genotype->ntfgenes;			
//    pos_g=floor(random);				
//    temp1 = &genotype->cisreg_seq[pos_g][0];			
//    offset=CISREG_LEN;
//
//    /* remove it from cis_seq */
//    for(i=0;i<CISREG_LEN;i++) 
//    {				
//        *temp1=*(temp1+offset);         /* move sequence to fill the gap */
//        temp1++;				
//    }   
//    
//    /* remove from binding seq */
//    temp1=&genotype->tf_seq[pos_g][0][0];
//    offset=TF_ELEMENT_LEN;
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp1=*(temp1+offset);
//        temp1++;
//    }
//    
//    /* remove from rc binding seq */    
//    temp1=&genotype->tf_seq_rc[pos_g][0][0];
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp1=*(temp1+offset);
//        temp1++;
//    }
//    
//    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
//    for(i=pos_g;i<genotype->ngenes;i++)
//    {
//        genotype->pic_disassembly[i][0]=genotype->pic_disassembly[i+1][0];             /* shift elements in the array*/
//        genotype->mRNAdecay[i]=genotype->mRNAdecay[i+1];
//        genotype->proteindecay[i]=genotype->proteindecay[i+1];
//        genotype->translation[i]=genotype->translation[i+1];
//        genotype->re_calc[i][0]=genotype->re_calc[i+1][0];
//        genotype->re_calc[i][1]=genotype->re_calc[i+1][1];
//        genotype->re_calc[i][2]=genotype->re_calc[i+1][2];
//        genotype->re_calc[i][3]=1;              /* copy info back to this site in clone_cell */
//        
//        /* now move the info about binding sites*/
//        for(j=0;j<genotype->binding_sites_num[i+1];j++)
//        {
//            genotype->all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i+1][j].tf_id;
//            genotype->all_binding_sites[i][j].Koff=genotype->all_binding_sites[i+1][j].Koff;
//            genotype->all_binding_sites[i][j].BS_pos=genotype->all_binding_sites[i+1][j].BS_pos;
//            genotype->all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i+1][j].N_hindered;
//        }
//        genotype->binding_sites_num[i]=genotype->binding_sites_num[i+1];
//        genotype->max_hindered_sites[i]=genotype->max_hindered_sites[i+1];
//        genotype->N_act_BS[i]=genotype->N_act_BS[i+1];
//        genotype->N_rep_BS[i]=genotype->N_rep_BS[i+1];
//        
////        for(j=0;j<genotype->max_N_rep_bound[i+1];j++)
////        {
////            genotype->N_configurations[i][j]=genotype->N_configurations[i+1][j];
////        }        
////        genotype->max_N_rep_bound[i]=genotype->max_N_rep_bound[i+1];
////        genotype->max_N_act_bound[i]=genotype->max_N_act_bound[i+1];
//    }   
//    
//    /* change ngenes and ntfgenes*/
//    genotype->ngenes--;
//    genotype->ntfgenes--;
//    
//    /* now change protein_pool*/
//    protein_id=genotype->which_protein[pos_g];
//    
//    if(genotype->protein_pool[protein_id][0][0]==1) /* if this is the only gene copy */
//    {    
//        for(i=0;i<genotype->nproteins-protein_id;i++)  /* then we need to remove this protein from protein_pool*/
//        {
//            genotype->protein_pool[protein_id][0][0]=genotype->protein_pool[protein_id+1][0][0];
//            for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
//            {
//                genotype->protein_pool[protein_id][1][j]=genotype->protein_pool[protein_id+1][1][j];
//            }
//            protein_id++;
//        }
//        
//        if(genotype->activating[protein_id] && (genotype->N_act!=0)) /* reduce the number of activator if necessary */
//            genotype->N_act--;
//        if(!genotype->activating[protein_id] && (genotype->N_rep!=0))
//            genotype->N_rep--;
//                   
//        for(i=0;i<genotype->nproteins-protein_id-2;i++) /* also remove it from activating */
//        {
//            genotype->activating[protein_id]=genotype->activating[protein_id+1];
//        }
//        
//        genotype->nproteins--;                /* reduce the number of protein*/
//    }  
//    else
//    {
//        i=0;
//        while(genotype->protein_pool[protein_id][1][i]!=pos_g)i++; /* find where is this pos_g */
//        for(;i<genotype->protein_pool[protein_id][0][0];i++)
//        {
//            genotype->protein_pool[protein_id][1][i]=genotype->protein_pool[protein_id][1][i+1];
//        }
//        genotype->protein_pool[protein_id][0][0]--;
//    }
//}

//void gene_duplication(Genotype *genotype,long int *seed) //only tf genes can be duplicated
//{
//    float random;
//    int pos_g, i,j, protein_id;
//    char *temp1, *temp2;
//    
//    random=ran1(seed)*genotype->ntfgenes;
//    pos_g=floor(random);
//    
//    /* copy the promoter*/
//    temp1=&genotype->cisreg_seq[pos_g][0]; /* points to the gene to be duplicated*/
//    temp2=&genotype->cisreg_seq[genotype->ngenes-1][CISREG_LEN-1]; /* points to the end of the 2nd selection gene */
//    
//    for(i=0;i<2*CISREG_LEN;i++) 
//    {
//        *(temp2+CISREG_LEN)=*temp2; /* shift the sequences of the selection genes CISREG_LEN bp */
//        temp2--;
//    }
//    
//    temp2=&genotype->cisreg_seq[genotype->ngenes-2][0]; /* point temp2 to the start of the original first selection gene */
//    
//    for(i=0;i<CISREG_LEN;i++) 
//    {
//        *temp2++=*temp1++;  /* put the duplicated gene at the original place of the first selection gene */
//    }
//    
//    /* copy the binding sequence*/
//    temp1=&genotype->tf_seq[pos_g][0][0];
//    temp2=&genotype->tf_seq[genotype->ntfgenes][0][0];   
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp2++=*temp1++;
//    }
//    
//    /* copy the rc binding sequence*/
//    temp1=&genotype->tf_seq_rc[pos_g][0][0];
//    temp2=&genotype->tf_seq_rc[genotype->ntfgenes][0][0];    
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp2++=*temp1++;
//    }  
//    
//    /* add it to protein_pool, but do not change nproteins*/    
//    protein_id=genotype->which_protein[pos_g];
//    
//    genotype->protein_pool[protein_id][1][genotype->protein_pool[protein_id][0][0]]=genotype->ntfgenes;
//    
//    genotype->protein_pool[protein_id][0][0]++; 
//    
//    /* add it to PIC_assembly, mRNAdecay, proteinDecay, and translation*/
//    for(i=genotype->ngenes;i>genotype->ntfgenes;i--)
//    {
//        genotype->pic_disassembly[i][0]=genotype->pic_disassembly[i-1][0]; /* shift the selection genes to make a slot*/
//        genotype->mRNAdecay[i]=genotype->mRNAdecay[i-1];
//        genotype->proteindecay[i]=genotype->proteindecay[i-1];
//        genotype->translation[i]=genotype->translation[i-1];
//        genotype->re_calc[i][0]=genotype->re_calc[i-1][0];
//        genotype->re_calc[i][1]=genotype->re_calc[i-1][1];
//        genotype->re_calc[i][2]=genotype->re_calc[i-1][2];
//        genotype->re_calc[i][3]=1;        /* copy info back to this site in clone_cell */                   
//        
//         /* now move the info about binding sites*/
//        for(j=0;j<genotype->binding_sites_num[i-1];j++)
//        {
//            genotype->all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i-1][j].tf_id;
//            genotype->all_binding_sites[i][j].Koff=genotype->all_binding_sites[i-1][j].Koff;
//            genotype->all_binding_sites[i][j].BS_pos=genotype->all_binding_sites[i-1][j].BS_pos;
//            genotype->all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i-1][j].N_hindered;
//        }
//        genotype->binding_sites_num[i]=genotype->binding_sites_num[i-1];
//        genotype->max_hindered_sites[i]=genotype->max_hindered_sites[i-1];
//        genotype->N_act_BS[i]=genotype->N_act_BS[i-1];
//        genotype->N_rep_BS[i]=genotype->N_rep_BS[i-1];
//        
////        for(j=0;j<genotype->max_N_rep_bound[i-1];j++)
////        {
////            genotype->N_configurations[i][j]=genotype->N_configurations[i-1][j];
////        }        
////        genotype->max_N_rep_bound[i]=genotype->max_N_rep_bound[i-1];
////        genotype->max_N_act_bound[i]=genotype->max_N_act_bound[i-1];
//    }
//    
//    genotype->pic_disassembly[genotype->ntfgenes][0]=genotype->pic_disassembly[pos_g][0];
//    genotype->mRNAdecay[genotype->ntfgenes]=genotype->mRNAdecay[pos_g];
//    genotype->proteindecay[genotype->ntfgenes]=genotype->proteindecay[pos_g];
//    genotype->translation[genotype->ntfgenes]=genotype->translation[pos_g];
//    genotype->re_calc[genotype->ntfgenes][0]=pos_g;        /* we can copy the tf distribution from pos_g */
//    genotype->re_calc[genotype->ntfgenes][1]=1;            /* we also copy the tf distribution from here */
//    genotype->re_calc[genotype->ntfgenes][2]=0;            /* we do not need to calculate the binding sites */
//    genotype->re_calc[pos_g][1]=1;                         /* we just duplicated pos_g, so we can copy the tf distribution from pos_g */
//    genotype->re_calc[genotype->ntfgenes][3]=1;            /* copy info back to this site in clone_cell */
//    
//    /* now copy the info about binding sites*/
//    for(j=0;j<genotype->binding_sites_num[pos_g];j++)
//        {
//            genotype->all_binding_sites[genotype->ntfgenes][j].tf_id=genotype->all_binding_sites[pos_g][j].tf_id;
//            genotype->all_binding_sites[genotype->ntfgenes][j].Koff=genotype->all_binding_sites[pos_g][j].Koff;
//            genotype->all_binding_sites[genotype->ntfgenes][j].BS_pos=genotype->all_binding_sites[pos_g][j].BS_pos;
//            genotype->all_binding_sites[genotype->ntfgenes][j].N_hindered=genotype->all_binding_sites[pos_g][j].N_hindered;
//        }
//        genotype->binding_sites_num[genotype->ntfgenes]=genotype->binding_sites_num[pos_g];
//        genotype->max_hindered_sites[genotype->ntfgenes]=genotype->max_hindered_sites[pos_g];
//        genotype->N_act_BS[genotype->ntfgenes]=genotype->N_act_BS[pos_g];
//        genotype->N_rep_BS[genotype->ntfgenes]=genotype->N_rep_BS[pos_g];
//        
////        for(j=0;j<genotype->max_N_rep_bound[pos_g];j++)
////        {
////            genotype->N_configurations[genotype->ntfgenes][j]=genotype->N_configurations[pos_g][j];
////        }        
////        genotype->max_N_rep_bound[genotype->ntfgenes]=genotype->max_N_rep_bound[pos_g];
////        genotype->max_N_act_bound[genotype->ntfgenes]=genotype->max_N_act_bound[pos_g];
//    
//    /* update gene numbers*/    
//    genotype->ntfgenes++;
//    genotype->ngenes++;
//}

/******************************************************************************/
/***originally, clone_cell copies TFBS back to genes that are mutated**********/
/***in order to avoid recomputing TFBS*****************************************/
/******************************************************************************/
//void clone_cell(Genotype *genotype_orig,                
//                Genotype *genotype_clone,
//                int clone_type)
//{
//    int i, j;
//        
//    if(clone_type!=4) /* not a mutation in rate constant*/
//    {
//        for(i=0; i< genotype_orig->ngenes;i++)
//        {
//            if(genotype_clone->re_calc[i][3]  /* only copy to places that mutated*/
//               || genotype_orig->re_calc[i][3]) /* this argument is used when copy new genotype to current genotype, in try_fixation*/                
//            {
//                /* we copy binding sites info here, because recalc these info is more expensive*/
//                for(j=0;j<genotype_orig->binding_sites_num[i];j++)
//                {
//                    genotype_clone->all_binding_sites[i][j].tf_id=genotype_orig->all_binding_sites[i][j].tf_id;
//                    genotype_clone->all_binding_sites[i][j].Koff=genotype_orig->all_binding_sites[i][j].Koff;
//                    genotype_clone->all_binding_sites[i][j].BS_pos=genotype_orig->all_binding_sites[i][j].BS_pos;
//                    genotype_clone->all_binding_sites[i][j].N_hindered=genotype_orig->all_binding_sites[i][j].N_hindered;
//                }
//                genotype_clone->binding_sites_num[i]=genotype_orig->binding_sites_num[i];
//                genotype_clone->max_hindered_sites[i]=genotype_orig->max_hindered_sites[i];
//                genotype_clone->N_act_BS[i]=genotype_orig->N_act_BS[i];
//                genotype_clone->N_rep_BS[i]=genotype_orig->N_rep_BS[i];
//        
////                for(j=0;j<genotype_orig->max_N_rep_bound[i];j++)
////                {
////                    genotype_clone->N_configurations[i][j]=genotype_orig->N_configurations[i][j];
////                }        
////                genotype_clone->max_N_rep_bound[i]=genotype_orig->max_N_rep_bound[i];
////                genotype_clone->max_N_act_bound[i]=genotype_orig->max_N_act_bound[i];
//                
//                if(clone_type!=3) /* if the mutation was not in binding sequence */
//                {
//                    memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_orig->cisreg_seq[i][0],CISREG_LEN*sizeof(char));
//                }
//            }
//        }
//    }
//    
//    if(clone_type!=1) /* not a substitution or indel*/
//    {
//        for (i=0; i < genotype_orig->ngenes; i++) 
//        {              
//            genotype_clone->mRNAdecay[i]=genotype_orig->mRNAdecay[i];
//            genotype_clone->proteindecay[i]=genotype_orig->proteindecay[i];
//            genotype_clone->translation[i]=genotype_orig->translation[i];
//            genotype_clone->re_calc[i][0]=genotype_orig->re_calc[i][0];
//            genotype_clone->re_calc[i][1]=genotype_orig->re_calc[i][1];
//            genotype_clone->pic_disassembly[i][0]=genotype_orig->pic_disassembly[i][0];
//            genotype_clone->which_protein[i]=genotype_orig->which_protein[i];
//        }
//        
//        for(i=0;i<genotype_orig->nproteins;i++)
//        {
//            genotype_clone->activating[i]= genotype_orig->activating[i];  
//            genotype_clone->protein_pool[i][0][0]=genotype_orig->protein_pool[i][0][0];
//
//            for(j=0;j<MAXALLOC;j++)
//            {
//                genotype_clone->protein_pool[i][1][j]=genotype_orig->protein_pool[i][1][j];
//            }         
//        }
//    }
//    
//    /* since there is no tag to mark which binding seq is mutated, we copy all*/
//    for (i=0; i < genotype_orig->ntfgenes; i++) 
//    {        
//        for(j=0;j<TF_ELEMENT_LEN;j++)
//        {    
//            genotype_clone->tf_seq[i][0][j]=genotype_orig->tf_seq[i][0][j];
//            genotype_clone->tf_seq_rc[i][0][j]=genotype_orig->tf_seq_rc[i][0][j];
//        }
//    }
//
//    /* these are easy, so do it everytime*/
//    genotype_clone->fitness=genotype_orig->fitness;
//    genotype_clone->ngenes=genotype_orig->ngenes;
//    genotype_clone->ntfgenes=genotype_orig->ntfgenes;
//    genotype_clone->nproteins=genotype_orig->nproteins;
//    genotype_clone->N_act=genotype_orig->N_act;
//    genotype_clone->N_rep=genotype_orig->N_rep;
//         
//    for(i=0;i<genotype_orig->ngenes;i++)
//    {
//        genotype_clone->re_calc[i][2]=0;   /* do not recalc binding sites unless mutation changes it*/      
//    }
//    
//    if(clone_type==6)/* if it is a fixation event*/
//    {
//        for(i=0;i<genotype_orig->ngenes;i++)
//            genotype_clone->re_calc[i][3]=1;   
//    }
//    else
//    {
//        for(i=0;i<genotype_orig->ngenes;i++)
//            genotype_clone->re_calc[i][3]=0;   /* do not copy info for this gene, unless mutation changes it*/
//    }
//}
/******************************************************************************/
/******************************************************************************/

/*this use max_N_rep_bound and *N_act_bound */
//float calc_ratio_act_to_rep(AllTFBindingSites *BS_info,
//                           int ntfgenes,
//                           int max_N_hindered_BS,
//                           int N_BS,
//                           int N_act_BS,
//                           int N_rep_BS, 
//                           int activating[NGENES],
//                           int max_N_rep_bound,
//                           int max_N_act_bound,
//                           int *N_act_bound,
//                           float protein_number[NGENES])
//{
//    double ratio_matrices[max_N_hindered_BS+1][max_N_rep_bound+1][max_N_act_bound+1];   
//    double transition_matrix[max_N_rep_bound+1][max_N_act_bound+1];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;    
//    float Kon[ntfgenes];      
//    
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,k,m,n;
//    int n_col1=1;
//    int n_col2=1;
//    int n_row=1;
//    
//    /*calc Kon based on TF concentration*/
//    for(i=0;i<ntfgenes;i++)
//    {
//        Kon[i]=kon*protein_number[i];
//    }
//    
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_N_hindered_BS+1;i++)
//    {
//        for(j=0;j<max_N_rep_bound+1;j++)
//        {
//            for(k=0;k<N_act_bound[j];k++)
//            {
//                ratio_matrices[i][j][k]=0.0;
//            }
//        }
//    }
//    
//    for(j=0;j<max_N_rep_bound+1;j++)
//    {
//        for(k=0;k<N_act_bound[j];k++)
//        {
//            transition_matrix[j][k]=0.0;            
//        }
//    }
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    
//    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
//    
//    if(activating[BS_info[0].tf_id]==1) // if a activator binds to this BS
//    {
//        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col2++;
//    }
//    else
//    {
//        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
//    }    
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
//
//        product_of_freq = Kon[BS_info[m].tf_id];
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//            {
//                product_of_freq*=BS_info[n].Koff;
//            }
//        }
//
//        switch(activating[BS_info[m].tf_id])
//        {
//            case 1: // a BS of activators
//                n_col2++;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
//                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//                  for(i=0;i<n_row;i++)
//                  {
//                      transition_matrix[i][0]=0.0;
//                      
//                      n_col1=(n_col2<N_act_bound[i])? n_col2:N_act_bound[i];
//
//                      for(j=1;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
//                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[0][1]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              
//
//                for(i=0;i<n_row;i++)
//                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                    
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//
//            case 0: // a BS of repressors
//                
//                n_row++;
//                
//                n_row=(n_row < max_N_rep_bound+1)?n_row:max_N_rep_bound+1;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1);                                 
//
//                  n_col1=(n_col2<N_act_bound[0])?n_col2:N_act_bound[0];
//                  
//                  for(j=0;j<n_col1;j++)
//                  {
//                      transition_matrix[0][j]=0.0;
//                  }
//                  
//                  for(i=1;i<n_row;i++)
//                  {
//                      n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                      
//                      for(j=0;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
//                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[1][0]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
//
//                for(i=0;i<n_row;i++)
//                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                     
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//        }
//    }
//
//    sum=0.0;
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        for(j=0;j<N_act_bound[i];j++)
//        {
//            sum+=ratio_matrices[pos_next_record][i][j];
//        }
//    }   
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
//        
//        for(;j<N_act_bound[i];j++)
//        {
//            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
//        }
//    } 
//
//    return (float)(prob_act_over_rep/sum);
//     // end of the forward algorithm   
//}


/* this one does not use N_act_bound*/
//float calc_ratio_act_to_rep(AllTFBindingSites *BS_info,
//                           int ntfgenes,
//                           int max_N_hindered_BS,
//                           int N_BS,
//                           int N_act_BS,
//                           int N_rep_BS, 
//                           int activating[NGENES],
//                           int max_N_rep_bound,
//                           int max_N_act_bound, 
//                           int *N_act_bound,
//                           float protein_number[NGENES])
//{
//    double ratio_matrices[max_N_hindered_BS+1][max_N_rep_bound+1][max_N_act_bound+1];   
//    double transition_matrix[max_N_rep_bound+1][max_N_act_bound+1];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;    
//    float Kon[ntfgenes];      
//    
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,k,m,n;
//    int n_col1=1;
//    int n_col2=1;
//    int n_row=1;
//    
//    /*calc Kon based on TF concentration*/
//    for(i=0;i<ntfgenes;i++)
//    {
//        Kon[i]=kon*protein_number[i];
//    }
//    
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_N_hindered_BS+1;i++)
//    {
//        for(j=0;j<max_N_rep_bound+1;j++)
//        {
//            for(k=0;k<max_N_act_bound+1;k++)
//            {
//                ratio_matrices[i][j][k]=0.0;
//            }
//        }
//    }
//    
//    for(j=0;j<max_N_rep_bound+1;j++)
//    {
//        for(k=0;k<max_N_act_bound+1;k++)
//        {
//            transition_matrix[j][k]=0.0;            
//        }
//    }
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    
//    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
//    
//    if(activating[BS_info[0].tf_id]==1) // if a activator binds to this BS
//    {
//        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col1++;
//    }
//    else
//    {
//        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
//    }    
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
//
//        product_of_freq = Kon[BS_info[m].tf_id];
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//            {
//                product_of_freq*=BS_info[n].Koff;
//            }
//        }
//
//        switch(activating[BS_info[m].tf_id])
//        {
//            case 1: // a BS of activators
//                n_col1++;
//                
//                n_col1=(n_col1<max_N_act_bound+1)? n_col1:max_N_act_bound+1;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
//                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//                  for(i=0;i<n_row;i++)
//                  {
//                      transition_matrix[i][0]=0.0;
//                      
////                      n_col1=(n_col2<N_act_bound[i])? n_col2:N_act_bound[i];
//
//                      for(j=1;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
////                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[0][1]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              
//
//                for(i=0;i<n_row;i++)
//                {
////                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                    
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//
//            case 0: // a BS of repressors
//                
//                n_row++;
//                
//                n_row=(n_row < max_N_rep_bound+1)?n_row:max_N_rep_bound+1;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1);                                 
//
////                  n_col1=(n_col2<N_act_bound[0])?n_col2:N_act_bound[0];
//                  
//                  for(j=0;j<n_col1;j++)
//                  {
//                      transition_matrix[0][j]=0.0;
//                  }
//                  
//                  for(i=1;i<n_row;i++)
//                  {
////                      n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                      
//                      for(j=0;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
////                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[1][0]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
//
//                for(i=0;i<n_row;i++)
//                {
////                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                     
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//        }
//    }
//
//    sum=0.0;
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        for(j=0;j<max_N_act_bound+1;j++)
//        {
//            sum+=ratio_matrices[pos_next_record][i][j];
//        }
//    }   
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
//        
//        for(;j<max_N_act_bound+1;j++)
//        {
//            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
//        }
//    } 
//
//    return (float)(prob_act_over_rep/sum);
//     // end of the forward algorithm   
//}


/* this one use neither max_N_rep_bound or *N_act_bound*/
//float calc_ratio_act_to_rep(AllTFBindingSites *BS_info,
//                           int ntfgenes,
//                           int max_N_hindered_BS,
//                           int N_BS,
//                           int N_act_BS,
//                           int N_rep_BS, 
//                           int activating[NGENES],
//                           int max_N_rep_bound,
//                           int *N_act_bound,
//                           float protein_number[NGENES])
//{
//    double ratio_matrices[max_N_hindered_BS+1][N_rep_BS+1][N_act_BS+1];   
//    double transition_matrix[N_rep_BS+1][N_act_BS+1];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;    
//    float Kon[ntfgenes];      
//    
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,k,m,n;
//    int n_col=1;    
//    int n_row=1;
//    
//    /*calc Kon based on TF concentration*/
//    for(i=0;i<ntfgenes;i++)
//    {
//        Kon[i]=kon*protein_number[i];
//    }
//    
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_N_hindered_BS+1;i++)
//    {
//        for(j=0;j<N_rep_BS+1;j++)
//        {
//            for(k=0;k<N_act_BS+1;k++)
//            {
//                ratio_matrices[i][j][k]=0.0;
//            }
//        }
//    }
//    
//    for(j=0;j<max_N_rep_bound+1;j++)
//    {
//        for(k=0;k<N_act_BS+1;k++)
//        {
//            transition_matrix[j][k]=0.0;            
//        }
//    }
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    
//    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
//    
//    if(activating[BS_info[0].tf_id]==1) // if a activator binds to this BS
//    {
//        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col++;
//    }
//    else
//    {
//        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
//    }    
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
//
//        product_of_freq = Kon[BS_info[m].tf_id];
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//            {
//                product_of_freq*=BS_info[n].Koff;
//            }
//        }
//
//        switch(activating[BS_info[m].tf_id])
//        {
//            case 1: // a BS of activators
//                n_col++;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
//                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//                  for(i=0;i<n_row;i++)
//                  {
//                      transition_matrix[i][0]=0.0;
//
//                      for(j=1;j<n_col;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
//                        for(j=0;j<n_col;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[0][1]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              
//
//                for(i=0;i<n_row;i++)
//                {
//                    for(j=0;j<n_col;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//
//            case 0: // a BS of repressors
//                
//                n_row++;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); 
//                  
//                  for(j=0;j<n_col;j++)
//                  {
//                      transition_matrix[0][j]=0.0;
//                  }
//                  
//                  for(i=1;i<n_row;i++)
//                  {                      
//                      for(j=0;j<n_col;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {                        
//                        for(j=0;j<n_col;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[1][0]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
//
//                for(i=0;i<n_row;i++)
//                {                     
//                    for(j=0;j<n_col;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//        }
//    }
//
//    sum=0.0;
//
//    for(i=0;i<N_rep_BS+1;i++)
//    {
//        for(j=0;j<n_col;j++)
//        {
//            sum+=ratio_matrices[pos_next_record][i][j];
//        }
//    }   
//
//    for(i=0;i<N_rep_BS+1;i++)
//    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
//        
//        for(;j<n_col;j++)
//        {
//            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
//        }
//    } 
//
//    return (float)(prob_act_over_rep/sum);
//     // end of the forward algorithm   
//}

/* calc the numbers of TF binding configurations:
 * determine the maximal number of repressors that can bind to a promoter
 * and the numbers of activators that can bind given the number of repressors that bind
 * Use these information to accelerate the calculation of rep-to-act ratios.
 */
//void calc_configurations(Genotype *genotype, int gene_id)
//{
//    int j,k,l;
//    int hindered_rep;
//    int hindered_act; 
//    int N_rep_BS=0;
//    int N_act_BS=0;
//    
//    genotype->N_configurations[gene_id][0]=1;
//    
//    genotype->max_N_rep_bound[gene_id]=0;   
//    
//    genotype->max_N_act_bound[gene_id]=0;
//
//    for(j=0;j<genotype->binding_sites_num[gene_id];j++)
//    {
//        hindered_rep=0;
//        hindered_act=0;
//
//        for(k=0;k<j;k++)
//        {
//            /* for BS within hindrance range*/
//           if((genotype->all_binding_sites[gene_id][k].BS_pos) > (genotype->all_binding_sites[gene_id][j].BS_pos-TF_ELEMENT_LEN-2*HIND_LENGTH))
//           {
//               if(genotype->activating[genotype->all_binding_sites[gene_id][k].tf_id]==0) hindered_rep++;
//               else hindered_act++;
//           }
//        }
//
//        if(genotype->activating[genotype->all_binding_sites[gene_id][j].tf_id]==0)
//        {
//            N_rep_BS++;
//            
//            if((N_rep_BS-hindered_rep)!=1)
//            {
//                genotype->max_N_rep_bound[gene_id]++; 
//                
//                for(l=1;l<genotype->max_N_rep_bound[gene_id]+1;l++)
//                {
//                    genotype->N_configurations[gene_id][l]=genotype->N_configurations[gene_id][l-1];
//                }
//            }
//        }
//        else
//        {
//            N_act_BS++;
//            
//            if((N_act_BS-hindered_act)!=1)
//            {
//                for(l=0;l<genotype->max_N_rep_bound[gene_id]-hindered_rep+1;l++)
//                    genotype->N_configurations[gene_id][l]++;
//            }
//        }
//    } 
//    
//    for(j=0;j<genotype->max_N_rep_bound[gene_id];j++)
//    {
//        genotype->max_N_act_bound[gene_id]=(genotype->max_N_act_bound[gene_id]>genotype->N_configurations[gene_id][j])?
//                                            genotype->max_N_act_bound[gene_id]:genotype->N_configurations[gene_id][j];
//    }
//    genotype->max_N_act_bound[gene_id]--;
//}

//void try_fixation(Genotype *current_genotype, Genotype *new_genotype, float ori_fitness, int *pfixation, int *OuterWhile, long int *seed)
//{	
//    if(*pfixation) return; // if other threads have reported fixation
//    else
//    {
//        float s, P_fix, ref;
//        
//        s = (new_genotype->fitness-ori_fitness)/ori_fitness;
//
//        if (fabs(s)<EPSILON){P_fix = 1/(float)POP_SIZE;}	
//        else{ P_fix = (1-exp(-s))/(1-exp(-s*(float)POP_SIZE)); }
//
//        ref=ran1(seed);
//
//        if(ref > P_fix) return;
//        
//        else 
//        {
//            int clone_type=5;            
//            clone_cell(new_genotype, current_genotype, clone_type);
//            *pfixation=1;
//            *OuterWhile+=1;
//            return;
//        }
//    }
//}

/*copy genotype from a fixed offspring back to ancestor*/
//void clone_cell_backward(Genotype *genotype_orig,                
//                        Genotype *genotype_clone,
//                        int clone_type)
//{
//    int i, j;
//          
//    if(clone_type!=MUT_KCONST) /* not a mutation in rate constant*/
//    {
//        for(i=0; i< genotype_orig->ngenes;i++)
//        {
//            genotype_clone->which_cluster[i]=genotype_orig->which_cluster[i];
//         
//            genotype_clone->recalc_TFBS[i]=0;
//            
//            if(genotype_orig->clone_info[i]==1)                
//            {               
//                if(clone_type!=MUT_TFSEQ) /* if the mutation was not in binding sequence */
//                {
//                    memcpy(&genotype_clone->cisreg_seq[i][0][0],&genotype_orig->cisreg_seq[i][0][0],CISREG_LEN*sizeof(char));
//                }
//                genotype_clone->recalc_TFBS[i]=1;
//            }
//            
//            genotype_clone->clone_info[i]=1; 
//        }
//        
//        /*reset clone's cisreg_cluster*/
//        i=0;
//        while(genotype_clone->cisreg_cluster[i][0]!=-1)
//        {
//            j=0;
//            while(genotype_clone->cisreg_cluster[i][j]!=-1)
//            {
//                genotype_clone->cisreg_cluster[i][j]=-1;
//                j++;
//            }
//            i++;
//        }        
//        /*then copy from orig*/
//        i=0;
//        while(genotype_orig->cisreg_cluster[i][0]!=-1)
//        {
//            j=0;
//            while(genotype_orig->cisreg_cluster[i][j]!=-1)
//            {
//                genotype_clone->cisreg_cluster[i][j]=genotype_orig->cisreg_cluster[i][j];
//                j++;
//            }
//            i++;
//        }
//    }
//    
//    if(clone_type!= MUT_N_GENE) /* no change in gene numbers*/
//    {
//        for (i=0; i < genotype_orig->ngenes; i++) 
//        {              
//            genotype_clone->mRNAdecay[i]=genotype_orig->mRNAdecay[i];
//            genotype_clone->proteindecay[i]=genotype_orig->proteindecay[i];
//            genotype_clone->translation[i]=genotype_orig->translation[i];           
//            genotype_clone->pic_disassembly[i][0]=genotype_orig->pic_disassembly[i][0];
//            genotype_clone->which_protein[i]=genotype_orig->which_protein[i];
//        }
//        
//        for(i=0;i<genotype_orig->nproteins;i++)
//        {
//            genotype_clone->activating[i]= genotype_orig->activating[i];  
//            genotype_clone->protein_pool[i][0][0]=genotype_orig->protein_pool[i][0][0];
//
//            for(j=0;j<MAXALLOC;j++)
//            {
//                genotype_clone->protein_pool[i][1][j]=genotype_orig->protein_pool[i][1][j];
//            }         
//        }
//    }
//    
//    /* since there is no tag to mark which binding seq is mutated, we copy all*/
//    for (i=0; i < genotype_orig->ntfgenes; i++) 
//    {        
//        for(j=0;j<TF_ELEMENT_LEN;j++)
//        {    
//            genotype_clone->tf_seq[i][0][j]=genotype_orig->tf_seq[i][0][j];
//            genotype_clone->tf_seq_rc[i][0][j]=genotype_orig->tf_seq_rc[i][0][j];
//        }
//    }
//
//    /* these are easy, so do it everytime*/
//    genotype_clone->fitness=genotype_orig->fitness;
//    genotype_clone->ngenes=genotype_orig->ngenes;
//    genotype_clone->ntfgenes=genotype_orig->ntfgenes;
//    genotype_clone->nproteins=genotype_orig->nproteins;
//    genotype_clone->N_act=genotype_orig->N_act;
//    genotype_clone->N_rep=genotype_orig->N_rep;     
//}
