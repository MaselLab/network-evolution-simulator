/* 
 * This file contains functions to initialize simulation with specified 
 * selection condition, and to output summary of genotypes and network structure.
 * 
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "netsim.h"
#include "cellular_activity.h"
#include "mutation.h"
#include "numerical.h"
#include "lib.h"
#include "RngStream.h"

#define INITIALIZATION -1

int MAX_TFBS_NUMBER=100; 

static const float PROB_ACTIVATING=0.62;
static const float MEAN_PROTEIN_DECAY_RATE=-1.88;
static const float SD_PROTEIN_DECAY_RATE=0.56;
static const float MEAN_ACT_TO_INT_RATE=1.27;
static const float SD_ACT_TO_INT_RATE=0.226;
static const float MEAN_MRNA_DECAY_RATE=-1.49;
static const float SD_MRNA_DECAY_RATE=0.267;
static const float SD_PROTEIN_SYN_RATE=0.416;
static const float MEAN_PROTEIN_SYN_RATE=0.322;
static const float SD_GENE_LENGTH=0.34;
static const float MIN_Kd=1.0e-9;
static const float MAX_Kd=1.0e-6;
static const float log_MIN_Kd=-9.0;
static const float log_MAX_Kd=-6.0;
static const float NS_Kd=1.0e-5;
const float KD2APP_KD=1.8e10;
const float MEAN_GENE_LENGTH=2.568; //log10(aa)

static const float MIN_SELECTION_COEFFICIENT=1.0e-8;
/*Bounds*/
const float MAX_ACT_TO_INT_RATE=64.6;
const float MIN_ACT_TO_INT_RATE=0.59;
const float MAX_MRNA_DECAY=0.54;
const float MIN_MRNA_DECAY=7.5e-4;
const float MAX_PROTEIN_DECAY=0.69;
const float MIN_PROTEIN_DECAY=4.5e-6;
const float MAX_PROTEIN_SYN_RATE=61.4;
const float MIN_PROTEIN_SYN_RATE=4.5e-3;
const float MAX_KD=1.0e-5;
const float MIN_KD=0.0;
const int MAX_GENE_LENGTH=5000; //aa
const int MIN_GENE_LENGTH= 50; //aa

/*fitness*/
const float sampling_interval=1.0;
static const float exp_cost_factor=1.0;

/******************************************************************************
 * 
 *                     Private function prototypes
 *
 *****************************************************************************/
static void initialize_sequence(char *, int, int, RngStream);

static void initialize_genotype_fixed(Genotype *, int, int, int, int, RngStream);

static void run_fitness_test(Genotype *, CellState *, GillespieRates *, Environment *, Phenotype *, float *, int *, float *, int, RngStream, int);

static void calc_avg_fitness(Genotype *, Selection *, Phenotype [N_REPLICATES], int [MAX_GENES], float [MAX_PROTEINS], RngStream [N_THREADS], float *, float *, int); 

static void copy_environment(Environment *, Environment *);

static void clone_genotype(Genotype *, Genotype *);

static void try_replacement(Genotype *, Genotype *, int *, float*);

static void summarize_binding_sites(Genotype *,int);

static void set_signal(CellState *, Environment *, float, RngStream, int);

static int evolve_N_steps(Genotype *, Genotype *,  Mutation *, Selection *, Output_buffer [OUTPUT_INTERVAL], int *, int *, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS], int);

static void run_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_PROTEINS], int, int, RngStream, RngStream [N_THREADS]);

static void continue_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS]);

static void calc_fitness_stats(Genotype *, Selection *, float (*)[N_REPLICATES], float (*)[N_REPLICATES], int);

static float calc_replicate_fitness(CellState *, Environment *);

static void replay_mutations(Genotype *, Mutation *, int);

static void classify_all_mutations(Genotype *, Mutation *, Selection *, int, int [MAX_GENES], float [MAX_GENES], RngStream [N_THREADS]);

static void modify_network(Genotype *, Genotype *, Mutation *, Selection *, int [MAX_GENES], float [MAX_GENES], RngStream [N_THREADS]);

static void find_motifs(Genotype *);

static int find_TFBS_of_A_on_B(Genotype *, int, int, int);

static void tidy_output_files(char*, char*);

static void print_motifs(Genotype *);

static void store_mutant_info(Genotype *, Mutation *, Output_buffer *, int, int);

static void store_resident_info(Genotype *, Mutation *, Output_buffer *, int , int, int, float, int);

static void output_mutant_info(Output_buffer *, int);

static void output_resident_info(Output_buffer [OUTPUT_INTERVAL], int, int);

static void sample_parameters(Genotype *, int, RngStream);

static void mark_genes_for_sampling(Genotype *, int, int, int, char);

static void who_regulates_effector(Genotype *, int , int [MAX_PROTEINS], int [MAX_PROTEINS], int *, int *);

static int determine_adaptive_TFBSs(Genotype *, Selection *, float, float, int [MAX_GENES], float [MAX_GENES], RngStream [N_THREADS], int, float [10], float);

static void sample_expression_lvl(Genotype *, Selection *, int [MAX_GENES], float [MAX_GENES], RngStream [N_THREADS]);

static void sample_effector_expression_lvl(Genotype *, Selection *, Mutation *, int [MAX_GENES], float [MAX_GENES], int, RngStream[N_THREADS]);
  
static void restore_all_TFBS(Genotype *);

static void ignore_TFBS(Genotype*, int, int, int, int*);

/***************************************************************************** 
 * 
 *                              Global functions
 *
 *****************************************************************************/
int evolve_under_selection(Genotype *resident, 
                            Genotype *mutant, 
                            Mutation *mut_record, 
                            Selection *burn_in, 
                            Selection *selection, 
                            int init_mRNA[MAX_GENES], 
                            float init_protein[MAX_GENES],
                            RngStream RS_main,
                            RngStream RS_parallel[N_THREADS])
{  
    int i; 
    int init_step;
    Phenotype phenotype[N_REPLICATES];
    FILE *fp; 

    /*create threads and rng streams*/
    omp_set_num_threads(N_THREADS);  
    
    /* continue a simulation from a previously saved state?*/    
    fp=fopen("saving_point.txt","r");   
    if(fp!=NULL) 
    {        
        int replay_N_steps=0, int_buffer;       
        fscanf(fp,"%d %d",&replay_N_steps,&int_buffer);
        fclose(fp);
        fp=fopen(setup_summary,"a+");
        fprintf(fp,"Continue simulation from step %d\n",replay_N_steps);
        fclose(fp);
        if(replay_N_steps!=0)                              
            continue_simulation(resident, 
                                mutant, 
                                mut_record,
                                burn_in,
                                selection,
                                replay_N_steps, 
                                init_mRNA,
                                init_protein,                                 
                                RS_main,
                                RS_parallel);           
    }
    else /* otherwise the simulation starts over from beginning*/
    {   
        /* record the initial network topology*/
        init_step=0;
        summarize_binding_sites(resident,init_step); /*snapshot of the initial (0) distribution binding sites */   
        find_motifs(resident); 
        print_motifs(resident);           

        /*calculate the fitness of the initial genotype*/      
        float GR1[HI_RESOLUTION_RECALC][N_REPLICATES],GR2[HI_RESOLUTION_RECALC][N_REPLICATES]; 
        if(burn_in->MAX_STEPS!=0)
        {
            for(i=0;i<HI_RESOLUTION_RECALC;i++)  
                calc_avg_fitness(resident, burn_in, phenotype, init_mRNA, init_protein, RS_parallel, GR1[i], GR2[i], RM_NONE);
            calc_fitness_stats(resident,burn_in,&(GR1[0]),&(GR2[0]),HI_RESOLUTION_RECALC); 
        }
        else
        {
            for(i=0;i<HI_RESOLUTION_RECALC;i++)  
                calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, GR1[i], GR2[i], RM_NONE);  
            calc_fitness_stats(resident,selection,&(GR1[0]),&(GR2[0]),HI_RESOLUTION_RECALC); 
        }
       
        /* make title of the output file*/
        fp=fopen(evo_summary,"w");
        fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step N_hit_bound accepted_mut selection_coeff avg_fitness fitness1 fitness2 se_avg_fitness se_fitness1 se_fitness2 N_genes N_effector_genes N_proteins\n");
        fprintf(fp,"0 0 0 0 na na %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d\n",  
                resident->avg_fitness,               
                resident->fitness1,
                resident->fitness2,           
                resident->SE_avg_fitness,
                resident->SE_fitness1,
                resident->SE_fitness2,        
                resident->ngenes,
                resident->n_output_genes,
                resident->N_node_families);
        fclose(fp);
        run_simulation( resident, 
                        mutant, 
                        mut_record, 
                        burn_in,
                        selection,                        
                        init_mRNA,
                        init_protein,
                        0, // this is the number of total mutations that have been tried
                        1, // this tells the program from which step the simulation begins
                        RS_main,
                        RS_parallel);     

    }
    print_mutatable_parameters(resident,1);
    
    /*delete rng streams*/
    for(i=0;i<N_THREADS;i++)
        RngStream_DeleteStream (&RS_parallel[i]);   
    
    return 1;	
}

#if NEUTRAL
void evolve_neutrally(Genotype *resident, Mutation *mut_record, Selection *burn_in, Selection *selection, RngStream RS_main)
{
    int i, output_counter=0;
    FILE *fp;    
    Output_buffer resident_info[OUTPUT_INTERVAL];
    
    /*Set resident fitness to 0*/
    resident->fitness1=0.0;
    resident->fitness2=0.0;
    resident->avg_fitness=0.0;
    resident->SE_avg_fitness=0.0;
    resident->SE_fitness1=0.0;
    resident->SE_fitness2=0.0;
    
    /*Create title for output files*/
    fp=fopen(evo_summary,"a+");
    fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step N_hit_bound accepted_mut selection_coeff avg_fitness fitness1 fitness2 se_avg_fitness se_fitness1 se_fitness2 N_genes N_proteins N_activator N_repressor\n");
    fprintf(fp,"0 0 0 na na 0.0 0.0 0.0 0.0 0.0 0.0 0 0 0 0 \n");
    fclose(fp); 
    
    /*record the initial network*/
    calc_all_binding_sites(resident,NMIN);
    summarize_binding_sites(resident,0);
    
    /*set BURN-IN conditions*/              
    DUPLICATION=burn_in->temporary_DUPLICATION;                 
    SILENCING=burn_in->temporary_SILENCING;
//    N_EFFECTOR_GENES=burn_in->temporary_N_effector_genes;
//    N_TF_GENES=burn_in->temporary_N_tf_genes; 
    miu_ACT_TO_INT_RATE=burn_in->temporary_miu_ACT_TO_INT_RATE; 
    miu_Kd=burn_in->temporary_miu_Kd;       
    miu_protein_syn_rate=burn_in->temporary_miu_protein_syn_rate; 
    
    /*burn in*/
    for(i=1;i<=burn_in->MAX_STEPS;i++)
    {  
        mutate(resident,RS_main,mut_record);  
        calc_all_binding_sites(resident,NMIN);
        find_motifs(resident);       
        store_resident_info(resident, mut_record, &(resident_info[output_counter]), i, 1, i, 0.0, 1); //magic number 1 means store everything
        output_counter++;
        /*output network topology every OUTPUT_INTERVAL steps*/ 
        if(i%OUTPUT_INTERVAL==0)
        {
            summarize_binding_sites(resident,i);
            output_resident_info(resident_info,output_counter,1); //magic number 1 means to output everything
            output_counter=0;
        }       
    } 
    
    /*set post burn-in condition*/
    DUPLICATION=selection->temporary_DUPLICATION;                 
    SILENCING=selection->temporary_SILENCING;
//    N_EFFECTOR_GENES=selection->temporary_N_effector_genes;
//    N_TF_GENES=selection->temporary_N_tf_genes; 
    miu_ACT_TO_INT_RATE=selection->temporary_miu_ACT_TO_INT_RATE; 
    miu_Kd=selection->temporary_miu_Kd;       
    miu_protein_syn_rate=selection->temporary_miu_protein_syn_rate;      
    
    for(;i<=selection->MAX_STEPS;i++)
    {       
        mutate(resident,RS_main,mut_record);   
        calc_all_binding_sites(resident,NMIN);
        find_motifs(resident);            
        store_resident_info(resident, mut_record, &(resident_info[output_counter]), i, 1, i, 0.0, 1); //magic number 1 means store everything
        output_counter++;
        /*output network topology every OUTPUT_INTERVAL steps*/        
        if(i%OUTPUT_INTERVAL==0 && i!=burn_in->MAX_STEPS)
        {
            summarize_binding_sites(resident,i);
            output_resident_info(resident_info,output_counter,1); //magic number 1 means to output everything
            output_counter=0;
        } 
    }
    print_mutatable_parameters(resident,1);    
}
#endif

#if PHENOTYPE
void show_phenotype(Genotype *resident, Mutation *mut_record, Selection *selection, int init_mRNA[MAX_GENES], float init_protein[MAX_GENES], int replay_N_steps, RngStream RS_parallel[N_THREADS])
{      
    omp_set_num_threads(N_THREADS); 
    /*show expression timecourse of all genes at a specific evolutionary step*/        
    if(SAMPLE_GENE_EXPRESSION)
    {
        replay_mutations(resident, mut_record, replay_N_steps); 
        remove("N_motifs.txt");
        remove("networks.txt");
        sample_expression_lvl(resident, selection, init_mRNA, init_protein, RS_parallel);
    }
    
    /*show expression timecourse of effector protein during the first N evolutionary steps*/   
    if(SAMPLE_EFFECTOR_EXPRESSION)
        sample_effector_expression_lvl(resident, selection, mut_record, init_mRNA, init_protein, replay_N_steps, RS_parallel);
}
#endif

#if PERTURB
void characterize_network(Genotype *resident,
                            Genotype *mutant, 
                            Mutation *mut_record,  
                            Selection *selection,
                            int max_evolutionary_steps,
                            int init_mRNA[MAX_GENES],
                            float init_protein[MAX_PROTEINS],
                            RngStream RS_parallel[N_THREADS])
{
    if(CLEAN_UP_NETWORK)
        modify_network(resident, mutant, mut_record, selection, init_mRNA, init_protein, RS_parallel); 
    
    if(CLASSIFY_MUTATION)
        classify_all_mutations(resident, mut_record, selection, max_evolutionary_steps, init_mRNA, init_protein, RS_parallel);
}
#endif

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

/*
 *Set default values and allocate space for variables in Genotype
 */
void initialize_cache(Genotype *genotype)
{
    int j,k; 
    /*Initialize variables that applies to loci*/
    for(j=0;j<MAX_GENES;j++)
    {  
        genotype->which_protein[j]=NA;                
        genotype->recalc_TFBS[j]=YES;
        genotype->which_cluster[j]=NA; 
        genotype->min_N_activator_to_transc[j]=MAX_BINDING+1; /*by default a gene cannot be turned on. 
                                                       *MAX_BINDING is the maximum number of tf that 
                                                       *can bind to a cis-reg sequence.*/        
        genotype->which_node_family[j]=NA;
        genotype->is_output[j]=NON_OUTPUT_PROTEIN;
        genotype->cisreg_cluster_pool[j][0][0]=0;
        for(k=0;k<MAX_GENES;k++)
            genotype->cisreg_cluster_pool[j][1][k]=NA;       
        for(k=0;k<MAX_PROTEINS;k++)
            genotype->locus_specific_TF_behavior[j][k]=NON_TF;
    } 
    /* initialize variables that applies to protein */
    for(j=0;j<MAX_PROTEINS;j++)
    {
        genotype->Kd[j]=-1.0;
        genotype->protein_pool[j][0][0]=0;
        genotype->node_family_pool[j][0][0]=0;
        genotype->protein_identity[j]=NA;
        for(k=0;k<MAX_GENES;k++)        
        {
            genotype->protein_pool[j][1][k]=NA;       
            genotype->node_family_pool[j][1][k]=NA;
        }
    }    
    for(j=0;j<MAX_OUTPUT_GENES;j++)        
        genotype->output_gene_ids[j]=NA;
    /* alloc space for binding sites*/
    genotype->N_allocated_elements=MAX_TFBS_NUMBER;
    for(j=0;j<MAX_GENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        if (!(genotype->all_binding_sites[j])) 
        {
#if MAKE_LOG
            LOG("Failed to allocate space\n");
#endif           
            exit(-1);
        }
    }
    /*Initialize binding sites summary*/
    for(j=N_SIGNAL_TF;j<MAX_GENES;j++)
    {
        genotype->N_act_BS[j]=0;
        genotype->N_rep_BS[j]=0;
        genotype->binding_sites_num[j]=0;
    }    
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype, 
                        int init_N_output_act, 
                        int init_N_output_rep,
                        int init_N_non_output_act,
                        int init_N_non_output_rep,
                        RngStream RS)
{ 
    int i,k;

    genotype->ngenes=N_SIGNAL_TF+init_N_output_act+init_N_output_rep+init_N_non_output_act+init_N_non_output_rep; /*including the signal genes and 1 selection gene*/  
    genotype->nproteins=genotype->ngenes;  /*at initialization, each protein is encoded by one copy of gene*/   
    genotype->n_output_genes=init_N_output_act+init_N_output_rep;
    genotype->N_node_families=genotype->nproteins;
    genotype->N_cisreg_clusters=genotype->ngenes;    
    genotype->flag_rm_which_TFBS=RM_NONE;
    /*at initialization, each copy of gene should have a unique cis-regulatory sequence*/
    for(i=0;i<genotype->ngenes;i++)
    {    
        genotype->which_cluster[i]=i; 
        genotype->cisreg_cluster_pool[i][0][0]=1;
        genotype->cisreg_cluster_pool[i][1][0]=i;
    }     
    /* initially, each protein has only one copy of gene*/    
    for(i=0;i<genotype->nproteins;i++)
    {
        genotype->protein_pool[i][0][0]=1;
        genotype->protein_pool[i][1][0]=i;
        genotype->which_protein[i]=i;
    } 
    for(i=0;i<genotype->N_node_families;i++)
    {
        genotype->which_node_family[i]=i;
        genotype->node_family_pool[i][0][0]=1;
        genotype->node_family_pool[i][1][0]=i;
    }
    /*initialize cis-reg sequences and tf binding sequences*/
    initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*MAX_GENES, genotype->ngenes, RS);  // initialize cis-reg sequence
    initialize_sequence((char *)genotype->tf_binding_seq, TF_ELEMENT_LEN*MAX_GENES, genotype->ngenes, RS);    //initialize binding sequence of TFs    
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
    initialize_genotype_fixed(genotype, init_N_output_act, init_N_output_rep, init_N_non_output_act, init_N_non_output_rep,RS);
    calc_all_binding_sites(genotype);
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
    int N_min_match;
    FILE *fperror;
    genotype->N_act_BS[gene_id]=0;
    genotype->N_rep_BS[gene_id]=0;
    genotype->max_hindered_sites[gene_id]=0;  
    //some helper pointer 
    char *tf_binding_seq;
    char *cis_seq;
    char *tf_binding_seq_rc; 
    cis_seq=&(genotype->cisreg_seq[gene_id][0]); 
  
    for(i=3;i<CISREG_LEN-TF_ELEMENT_LEN-3;i++) /* scan promoter */
    {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(genotype->all_binding_sites[gene_id][j].BS_pos>i-TF_ELEMENT_LEN-2*HIND_LENGTH)
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
        for(k=start_TF;k<genotype->nproteins;k++) 
        {            
            N_min_match=NMIN; //default value of the number of minimal matches
            if(genotype->flag_effector_is_TF==0) //if effector is not a TF
            {
                while(genotype->is_output[genotype->protein_pool[k][1][0]]==OUTPUT_PROTEIN) //and if k is an effector
                {
                    k++;                                   //then we check next protein  
                    if(k==genotype->nproteins)             //unless k is the last protein
                        break;                             //in which case we quit this for-loop                 
                }   
                if(k==genotype->nproteins)
                    break;                
#if PERTURB      
                ignore_TFBS(genotype, gene_id, k, start_TF, &N_min_match);                
#endif
            }            
            else if(genotype->flag_signal_ctrl_repressor==0) //repressors cannot be bound by signal
            {
                if(k==start_TF && genotype->protein_identity[genotype->which_protein[gene_id]]==REPRESSOR)
                    k++;
#if PERTURB
                ignore_TFBS(genotype, gene_id, k, start_TF, &N_min_match); 
#endif      
            }           
#if PERTURB
            else
            {
                ignore_TFBS(genotype, gene_id, k, start_TF, &N_min_match);
            }
#endif            
            tf_binding_seq=&(genotype->tf_binding_seq[k][0]);
            tf_binding_seq_rc=&(genotype->tf_binding_seq_rc[k][0]);            
            /*find BS on the template strand*/
            match=0;
            for(j=i;j<i+TF_ELEMENT_LEN;j++) /*calculate the number of nucleotides that match in each [i,i+TF_ELEMENT_LEN] window. The window slides by 1 each time when scanning the promoter*/
                if(cis_seq[j]==tf_binding_seq[j-i]) match++; 
                        
            if(match>=N_min_match)
            {
                if(N_binding_sites+1>=genotype->N_allocated_elements) 
                {  
                    while(genotype->N_allocated_elements<=N_binding_sites+1)
                        genotype->N_allocated_elements+=100;
                   
                    for(j=0;j<MAX_GENES;j++)
                    {
                        genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                        if (!genotype->all_binding_sites[j]) 
                        {
#if MAKE_LOG
                            LOG("error in calc_all_binding_sites_copy\n");
#endif
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
                for(j=i;j<i+TF_ELEMENT_LEN;j++)                
                    if(cis_seq[j]==tf_binding_seq_rc[j-i]) match_rc++;

                if(match_rc>=N_min_match)
                {
                    /**********************************************************************/     
                    if(N_binding_sites+1>=genotype->N_allocated_elements) 
                    {  
                        while(genotype->N_allocated_elements<=N_binding_sites+1)
                            genotype->N_allocated_elements+=100;

                        for(j=0;j<MAX_GENES;j++)
                        {
                            genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                            if (!genotype->all_binding_sites[j]) 
                            {
#if MAKE_LOG
                                LOG("error in calc_all_binding_sites_copy\n");
#endif
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
                    if(genotype->locus_specific_TF_behavior[gene_id][k]==ACTIVATOR) 
                        genotype->N_act_BS[gene_id]++;
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
    
    /* calculate max_unhindered_sites */
    /* max_unhindered_sites is maximum number of TFs that can bind to a cis-reg sequence at the same time*/
    /* We use it to faciliate the calculation of Pact and Prep. See calc_TF_dist_from_all_BS for its usage.*/
    int act_BS[MAX_TFBS_NUMBER][2],rep_BS[MAX_TFBS_NUMBER][2];
    int N_act_BS,N_rep_BS;    
    N_act_BS=1;
    N_rep_BS=1;
    for(i=0;i<N_binding_sites;i++) /* make lists BS by their types*/    
    {
        if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id]==ACTIVATOR)
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
    if(genotype->N_allocated_elements<MAX_TFBS_NUMBER)
    {
        for(gene_id=0;gene_id<MAX_GENES;gene_id++)
            genotype->all_binding_sites[gene_id]=realloc(genotype->all_binding_sites[gene_id], MAX_TFBS_NUMBER*sizeof(AllTFBindingSites));
        genotype->N_allocated_elements=MAX_TFBS_NUMBER;
    }
    
    for(gene_id=N_SIGNAL_TF;gene_id < genotype->ngenes;gene_id++)
    {        
        if(genotype->recalc_TFBS[gene_id]) /* do not calculate the binding sites if there's no mutation in the promoter or in TF binding seq*/
        {            
            calc_all_binding_sites_copy(genotype, gene_id);          
            genotype->recalc_TFBS[gene_id]=0;
        }
    }
}

/*****************************************************************************
 * 
 *                           Private functions
 *
 ****************************************************************************/

static void initialize_sequence(char *Seq, 
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
static void initialize_genotype_fixed(Genotype *genotype, 
                                    int init_N_output_act, 
                                    int init_N_output_rep,
                                    int init_N_non_output_act,
                                    int init_N_non_output_rep,
                                    RngStream RS)
{
    int i, j;
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
        /* mRNA decay */
        /* mRNA decay */
        genotype->mRNA_decay_rate[i] = pow(10.0,SD_MRNA_DECAY_RATE*gasdev(RS)+MEAN_MRNA_DECAY_RATE);
        if(genotype->mRNA_decay_rate[i]>MAX_MRNA_DECAY)
            genotype->mRNA_decay_rate[i]=MAX_MRNA_DECAY;
        if(genotype->mRNA_decay_rate[i]<MIN_MRNA_DECAY)
            genotype->mRNA_decay_rate[i]=MIN_MRNA_DECAY;
        /* protein decay */
        genotype->protein_decay_rate[i] = pow(10.0,SD_PROTEIN_DECAY_RATE*gasdev(RS)+MEAN_PROTEIN_DECAY_RATE);  
        if(genotype->protein_decay_rate[i]>MAX_PROTEIN_DECAY)
            genotype->protein_decay_rate[i]=MAX_PROTEIN_DECAY;
        if(genotype->protein_decay_rate[i]<MIN_PROTEIN_DECAY)
            genotype->protein_decay_rate[i]=MIN_PROTEIN_DECAY;
        /* translation rate */
        genotype->protein_syn_rate[i] = pow(10.0,SD_PROTEIN_SYN_RATE*gasdev(RS)+MEAN_PROTEIN_SYN_RATE);  
        if(genotype->protein_syn_rate[i]>MAX_PROTEIN_SYN_RATE)
            genotype->protein_syn_rate[i]=MAX_PROTEIN_SYN_RATE;
        if(genotype->protein_syn_rate[i]<MIN_PROTEIN_SYN_RATE)
            genotype->protein_syn_rate[i]=MIN_PROTEIN_SYN_RATE;
        /*ACT to INT rate*/
        genotype->active_to_intermediate_rate[i]=pow(10.0,SD_ACT_TO_INT_RATE*gasdev(RS)+MEAN_ACT_TO_INT_RATE);  
        if(genotype->active_to_intermediate_rate[i]>MAX_ACT_TO_INT_RATE)
            genotype->active_to_intermediate_rate[i]=MAX_ACT_TO_INT_RATE;
        if(genotype->active_to_intermediate_rate[i]<MIN_ACT_TO_INT_RATE)
            genotype->active_to_intermediate_rate[i]=MIN_ACT_TO_INT_RATE;        
        /*locus length*/
        genotype->locus_length[i]=(int)round(pow(10.0,SD_GENE_LENGTH*gasdev(RS)+MEAN_GENE_LENGTH));
        if(genotype->locus_length[i]>MAX_GENE_LENGTH)
            genotype->locus_length[i]=MAX_GENE_LENGTH;
        if(genotype->locus_length[i]<MIN_GENE_LENGTH)
            genotype->locus_length[i]=MIN_GENE_LENGTH;        
        genotype->total_loci_length+=genotype->locus_length[i];       
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
        genotype->protein_identity[i]=ACTIVATOR;
        genotype->is_output[i]=NON_OUTPUT_PROTEIN;
    }
    for(;i<init_N_non_output_act+N_SIGNAL_TF+init_N_non_output_rep;i++)
    {
        genotype->protein_identity[i]=REPRESSOR;
        genotype->is_output[i]=NON_OUTPUT_PROTEIN;
    }
    j=0;
    for(;i<init_N_non_output_rep+init_N_non_output_act+init_N_output_act+N_SIGNAL_TF;i++)
    {
        genotype->protein_identity[i]=ACTIVATOR;
        genotype->is_output[i]=OUTPUT_PROTEIN;
        genotype->output_gene_ids[j]=i;
        j++;
    }
    for(;i<init_N_non_output_rep+init_N_non_output_act+init_N_output_act+init_N_output_rep+N_SIGNAL_TF;i++)
    {
        genotype->protein_identity[i]=REPRESSOR;
        genotype->is_output[i]=OUTPUT_PROTEIN;
        genotype->output_gene_ids[j]=i;
        j++;
    }
    genotype->N_act=init_N_non_output_act+init_N_output_act;
    genotype->N_rep=init_N_non_output_rep+init_N_output_rep;
    /* parameterize sensor TF*/
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        genotype->mRNA_decay_rate[i]=0.0; // we assume environmental signal toggles the state of sensor TF between active and inactive 
        genotype->protein_decay_rate[i]=0.0; // the concentration of sensor TF is constant.
        genotype->protein_syn_rate[i]=0.0;
        genotype->active_to_intermediate_rate[i]=0.0; 
        genotype->protein_identity[i]=ACTIVATOR; /*make sensor TF an activator*/
        genotype->is_output[i]=NON_OUTPUT_PROTEIN;
        genotype->N_act++;
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
    } 
}

/*
 * Set how the environmental signal should change
 */
static void set_signal(CellState *state, Environment *env, float t_burn_in, RngStream RS, int thread_ID)
{
    float t=0.0;     
    char flag;   
 
    env->external_signal=NULL; // signal is not specified by an external file
    
    if(env->external_signal==NULL)   
    {
        /*always start a burn-in with signal off*/       
        state->protein_number[N_SIGNAL_TF-1]=env->signal_strength_stage1;
        state->gene_specific_protein_number[N_SIGNAL_TF-1]=env->signal_strength_stage1;
        if(t_burn_in!=0.0)
            t=t+t_burn_in; //the completion of burn_in is a fixed event, which is added in initialize_cell       
        /*after burn-in, signal should be turned "o"n*/
        flag='f'; 
    
        while(t<env->t_development+t_burn_in)
        {
            if(flag=='o')
            {                
                if(env->t_stage2!=0.0) 
                {
                    /*add a fixed event to start stage 1 again (We never do this).*/
                    add_fixed_event(-1,t+env->t_stage2,&(state->signal_off_head),&(state->signal_off_tail));
                    t=t+env->t_stage2; 
                }
                flag='f';                                  
            }    
            else
            {
                if(env->t_stage1!=0.0)
                {
                    /*add when to stage 2*/
                    add_fixed_event(-1,t+env->t_stage1,&(state->signal_on_head),&(state->signal_on_tail));
                    t=t+env->t_stage1;
                }
                flag='o';                
            }
        } 
    }
    else
    {
        int time_point=1;
        state->protein_number[N_SIGNAL_TF-1]=env->external_signal[0];
        t=1.0;
        while(t<env->t_development)
        {
            add_fixed_event(time_point,t,&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            time_point++;
            t+=1.0;
        } 
    }
}

/*copy genotype from the acestor to offsprings*/
static void clone_genotype(Genotype *genotype_templet, Genotype *genotype_clone)
{
    int i, j;           
    /*reset which_cluster for the clone*/
    for(i=0;i<MAX_GENES;i++)
        genotype_clone->which_cluster[i]=NA;
    /*copy which_cluster and cis-reg sequence*/
    for(i=0;i< genotype_templet->ngenes;i++)
    {
        genotype_clone->which_cluster[i]=genotype_templet->which_cluster[i];            
        memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_templet->cisreg_seq[i][0],CISREG_LEN*sizeof(char));                    
        genotype_clone->recalc_TFBS[i]=YES;                
    }    
    /*reset clone's cisreg_cluster*/
    for(i=N_SIGNAL_TF;i<MAX_GENES;i++)
    {
        genotype_clone->cisreg_cluster_pool[i][0][0]=0;
        for(j=0;j<MAX_GENES;j++)
            genotype_clone->cisreg_cluster_pool[i][1][j]=NA;
    }
    /*then copy from templet*/
    for(i=N_SIGNAL_TF;i<genotype_templet->N_cisreg_clusters;i++)
    {
        genotype_clone->cisreg_cluster_pool[i][0][0]=genotype_templet->cisreg_cluster_pool[i][0][0];
        for(j=0;j<genotype_templet->cisreg_cluster_pool[i][0][0];j++)
            genotype_clone->cisreg_cluster_pool[i][1][j]=genotype_templet->cisreg_cluster_pool[i][1][j];
    }
       
    /*reset clone's information*/
    for(i=0;i<MAX_GENES;i++)
    {
        genotype_clone->which_protein[i]=NA;
        genotype_clone->min_N_activator_to_transc[i]=MAX_BINDING+1;
        genotype_clone->is_output[i]=NON_OUTPUT_PROTEIN;
        genotype_clone->which_node_family[i]=NA;
    }    
    /*reset clone's protein_pool*/
    for(i=0;i<MAX_PROTEINS;i++)
    {            
        for(j=0;j<MAX_GENES;j++)
            genotype_clone->protein_pool[i][1][j]=NA;
        genotype_clone->protein_pool[i][0][0]=0;            
    }
    /*reset clone's output_protein_idss*/
    for(i=0;i<MAX_OUTPUT_GENES;i++)
        genotype_clone->output_gene_ids[i]=NA;
    /*reset clone's node_family_pool and which_node_family*/
    for(i=0;i<MAX_PROTEINS;i++)
    {        
        for(j=0;j<MAX_GENES;j++)
            genotype_clone->node_family_pool[i][1][j]=NA;
        genotype_clone->node_family_pool[i][0][0]=0;      
    }
    /*copy from templet's node_family_pool*/
    for(i=0;i<genotype_templet->N_node_families;i++)
    {
        genotype_clone->node_family_pool[i][0][0]=genotype_templet->node_family_pool[i][0][0];
        for(j=0;j<genotype_templet->node_family_pool[i][0][0];j++)
            genotype_clone->node_family_pool[i][1][j]=genotype_templet->node_family_pool[i][1][j];
    }
    /*copy from templet's which_node_family and is_output*/
    for(i=0;i<genotype_templet->ngenes;i++)
    {
        genotype_clone->which_node_family[i]=genotype_templet->which_node_family[i];
        genotype_clone->is_output[i]=genotype_templet->is_output[i];
    }
    /*copy from templet's protein_pool*/
    for(i=0;i<genotype_templet->nproteins;i++)
    {            
        genotype_clone->protein_pool[i][0][0]=genotype_templet->protein_pool[i][0][0];        
        for(j=0;j<genotype_templet->protein_pool[i][0][0];j++)
            genotype_clone->protein_pool[i][1][j]=genotype_templet->protein_pool[i][1][j];                     
    }   
    /*copy from templet's output_protein_ids*/
    for(i=0;i<genotype_templet->n_output_genes;i++)
        genotype_clone->output_gene_ids[i]=genotype_templet->output_gene_ids[i];           
    /* copy binding sites' sequences*/  
    for(i=0;i<genotype_templet->ngenes;i++) 
    {          
        for(j=0;j<TF_ELEMENT_LEN;j++)
        {    
            genotype_clone->tf_binding_seq[i][j]=genotype_templet->tf_binding_seq[i][j];
            genotype_clone->tf_binding_seq_rc[i][j]=genotype_templet->tf_binding_seq_rc[i][j];
        }
    }
    /*reset locus_specific_TF_behavior*/
    for(i=0;i<MAX_GENES;i++)
        for(j=0;j<MAX_PROTEINS;j++)
            genotype_clone->locus_specific_TF_behavior[i][j]=NON_TF;
    /*copy kinetic constants*/
    for(i=0;i<genotype_templet->ngenes;i++) 
    {            
        genotype_clone->mRNA_decay_rate[i]=genotype_templet->mRNA_decay_rate[i];
        genotype_clone->protein_decay_rate[i]=genotype_templet->protein_decay_rate[i];
        genotype_clone->protein_syn_rate[i]=genotype_templet->protein_syn_rate[i];            
        genotype_clone->active_to_intermediate_rate[i]=genotype_templet->active_to_intermediate_rate[i];
        genotype_clone->which_protein[i]=genotype_templet->which_protein[i];
        genotype_clone->min_N_activator_to_transc[i]=genotype_templet->min_N_activator_to_transc[i];  
        genotype_clone->locus_length[i]=genotype_templet->locus_length[i];
        for(j=0;j<genotype_templet->nproteins;j++)
            genotype_clone->locus_specific_TF_behavior[i][j]=genotype_templet->locus_specific_TF_behavior[i][j];
    } 
    /* copy TF information*/
    for(i=0;i<MAX_PROTEINS;i++)
    {
        genotype_clone->protein_identity[i]=genotype_templet->protein_identity[i];
        genotype_clone->Kd[i]=genotype_templet->Kd[i];
    }    
    /* copy gene and protein numbers*/
    genotype_clone->ngenes=genotype_templet->ngenes; 
    genotype_clone->nproteins=genotype_templet->nproteins;
    genotype_clone->n_output_genes=genotype_templet->n_output_genes;
    genotype_clone->N_node_families=genotype_templet->N_node_families;
    genotype_clone->N_cisreg_clusters=genotype_templet->N_cisreg_clusters;
    genotype_clone->N_act=genotype_templet->N_act;
    genotype_clone->N_rep=genotype_templet->N_rep;
    genotype_clone->total_loci_length=genotype_templet->total_loci_length;
    genotype_clone->flag_effector_is_TF=genotype_templet->flag_effector_is_TF;
    genotype_clone->flag_rm_which_TFBS=genotype_templet->flag_rm_which_TFBS;
    genotype_clone->flag_signal_ctrl_repressor=genotype_templet->flag_signal_ctrl_repressor;
}


/***/
static void run_fitness_test(Genotype *genotype,
                            CellState *state,
                            GillespieRates *rate,
                            Environment *Env,
                            Phenotype *timecourse,
                            float *f,
                            int *mRNA,
                            float *protein,
                            int N_replicates_per_thread,
                            RngStream RS,
                            int thread_ID)
{
    int i;
    float t_burn_in;
    
    for(i=0;i<N_replicates_per_thread;i++) /* env 1, usually a constant signal that matches env*/
    {  
        /*make a t_burn_in before turning on signal*/
        do
            t_burn_in=Env->avg_duration_of_burn_in_growth_rate*expdev(RS);
        while(t_burn_in>Env->max_duration_of_burn_in_growth_rate);
        /*initialize mRNA and protein numbers, and gene states etc.*/
        initialize_cell(genotype, state, Env, t_burn_in, mRNA, protein);

        /*set how the signal should change during simulation*/
        set_signal(state, Env, t_burn_in, RS, thread_ID);

        /*calcualte the rates of cellular activity based on the initial cellular state*/
        calc_all_rates(genotype, state, rate, Env, timecourse, t_burn_in, INITIALIZATION);             

        /*run developmental simulation until tdevelopment or encounter an error*/
        while(state->t<Env->t_development+t_burn_in) 
            do_single_timestep(genotype, state, rate, Env, t_burn_in, &(timecourse[i]), RS);            
#if !PHENOTYPE
        f[i]=calc_replicate_fitness(state,Env);
#endif
        /*free linked tables*/
        free(state->sampled_response);
        free_fixedevent(state);           
    }  
}

static void copy_environment(Environment *env_source, Environment *env_copy)
{
    env_copy->max_t_mid=env_source->max_t_mid;
    env_copy->t_development=env_source->t_development;
    env_copy->signal_strength_stage1=env_source->signal_strength_stage1;
    env_copy->signal_strength_stage2=env_source->signal_strength_stage2;
    env_copy->t_stage1=env_source->t_stage1;
    env_copy->t_stage2=env_source->t_stage2;
    env_copy->initial_effect_of_effector=env_source->initial_effect_of_effector;
    env_copy->effect_of_effector_aft_burn_in=env_source->effect_of_effector_aft_burn_in;
    env_copy->fixed_effector_effect=env_source->fixed_effector_effect;
    env_copy->max_duration_of_burn_in_growth_rate=env_source->max_duration_of_burn_in_growth_rate;
    env_copy->avg_duration_of_burn_in_growth_rate=env_source->avg_duration_of_burn_in_growth_rate; 
    env_copy->min_reduction_relative_to_peak=env_source->min_reduction_relative_to_peak;
    env_copy->fitness_decay_constant=env_source->fitness_decay_constant;
    env_copy->effector_level_before_stage2=env_source->effector_level_before_stage2;
    env_copy->window_size=env_source->window_size;  
    env_copy->opt_peak_response=env_source->opt_peak_response;
    env_copy->w1=env_source->w1;
    env_copy->w2=env_source->w2;
    env_copy->w3=env_source->w3;
    env_copy->w4=env_source->w4;
    env_copy->is_burn_in=env_source->is_burn_in;
}

/**
 *Calculate the fintess of a given genotype.
 *Essentially calling do_single_timestep until tdevelopment and calculate 
 *average growth rate over tdevelopment.
 */
static void calc_avg_fitness( Genotype *genotype,
                                Selection *selection, 
                                Phenotype phenotype[N_REPLICATES],
                                int init_mRNA[MAX_GENES],
                                float init_protein_number[MAX_PROTEINS],
                                RngStream RS_parallel[N_THREADS],           
                                float Fitness1[N_REPLICATES],
                                float Fitness2[N_REPLICATES],                        
                                int RM_WHAT)         
{        
    /*Making clones of a genotype, and have the clones run in parallel*/
    #pragma omp parallel num_threads(N_THREADS) 
    {
        int thread_ID=omp_get_thread_num();   
        int i,j,k;
        int N_replicates_per_thread=N_REPLICATES/N_THREADS;  
        Genotype genotype_clone;
        CellState state_clone;
        GillespieRates rate_clone;
        int init_mRNA_clone[MAX_GENES];      
        float init_protein_number_clone[MAX_GENES];
        float f1[N_replicates_per_thread];

        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];         
        Selection selection_copy;

        /*alloc space for linked tables and set default values for parameters, in genotype*/
        initialize_cache(&genotype_clone);  
        
        /*clone genotype and initial mRNA and protein numbers*/
        #pragma omp critical
        {
            genotype_clone.ngenes=genotype->ngenes;
            genotype_clone.n_output_genes=genotype->n_output_genes;         
            genotype_clone.nproteins=genotype->nproteins;
            genotype_clone.N_node_families=genotype->N_node_families;
            genotype_clone.N_cisreg_clusters=genotype->N_cisreg_clusters;            
            clone_genotype(genotype, &genotype_clone);  
            copy_environment(&(selection->env1),&(selection_copy.env1));

            for(j=0; j < MAX_GENES; j++) 
            {  
                init_mRNA_clone[j] = init_mRNA[j];
                init_protein_number_clone[j] = init_protein_number[j];
            } 
        } 

        genotype_clone.flag_rm_which_TFBS=RM_WHAT;
        calc_all_binding_sites(&genotype_clone); 

        /*Set initial mRNA and protein number using given values*/
        for(j=N_SIGNAL_TF; j < genotype_clone.ngenes; j++)        
            mRNA[j] = init_mRNA_clone[j];                       
        for(j=N_SIGNAL_TF; j<genotype_clone.nproteins;j++)
        {
            for(k=0;k<genotype_clone.protein_pool[j][0][0];k++)
                protein[genotype_clone.protein_pool[j][1][k]]=(float)init_protein_number_clone[j]/genotype_clone.protein_pool[j][0][0]; //split the initial protein number equally to different copies
                                                                                                                                        //this is to make sure all proteins have equal initial numbers        
        }    
        
        /* now calc growth rate under two environments*/
        /********************************************************************** 
         * 
         *                              TEST1: env1 usually a constant signal that matches env 
         *
         *********************************************************************/
        run_fitness_test(&genotype_clone, 
                        &state_clone, 
                        &rate_clone, 
                        &(selection_copy.env1), 
                        &(phenotype[thread_ID*N_replicates_per_thread]),
                        &(f1[0]),
                        &(mRNA[0]),
                        &(protein[0]),
                        N_replicates_per_thread,
                        RS_parallel[thread_ID],
                        thread_ID);      

        /*free linked tables*/
        for(j=0;j<MAX_GENES;j++)
            free(genotype_clone.all_binding_sites[j]);
            
        /*pool fitness from each thread*/
        #pragma omp critical
        {
            j=0;
            for(i=thread_ID*N_replicates_per_thread;i<(thread_ID+1)*N_replicates_per_thread;i++)
            {
                Fitness1[i]=f1[j];
                Fitness2[i]=0.0;            
                j++;
            }
        }
    }  
}

static float calc_replicate_fitness(CellState *state, Environment *env)
{
    float fitness; 
    float avg_response_bf_peak, avg_response_aft_peak; 
    float peak_response, half_response; 
    float max_response_before_peak, max_response_after_peak;
    int pos_peak;   
    float pos_half_response_bf_peak;   
    int i;
    float temp;
    
    fitness=0.0;
    /*find peak*/
    find_max(&(state->sampled_response[0]),env->window_size-1,state->N_samples,&peak_response,&pos_peak,-1.0); 

    /*calculate the average signal strength at the end of simulation*/
    avg_response_aft_peak=0.0;
    for(i=state->N_samples-1;i>state->N_samples-1-env->window_size;i--)
        avg_response_aft_peak+=state->sampled_response[i];
    avg_response_aft_peak/=(float)env->window_size;   

    /*calculate the average signal strength before the signal change*/
    avg_response_bf_peak=0.0;
    for(i=0;i<env->window_size;i++)
        avg_response_bf_peak+=state->sampled_response[i];
    avg_response_bf_peak/=(float)env->window_size; 
       
    if(peak_response==0.0 || peak_response<=state->sampled_response[env->window_size-1])
        fitness=0.0-exp_cost_factor*state->cumulative_cost; 
    else
    {
        /*select for peak level*/
        fitness+=env->w1*exp(-log(peak_response/env->opt_peak_response)*log(peak_response/env->opt_peak_response)/env->fitness_decay_constant);

        /* look for mid point bf the peak*/
        half_response=(peak_response+state->sampled_response[env->window_size-1])*0.5;      
        find_x(&(state->sampled_response[0]),env->window_size-1,pos_peak,half_response,&pos_half_response_bf_peak,0); 

        /*faster is better*/
        temp=(env->t_development-env->t_stage1+(float)env->window_size-pos_half_response_bf_peak-1.0)/(env->t_development-env->t_stage1-env->max_t_mid);
        temp=(temp>1.0)?1.0:temp;
        fitness+=env->w2*temp;  

        /*select for low ss level bf signal change*/
        max_response_before_peak=env->effector_level_before_stage2*peak_response;                 
        temp=(avg_response_bf_peak<max_response_before_peak)?1.0:(peak_response-avg_response_bf_peak)/(peak_response-max_response_before_peak);
        fitness+=env->w3*temp;

        /*low ss level after signal change*/
        max_response_after_peak=(1.0-env->min_reduction_relative_to_peak)*peak_response;
        temp=(avg_response_aft_peak<max_response_after_peak)?1.0:(peak_response-avg_response_aft_peak)/(peak_response-max_response_after_peak);
        fitness+=env->w4*temp;

        /*adds up*/
        fitness=fitness/(env->w1+env->w2+env->w3+env->w4)-exp_cost_factor*state->cumulative_cost;               
    }
    return fitness;
}

/**
 * Given the fitness of the resident and a mutant, decide whether the mutant can replace the resident
 */
static void try_replacement(Genotype *resident, Genotype *mutant, int *flag_replaced, float *selection_coefficient)
{     
    *selection_coefficient=(mutant->avg_fitness-resident->avg_fitness)/fabs(resident->avg_fitness);
    if(*selection_coefficient>=MIN_SELECTION_COEFFICIENT)
        *flag_replaced=1;
    else          
        *flag_replaced=0;    
}

static void replay_mutations(Genotype *resident, Mutation *mut_record, int replay_N_steps)
{
    int i, j, output_counter;
    Output_buffer resident_info[OUTPUT_INTERVAL];
    FILE *fp;

    /*load mutation record*/
    fp=fopen(mutation_file,"r");    
    if(fp!=NULL)        
        printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
    else
    {
        printf("Loading mutation record failed! Quit program!");
#if MAKE_LOG
        LOG("Loading mutation record failed!");
#endif
        exit(-2);
    }
    
    /*remove the old file*/
    remove("networks.txt");  
    remove("N_motifs.txt");
    calc_all_binding_sites(resident);
    summarize_binding_sites(resident,0); //make new file and record initial network
      
    output_counter=0;   
    for(i=1;i<=replay_N_steps;i++)
    {  
        fscanf(fp,"%c %d %d %d %s %d %a\n",&(mut_record->mut_type),
                                            &(mut_record->which_gene),                                                    
                                            &(mut_record->which_nucleotide), 
                                            &(mut_record->which_protein),
                                            mut_record->nuc_diff,               
                                            &(mut_record->kinetic_type),
                                            &(mut_record->kinetic_diff));
        reproduce_mutate(resident,mut_record);
        calc_all_binding_sites(resident);
        find_motifs(resident);        
        store_resident_info(resident,NULL,&(resident_info[output_counter]),NA,NA,NA,(float)NA,-1);  
        output_counter++;
        if(i%OUTPUT_INTERVAL==0)
        {
            summarize_binding_sites(resident,i);
            output_resident_info(resident_info, output_counter,-1);
            output_counter=0;
        }        
    }
    
    /*call output_resident_info again, just in case buffer stores less than OUTPUT_INTERVAL cycles*/
    if(output_counter<OUTPUT_INTERVAL && output_counter!=0)
        output_resident_info(resident_info, output_counter, -1);
    
    printf("Reproduce mutations successfully!\n");
    fclose(fp);
}

#if PERTURB
static void modify_network(Genotype *resident,
                    Genotype *mutant, 
                    Mutation *mut_record,  
                    Selection *selection,                   
                    int init_mRNA[MAX_GENES],
                    float init_protein[MAX_PROTEINS],
                    RngStream RS_parallel[N_THREADS])
{
    int i,j;  
    char buffer[600],char_buffer,char_buffer2[2];
    int int_buffer;
    float float_buffer;
    float fitness_after_perturbation[10],fitness_bf_perturbation,se_fitness_bf_perturbation;
    FILE *file_mutation,*f_aft_perturbation,*fitness_record;
    
    /*remove previous output*/  
    remove("N_motifs.txt");
    remove("networks.txt");
    
    /*load mutation record*/
    file_mutation=fopen(mutation_file,"r");    
    if(file_mutation!=NULL)        
        printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
    else
    {
        printf("Loading mutation record failed! Quit program!");
#if LOG
        LOG("Loading mutation record failed!");
#endif
        exit(-2);
    }     
    
    /*load the original fitness*/
    fitness_record=fopen(evo_summary,"r");
    if(fitness_record!=NULL)        
        printf("LOAD FITNESS SUCCESSFUL!\n");
    else
    {
        printf("Loading fitness record failed! Quit program!");
#if LOG
        LOG("Loading fitness record failed!");
#endif
        exit(-2);
    } 
    fgets(buffer,600,fitness_record); // skip the header    
    
    /*create threads*/
    omp_set_num_threads(N_THREADS);
        
    /*determine non-adaptive weak TFBSs in the initial genotype*/
    fscanf(fitness_record,"%d %d %d %d %s %s %f %f %f %f %f %f %d %d %d\n",
            &int_buffer,
            &int_buffer,
            &int_buffer,
            &int_buffer,
            &(char_buffer2[0]),
            &(char_buffer2[0]),
            &fitness_bf_perturbation,     // load the original fitness           
            &float_buffer,
            &float_buffer,
            &se_fitness_bf_perturbation, // and the se of fitness
            &float_buffer,
            &float_buffer,
            &int_buffer,
            &int_buffer,
            &int_buffer);  
    
    clone_genotype(resident,mutant); // do perturbation on a copy of the initial genotype
    mutant->flag_rm_which_TFBS=determine_adaptive_TFBSs(mutant,
                                                        selection,
                                                        fitness_bf_perturbation,
                                                        se_fitness_bf_perturbation,
                                                        init_mRNA,
                                                        init_protein,
                                                        RS_parallel,
                                                        HI_RESOLUTION_RECALC,
                                                        fitness_after_perturbation,
                                                        CUT_OFF_NONADAPTIVE_TFBS); 
    /*output the result*/
    f_aft_perturbation=fopen("after_perturbation.txt","a+");
    fprintf(f_aft_perturbation,"0 %d ", mutant->flag_rm_which_TFBS);            
    for(j=0;j<10;j++)           
        fprintf(f_aft_perturbation,"%.10f ", fitness_after_perturbation[j]);                               
    fprintf(f_aft_perturbation,"\n");
    fclose(f_aft_perturbation); 

    /*use the perturbation result to exclude non-adaptive 2-mismatch TFBSs and count motifs*/
    calc_all_binding_sites(mutant);
    find_motifs(mutant);
    print_motifs(mutant);
    summarize_binding_sites(mutant,0); 
    
    /*repeat above for the evolved genotypes*/
    for(i=1;i<=selection->MAX_STEPS;i++)
    {  
        /*reproduce mutation*/
        clone_genotype(resident,mutant);
        fscanf(file_mutation,"%c %d %d %d %s %d %a\n",&(mut_record->mut_type),
                                                        &(mut_record->which_gene),
                                                        &(mut_record->which_nucleotide), 
                                                        &(mut_record->which_protein),
                                                        mut_record->nuc_diff,               
                                                        &(mut_record->kinetic_type),
                                                        &(mut_record->kinetic_diff));
        reproduce_mutate(mutant,mut_record); 
        clone_genotype(mutant,resident);
        
        /*load the orignal fitness*/
        fscanf(fitness_record,"%d %d %d %d %c %f %f %f %f %f %f %f %d %d %d\n",
                &int_buffer,
                &int_buffer,
                &int_buffer,
                &int_buffer,
                &char_buffer,
                &float_buffer,
                &fitness_bf_perturbation,                
                &float_buffer,
                &float_buffer,
                &se_fitness_bf_perturbation,
                &float_buffer,
                &float_buffer,
                &int_buffer,
                &int_buffer,
                &int_buffer);        
       
        /*perturb and determine adaptive 2-mismatch TFBSs*/
        mutant->flag_rm_which_TFBS=determine_adaptive_TFBSs(mutant,
                                                            selection,
                                                            fitness_bf_perturbation,
                                                            se_fitness_bf_perturbation,
                                                            init_mRNA,
                                                            init_protein,
                                                            RS_parallel,
                                                            HI_RESOLUTION_RECALC,
                                                            fitness_after_perturbation,
                                                            CUT_OFF_NONADAPTIVE_TFBS);

        /*output the result*/
        f_aft_perturbation=fopen("after_perturbation.txt","a+");
        fprintf(f_aft_perturbation,"%d %d ",i, mutant->flag_rm_which_TFBS);            
        for(j=0;j<10;j++)           
            fprintf(f_aft_perturbation,"%.10f ", fitness_after_perturbation[j]);                               
        fprintf(f_aft_perturbation,"\n");
        fclose(f_aft_perturbation); 

        /*use the perturbation result to exclude non-adaptive 2-mismatch TFBSs and count motifs*/
        for(j=0;j<mutant->ngenes;j++)
            mutant->recalc_TFBS[j]=YES;
        calc_all_binding_sites(mutant);
        find_motifs(mutant);
        print_motifs(mutant);
        if(i%OUTPUT_INTERVAL==0 && i!=0)
            summarize_binding_sites(mutant,i);       
    } 
    rename("networks.txt","networks_clean.txt");
    rename("N_motifs.txt","N_motifs_clean.txt");
    fclose(file_mutation);
    fclose(fitness_record);
}

static void classify_all_mutations(Genotype *resident, 
                            Mutation *mut_record,  
                            Selection *selection, 
                            int max_evolutionary_steps,
                            int init_mRNA[MAX_GENES], 
                            float init_protein[MAX_GENES], 
                            RngStream RS_parallel[N_THREADS])
{
    int i, j, k, output_counter;
    FILE *fp, *fp_all_mut, *fp_all_mut_fitness, *fp_which_TFBS, *fp_output; 
    int which_TFBS;
    int int_buffer,evo_step,evo_step_copy,N_total_mutations;
    char buffer[600],char_buffer;
    float float_buffer;
    float initial_fitness,initial_fitness_se;
    Mutation mut_record_copy;
    Genotype resident_copy;
    initialize_cache(&resident_copy);
    int N_I1_bf_mut,N_I1_aft_mut,N_NFBL_bf_mut,N_NFBL_aft_mut;
    int N_II_bf_mut,N_II_aft_mut,N_IN_bf_mut,N_IN_aft_mut;
    int output_buffer[50][6];
    float output_buffer2[50][10];   
     
    /*load mutation record*/    
    fp=fopen(evo_summary,"r");
    for(i=0;i<(max_evolutionary_steps+1);i++)   
        fgets(buffer,600,fp); // skip header, the initial genotype, and 1 to (max_evolutionary_steps-1) evolved genotype 
    
    /* get total number of mutations that were trialed during the first max_evolutionary_steps */
    fscanf(fp,"%d %d %d %d %c %f %f %f %f %f %f %f %d %d %d\n",
                &int_buffer,
                &N_total_mutations, 
                &int_buffer,
                &int_buffer,
                &char_buffer,
                &float_buffer,
                &float_buffer,
                &float_buffer,
                &float_buffer,
                &float_buffer,
                &float_buffer,
                &float_buffer,            
                &int_buffer,
                &int_buffer,
                &int_buffer); 
    fclose(fp);
          
    /*load other records*/
    fp_all_mut=fopen("all_mutations.txt","r");
    fp_all_mut_fitness=fopen("fitness_all_mutants.txt","r");
    fp_which_TFBS=fopen("after_perturbation.txt","r");   // needs to run perturbation analysis to generate this file
    
    /*set threads*/
    omp_set_num_threads(N_THREADS);  

    /*determine adaptive motifs in the initial network*/
    fscanf(fp_which_TFBS,"%d %d %f %f %f %f %f %f %f %f %f %f \n",
                        &int_buffer,                
                        &(resident->flag_rm_which_TFBS), // load non-adaptive TFBSs
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer, 
                        &float_buffer,
                        &float_buffer,
                        &float_buffer); 
    for(j=0;j<resident->ngenes;j++)
        resident->recalc_TFBS[j]=YES;
    calc_all_binding_sites(resident);
    find_motifs(resident); 
    N_I1_bf_mut=resident->N_motifs[3]+resident->N_motifs[4]+resident->N_motifs[11]+resident->N_motifs[12];           
    N_NFBL_bf_mut=resident->N_motifs[5]+resident->N_motifs[6]+resident->N_motifs[13]+resident->N_motifs[14];   
    N_II_bf_mut=resident->N_motifs[7]+resident->N_motifs[8]+resident->N_motifs[15]+resident->N_motifs[16];           
    N_IN_bf_mut=resident->N_motifs[1]+resident->N_motifs[2]+resident->N_motifs[9]+resident->N_motifs[10]; 
    restore_all_TFBS(resident);   
    
    /*reset output counter*/
    output_counter=0;    
    
    /*load the first mutation*/
    fscanf(fp_all_mut,"%d %d %c %d %d %d %s %d %a\n",&evo_step,
                                            &int_buffer,
                                            &(mut_record->mut_type),
                                            &(mut_record->which_gene),                                                    
                                            &(mut_record->which_nucleotide), 
                                            &(mut_record->which_protein),
                                            mut_record->nuc_diff,               
                                            &(mut_record->kinetic_type),
                                            &(mut_record->kinetic_diff));
    evo_step_copy=evo_step;
    mut_record_copy.mut_type=mut_record->mut_type;
    mut_record_copy.which_gene=mut_record->which_gene;
    mut_record_copy.which_nucleotide=mut_record->which_nucleotide;
    mut_record_copy.which_protein=mut_record->which_protein;
    mut_record_copy.nuc_diff[0]=mut_record->nuc_diff[0];
    mut_record_copy.nuc_diff[1]=mut_record->nuc_diff[1];
    mut_record_copy.nuc_diff[2]=mut_record->nuc_diff[2];
    mut_record_copy.kinetic_type=mut_record->kinetic_type;
    mut_record_copy.kinetic_diff=mut_record->kinetic_diff; 
    
    /*replay mutations*/
    for(i=1;i<=N_total_mutations;i++)  
    { 
        if(i<N_total_mutations) // load mutation
            fscanf(fp_all_mut,"%d %d %c %d %d %d %s %d %a\n",&evo_step,
                                                            &int_buffer,
                                                            &(mut_record->mut_type),
                                                            &(mut_record->which_gene),                                                    
                                                            &(mut_record->which_nucleotide), 
                                                            &(mut_record->which_protein),
                                                            mut_record->nuc_diff,               
                                                            &(mut_record->kinetic_type),
                                                            &(mut_record->kinetic_diff));
        /*load fitness of mutant*/
        fscanf(fp_all_mut_fitness,"%f %f %f %f %f %f\n", &initial_fitness,
                                                        &float_buffer,
                                                        &float_buffer,
                                                        &initial_fitness_se,
                                                        &float_buffer,
                                                        &float_buffer);
       
        /*col 1 of all mutation is the evolutionary step which is shared by all mutations tried.
         *Here we compare the evo step (evo_step) of a mutation with that (evo_step_copy) of the last mutation.
         * 
         */
        if(evo_step!=evo_step_copy || i==N_total_mutations) // this means an accepted mutation
        {                         
            reproduce_mutate(resident,&mut_record_copy);  
            fscanf(fp_which_TFBS,"%d %d %f %f %f %f %f %f %f %f %f %f \n",
                        &int_buffer,                
                        &which_TFBS, // load non-adaptive TFBSs
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer,
                        &float_buffer, 
                        &float_buffer,
                        &float_buffer,
                        &float_buffer);  
            resident->flag_rm_which_TFBS=which_TFBS;
            for(j=0;j<resident->ngenes;j++)
                resident->recalc_TFBS[j]=YES;
            calc_all_binding_sites(resident);
            find_motifs(resident); 
            N_I1_aft_mut=(resident->N_motifs[3]+resident->N_motifs[4]+resident->N_motifs[11]+resident->N_motifs[12]);
            N_NFBL_aft_mut=(resident->N_motifs[5]+resident->N_motifs[6]+resident->N_motifs[13]+resident->N_motifs[14]);
            N_II_aft_mut=resident->N_motifs[7]+resident->N_motifs[8]+resident->N_motifs[15]+resident->N_motifs[16];           
            N_IN_aft_mut=resident->N_motifs[1]+resident->N_motifs[2]+resident->N_motifs[9]+resident->N_motifs[10]; 
            restore_all_TFBS(resident);     
            
            /*reset output buffer*/
            for(j=0;j<6;j++)
                output_buffer[output_counter][j]=0;
            
            if(N_I1_bf_mut==0 && N_I1_aft_mut>0) // I1-creating mutation
                output_buffer[output_counter][2]=1;            
            else if(N_I1_bf_mut>0 && N_I1_aft_mut==0) // I1-destroying mutation
                output_buffer[output_counter][2]=-1;   
            
            if(N_NFBL_bf_mut==0 && N_NFBL_aft_mut>0)   // NFBL-creating mutation       
                output_buffer[output_counter][3]=1;   
            else if(N_NFBL_bf_mut>0 && N_NFBL_aft_mut==0)    // NFBL-destroying mutation      
                output_buffer[output_counter][3]=-1;  
            
            if(N_II_bf_mut==0 && N_II_aft_mut>0) // overlapping-I1-creating mutation
                output_buffer[output_counter][4]=1;            
            else if(N_II_bf_mut>0 && N_II_aft_mut==0) // overlapping-I1-destroying mutation
                output_buffer[output_counter][4]=-1;   
            
            if(N_IN_bf_mut==0 && N_IN_aft_mut>0)   // I+N-conjugate-creating mutation       
                output_buffer[output_counter][5]=1;   
            else if(N_IN_bf_mut>0 && N_IN_aft_mut==0)    // I+N-conjugate-destroying mutation      
                output_buffer[output_counter][5]=-1;  
            
            N_I1_bf_mut=N_I1_aft_mut;
            N_NFBL_bf_mut=N_NFBL_aft_mut;
            N_II_bf_mut=N_II_aft_mut;
            N_IN_bf_mut=N_IN_aft_mut;
            output_buffer[output_counter][0]=1;    // means accepted mutation 
            output_buffer[output_counter][1]=which_TFBS; 
            /*already recorded the effect of perturbing accepted mutations in purbation analysis, so set fitness record to -1.0 here*/
            for(j=0;j<10;j++)                     
                output_buffer2[output_counter][j]=-1.0;           
        }
        else // mutations not accepted
        {            
            clone_genotype(resident,&resident_copy);    
            reproduce_mutate(&resident_copy,&mut_record_copy);
            which_TFBS=determine_adaptive_TFBSs(&resident_copy,
                                                selection,      
                                                initial_fitness,
                                                initial_fitness_se,
                                                init_mRNA,
                                                init_protein,  
                                                RS_parallel,
                                                LOW_RESOLUTION_RECALC,
                                                output_buffer2[output_counter],
                                                CUT_OFF_NONADAPTIVE_TFBS_2); // use a different cutoff because fitness resolution is low     
            for(j=0;j<resident_copy.ngenes;j++)
                resident_copy.recalc_TFBS[j]=YES;
            calc_all_binding_sites(&resident_copy);
            find_motifs(&resident_copy); 
            N_I1_aft_mut=(resident_copy.N_motifs[3]+resident_copy.N_motifs[4]+resident_copy.N_motifs[11]+resident_copy.N_motifs[12]);
            N_NFBL_aft_mut=(resident_copy.N_motifs[5]+resident_copy.N_motifs[6]+resident_copy.N_motifs[13]+resident_copy.N_motifs[14]);  
            N_II_aft_mut=(resident_copy.N_motifs[7]+resident_copy.N_motifs[8]+resident_copy.N_motifs[15]+resident_copy.N_motifs[16]);
            N_IN_aft_mut=(resident_copy.N_motifs[1]+resident_copy.N_motifs[2]+resident_copy.N_motifs[9]+resident_copy.N_motifs[10]);

            for(j=0;j<6;j++)
                output_buffer[output_counter][j]=0;

            if(N_I1_bf_mut==0 && N_I1_aft_mut>0)
                output_buffer[output_counter][2]=1;            
            else if(N_I1_bf_mut>0 && N_I1_aft_mut==0) 
                output_buffer[output_counter][2]=-1;   

            if(N_NFBL_bf_mut==0 && N_NFBL_aft_mut>0)          
                output_buffer[output_counter][3]=1;   
            else if(N_NFBL_bf_mut>0 && N_NFBL_aft_mut==0)          
                output_buffer[output_counter][3]=-1;  

            if(N_II_bf_mut==0 && N_II_aft_mut>0) // overlapping-I1-creating mutation
                output_buffer[output_counter][4]=1;            
            else if(N_II_bf_mut>0 && N_II_aft_mut==0) // overlapping-I1-destroying mutation
                output_buffer[output_counter][4]=-1;   
            
            if(N_IN_bf_mut==0 && N_IN_aft_mut>0)   // I+N-conjugate-creating mutation       
                output_buffer[output_counter][5]=1;   
            else if(N_IN_bf_mut>0 && N_IN_aft_mut==0)    // I+N-conjugate-destroying mutation      
                output_buffer[output_counter][5]=-1;
            
            output_buffer[output_counter][0]=0;             
            output_buffer[output_counter][1]=which_TFBS;    
        }     
        evo_step_copy=evo_step;
        mut_record_copy.mut_type=mut_record->mut_type;
        mut_record_copy.which_gene=mut_record->which_gene;
        mut_record_copy.which_nucleotide=mut_record->which_nucleotide;
        mut_record_copy.which_protein=mut_record->which_protein;
        mut_record_copy.nuc_diff[0]=mut_record->nuc_diff[0];
        mut_record_copy.nuc_diff[1]=mut_record->nuc_diff[1];
        mut_record_copy.nuc_diff[2]=mut_record->nuc_diff[2];
        mut_record_copy.kinetic_type=mut_record->kinetic_type;
        mut_record_copy.kinetic_diff=mut_record->kinetic_diff;   
        output_counter++;

        /*output OUTPUT_INTERVAL mutations at a time*/
        if(output_counter%OUTPUT_INTERVAL==0 && output_counter!=0)
        {
            fp_output=fopen("after_perturbation_all_mutations.txt","a+");
            for(j=0;j<OUTPUT_INTERVAL;j++)   
            {   
                for(k=0;k<10;k++)                
                    fprintf(fp_output,"%.10f ", output_buffer2[j][k]);                
            }
            fclose(fp_output);
            
            fp_output=fopen("mutation_classification.txt","a+");
            for(j=0;j<OUTPUT_INTERVAL;j++)                             
                fprintf(fp_output,"%d %d %d %d %d %d\n",output_buffer[j][0],output_buffer[j][1],output_buffer[j][2],output_buffer[j][3],output_buffer[j][4],output_buffer[j][5]);
            fclose(fp_output);            

            output_counter=0;
        }  
    } 
    fclose(fp_all_mut);
    fclose(fp_all_mut_fitness);
    fclose(fp_which_TFBS);
    /*flush out the remaining record*/
    if(output_counter!=0)
    {
        fp_output=fopen("fitness_aft_perturbation2.txt","a+");
        for(j=0;j<output_counter;j++)   
        {   
            for(k=0;k<10;k++)
                fprintf(fp_output,"%.10f ", output_buffer2[j][k]);            
        }
        fclose(fp_output);

        fp_output=fopen("mutation_classification.txt","a+");
        for(j=0;j<output_counter;j++)                             
            fprintf(fp_output,"%d %d %d %d\n",output_buffer[j][0],output_buffer[j][1],output_buffer[j][2],output_buffer[j][3],output_buffer[j][4],output_buffer[j][5]);
        fclose(fp_output);        
    }   
}
#endif

static void sample_expression_lvl(Genotype *resident, 
                                    Selection *selection,                                    
                                    int init_mRNA[MAX_GENES], 
                                    float init_protein[MAX_GENES], 
                                    RngStream RS_parallel[N_THREADS])
{
    int i,j,k;
    float fitness1[N_REPLICATES],fitness2[N_REPLICATES];   
    Phenotype phenotype[N_REPLICATES];
    FILE *fp;
    char filename[32];
   
    /*initialize phenotype*/
    for (i=0;i<N_REPLICATES;i++)
    {
        phenotype[i].total_time_points=(int)(selection->env1.t_development/sampling_interval)+1;
        phenotype[i].gene_specific_concentration=(float *)calloc(phenotype[i].total_time_points,resident->ngenes*sizeof(float));
        phenotype[i].protein_concentration=(float *)calloc(phenotype[i].total_time_points,resident->N_node_families*sizeof(float));    
        phenotype[i].timepoint=0;         
    }
    
    /*generate phenotype*/    
    calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1, fitness2, RM_NONE);
  
    /*output phenotype*/    
    /*protein concentration: each protein has its own file, in which each row is a replicate*/    
    for(i=0;i<resident->N_node_families;i++)
    {
        snprintf(filename,sizeof(char)*32,"protein_%i.txt",i);
        fp=fopen(filename,"w");
        for(j=0;j<N_REPLICATES;j++)
        {
            for(k=0;k<phenotype[j].total_time_points;k++)                
                fprintf(fp,"%f ",phenotype[j].protein_concentration[k+i*phenotype[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    } 
    /*gene-specific concentration: each gene has its own file, in which each row is a replicate*/
    for(i=0;i<resident->ngenes;i++)
    {
        snprintf(filename,sizeof(char)*32,"gene_%i.txt",i);
        fp=fopen(filename,"w");
        for(j=0;j<N_REPLICATES;j++)
        {
            for(k=0;k<phenotype[j].total_time_points;k++)                
                fprintf(fp,"%f ",phenotype[j].gene_specific_concentration[k+i*phenotype[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    /*output which gene is which protein*/
    fp=fopen("gene_and_protein.txt","w");
    fprintf(fp,"gene protein is_effector\n");
    for(i=0;i<resident->ngenes;i++)
    {
        if(resident->is_output[i]==OUTPUT_PROTEIN)
            fprintf(fp,"%d %d Y\n",i,resident->which_node_family[i]);
        else
            fprintf(fp,"%d %d N\n",i,resident->which_node_family[i]);
    }
    fclose(fp);   
    /*release memory*/
    for(i=0;i<N_THREADS;i++)
    {
        free(phenotype[i].gene_specific_concentration);    
        free(phenotype[i].protein_concentration);        
    }  
}

static void sample_effector_expression_lvl(Genotype *resident, 
                                            Selection *selection, 
                                            Mutation *mut_record, 
                                            int init_mRNA[MAX_GENES], 
                                            float init_protein[MAX_GENES], 
                                            int replay_N_steps,
                                            RngStream RS_parallel[N_THREADS])
{
    int i, j, k;
    float fitness1[N_REPLICATES],fitness2[N_REPLICATES];  
    Phenotype phenotype[N_REPLICATES];
    float mean_expression_lvl[(int)(selection->env1.t_development/sampling_interval)+1];
    FILE *fp,*fp2;

    /*load mutation record*/
    fp=fopen(mutation_file,"r");    
    if(fp!=NULL)        
        printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
    else
    {
        printf("Loading mutation record failed! Quit program!");
#if MAKE_LOG
        LOG("Loading mutation record failed!");
#endif
        exit(-2);
    } 
    
    /*initialize phenotype*/
    for (i=0;i<N_REPLICATES;i++)
    {
        phenotype[i].total_time_points=(int)(selection->env1.t_development/sampling_interval)+1;       
        phenotype[i].protein_concentration=(float *)malloc(phenotype[i].total_time_points*sizeof(float));    
        phenotype[i].timepoint=0;         
    }    
   
    /*expression in the initial genotype*/
    calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1, fitness2, RM_NONE);    
    for(i=0;i<(int)(selection->env1.t_development/sampling_interval)+1;i++)
    {
        mean_expression_lvl[i]=0.0;
        for(j=0;j<N_REPLICATES;j++)
            mean_expression_lvl[i]+=phenotype[j].protein_concentration[i];
        mean_expression_lvl[i]/=(float)N_REPLICATES;
    }    
    fp2=fopen("effector_expression_lvl.txt","w");
    fprintf(fp2,"0 ");
    for(i=0;i<(int)(selection->env1.t_development/sampling_interval)+1;i++)
        fprintf(fp2,"%.10f ",mean_expression_lvl[i]);
    fprintf(fp2,"\n");
    fclose(fp2);
    /*the evolved genotypes*/
    for(i=1;i<=replay_N_steps;i++)
    {  
        for(j=0;j<N_REPLICATES;j++)
            phenotype[j].timepoint=0;
        fscanf(fp,"%c %d %d %d %s %d %a\n",&(mut_record->mut_type),
                                            &(mut_record->which_gene),                                                    
                                            &(mut_record->which_nucleotide), 
                                            &(mut_record->which_protein),
                                            mut_record->nuc_diff,               
                                            &(mut_record->kinetic_type),
                                            &(mut_record->kinetic_diff));
        reproduce_mutate(resident,mut_record);
        calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1, fitness2, RM_NONE);
        for(j=0;j<(int)(selection->env1.t_development/sampling_interval)+1;j++)
        {
            mean_expression_lvl[j]=0.0;
            for(k=0;k<N_REPLICATES;k++)
                mean_expression_lvl[j]+=phenotype[k].protein_concentration[j];
            mean_expression_lvl[j]/=(float)N_REPLICATES;
        }    
        fp2=fopen("effector_expression_lvl.txt","a+");
        fprintf(fp2,"%d ",i);
        for(j=0;j<(int)(selection->env1.t_development/sampling_interval)+1;j++)
            fprintf(fp2,"%.10f ",mean_expression_lvl[j]);
        fprintf(fp2,"\n");
        fclose(fp2);
    }
    fclose(fp); 
}

static int determine_adaptive_TFBSs(Genotype *resident, 
                                        Selection *selection,
                                        float original_fitness,
                                        float original_fitness_se, 
                                        int init_mRNA[MAX_GENES], 
                                        float init_protein[MAX_GENES], 
                                        RngStream RS_parallel[N_THREADS],
                                        int N_rep,
                                        float fitness[10],
                                        float cutoff_of_nonadaptive_TFBS)
{
    float fitness1[HI_RESOLUTION_RECALC][N_REPLICATES],fitness2[HI_RESOLUTION_RECALC][N_REPLICATES]; 
    Phenotype phenotype[N_REPLICATES];
    int rm_what=RM_NONE;
    int j;
       
    fitness[0]=original_fitness;
    fitness[5]=original_fitness_se;
    
    /*remove 2-mismatch TFBSs of the signal on non-effector transcription factor genes*/
    for(j=0;j<N_rep;j++)
        calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1[j], fitness2[j], RM_S2T);
    calc_fitness_stats( resident, selection, &(fitness1[0]), &(fitness2[0]), N_rep); 
    fitness[1]=resident->avg_fitness;
    fitness[6]=resident->SE_avg_fitness;
    if((resident->avg_fitness-original_fitness)>=(cutoff_of_nonadaptive_TFBS*original_fitness))
        rm_what+=RM_S2T;
    
    /*remove 2-mismatch TFBSs of the effector on non-effector transcription factor genes*/
    for(j=0;j<N_rep;j++)
        calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1[j], fitness2[j], RM_E2T);
    calc_fitness_stats( resident, selection, &(fitness1[0]), &(fitness2[0]),N_rep); 
    fitness[2]=resident->avg_fitness;
    fitness[7]=resident->SE_avg_fitness;
    if((resident->avg_fitness-original_fitness)>=(cutoff_of_nonadaptive_TFBS*original_fitness))
        rm_what+=RM_E2T;
    
    /*remove 2-mismatch TFBSs of the effector on the effector genes*/
    for(j=0;j<N_rep;j++)
        calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1[j], fitness2[j], RM_E2E);
    calc_fitness_stats( resident, selection, &(fitness1[0]), &(fitness2[0]),N_rep); 
    fitness[3]=resident->avg_fitness;
    fitness[8]=resident->SE_avg_fitness;
    if((resident->avg_fitness-original_fitness)>=(cutoff_of_nonadaptive_TFBS*original_fitness))
        rm_what+=RM_E2E;
    
    /*remove 2-mismatch TFBSs of the non-effector TFs on the effector genes*/
    for(j=0;j<N_rep;j++)
        calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1[j], fitness2[j], RM_T2E);
    calc_fitness_stats( resident, selection, &(fitness1[0]), &(fitness2[0]), N_rep); 
    fitness[4]=resident->avg_fitness;
    fitness[9]=resident->SE_avg_fitness;
    if((resident->avg_fitness-original_fitness)>=(cutoff_of_nonadaptive_TFBS*original_fitness))
        rm_what+=RM_T2E;
    
    return rm_what;
}

static void restore_all_TFBS(Genotype *resident)
{
    int i;    
    resident->flag_rm_which_TFBS=RM_NONE;
    for(i=N_SIGNAL_TF;i<resident->ngenes;i++)
        resident->recalc_TFBS[i]=YES;
    calc_all_binding_sites(resident);
}

static void run_simulation( Genotype *resident, 
                            Genotype *mutant, 
                            Mutation *mut_record,
                            Selection *burn_in,
                            Selection *selection, 
                            int init_mRNA[MAX_GENES],  
                            float init_protein[MAX_PROTEINS],
                            int init_N_tot_mutations,   //this is the 
                            int init_step,              //init_step is either 0 or loaded from a saving point
                            RngStream RS_main,
                            RngStream RS_parallel[N_THREADS])
{
    FILE *fp;
    int i;
    int flag_burn_in,N_tot_trials,first_step; 
    Output_buffer resident_info[OUTPUT_INTERVAL];
    Phenotype phenotype[N_REPLICATES];
    first_step=init_step;
    N_tot_trials=init_N_tot_mutations;
    
    /* first, run burn-in */
    if(burn_in->MAX_STEPS!=0)
    {
        flag_burn_in=1; 
        DUPLICATION=burn_in->temporary_DUPLICATION;                 
        SILENCING=burn_in->temporary_SILENCING; 
        miu_ACT_TO_INT_RATE=burn_in->temporary_miu_ACT_TO_INT_RATE; 
        miu_Kd=burn_in->temporary_miu_Kd;       
        miu_protein_syn_rate=burn_in->temporary_miu_protein_syn_rate; 
        float fitness1[HI_RESOLUTION_RECALC][N_REPLICATES],fitness2[HI_RESOLUTION_RECALC][N_REPLICATES];      
        
        if(evolve_N_steps(  resident, 
                            mutant,
                            mut_record, 
                            burn_in,
                            resident_info,
                            &first_step,                     
                            &N_tot_trials, 
                            init_mRNA, 
                            init_protein,
                            RS_main,
                            RS_parallel,
                            flag_burn_in)==-1)   
            
            return;     
        
        /*Calculate fitness of the current genotype under post burn_in condition*/
        /*The "if" is always true when the simulation is run from the beginning,
         *i.e. when init_step=0. But when continuing a simulation from a saving  
         *point after the burn-in, the "if" is always false*/
        if(init_step<burn_in->MAX_STEPS)
        {            
            for(i=0;i<HI_RESOLUTION_RECALC;i++) 
                calc_avg_fitness(   resident, 
                                    selection,
                                    phenotype,
                                    init_mRNA,
                                    init_protein,
                                    RS_parallel,                                        
                                    fitness1[i],
                                    fitness2[i],
                                    RM_NONE); 
            calc_fitness_stats(resident,selection,&(fitness1[0]),&(fitness2[0]),HI_RESOLUTION_RECALC);   
            /*calculate the number of c1-ffls*/
            find_motifs(resident);
            /*save resident status to output buffer*/                  
            store_resident_info(resident, 
                                NULL, 
                                &(resident_info[OUTPUT_INTERVAL-1]), 
                                NA, 
                                NA, 
                                NA,
                                (float)NA,
                                0);      //magic number 0 means to store only fitness and number of motifs
            
            output_resident_info(resident_info, OUTPUT_INTERVAL, 1); //magic number 1 means to output everything
            
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
            /* marks the last step at which all state of the program has been output*/
            fp=fopen("saving_point.txt","w");
            fprintf(fp,"%d %d\n",burn_in->MAX_STEPS,N_tot_trials);
            fclose(fp);
        }
    }    
    
    /* post-burn-in simulations*/
    flag_burn_in=0;    
    DUPLICATION=selection->temporary_DUPLICATION;                 
    SILENCING=selection->temporary_SILENCING; 
    miu_ACT_TO_INT_RATE=selection->temporary_miu_ACT_TO_INT_RATE; 
    miu_Kd=selection->temporary_miu_Kd;       
    miu_protein_syn_rate=selection->temporary_miu_protein_syn_rate; 
    
    if(evolve_N_steps(  resident, 
                        mutant,
                        mut_record, 
                        selection,
                        resident_info,
                        &first_step,                   
                        &N_tot_trials,   
                        init_mRNA,
                        init_protein,
                        RS_main,
                        RS_parallel,
                        flag_burn_in)==-1);     
    return;
}


static void continue_simulation(Genotype *resident, 
                                Genotype *mutant, 
                                Mutation *mut_record, 
                                Selection *burn_in,
                                Selection *selection,
                                int replay_N_steps, 
                                int init_mRNA[MAX_GENES],
                                float init_protein[MAX_PROTEINS],                            
                                RngStream RS_main,
                                RngStream RS_parallel[N_THREADS])
{
    int i,j,N_tot_mutations;    
    unsigned long rng_seeds[N_THREADS+1][6];
    char buffer[400]; 
    FILE *fp;

    /*delete the incomplete lines in the output files*/
    tidy_output_files(evo_summary,mutation_file);
    
    /* set genotype based on previous steps*/   
    replay_mutations(resident, mut_record, replay_N_steps); 

    /* set random number seeds*/
    fp=fopen("RngSeeds.txt","r");
    if(fp!=NULL)
    {
        for(i=0;i<replay_N_steps/OUTPUT_INTERVAL;i++)
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
#if MAKE_LOG
        LOG("cannot open RngSeeds.txt\n");     
#endif
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
            fgets(buffer,400,fp);
        fscanf(fp,"%d %d %a %a %a %a %a %a\n",&N_tot_mutations, 
                                                        &(mut_record->N_hit_bound),
                                                        &(resident->avg_fitness),                                            
                                                        &(resident->fitness1), 
                                                        &(resident->fitness2), 
                                                        &(resident->SE_avg_fitness),
                                                        &(resident->SE_fitness1),
                                                        &(resident->SE_fitness2));
    }
    else
    {
#if MAKE_LOG
        LOG("cannot open precise_fitness.txt\n");       
#endif
        exit(-2);
    }        
    fclose(fp);
    /*continue running simulation*/
    run_simulation( resident, 
                    mutant,
                    mut_record,
                    burn_in,
                    selection,
                    init_mRNA,
                    init_protein,
                    N_tot_mutations,
                    replay_N_steps+1, 
                    RS_main,
                    RS_parallel);   
}

static void calc_fitness_stats(Genotype *genotype,
                                    Selection *selection,
                                    float (*f1)[N_REPLICATES],
                                    float (*f2)[N_REPLICATES],                                  
                                    int N_recalc_fitness)
{
    float avg_f1=0.0;
    float avg_f2=0.0;  
    float sum_sq_diff_f1=0.0;
    float sum_sq_diff_f2=0.0;   
    float sum_sq_diff_mean_f=0.0;
    float diff_f1,diff_f2,sq_SE_f1,sq_SE_f2;  
    int counter=0;
    int i,j;

    for(i=0;i<N_recalc_fitness;i++)
    {
        for(j=0;j<N_REPLICATES;j++)
        {
            avg_f1+=f1[i][j];
            avg_f2+=f2[i][j];  
            genotype->fitness_measurement[counter]=selection->env1_weight*f1[i][j]+selection->env2_weight*f2[i][j];
            counter++;
        }
    }
    avg_f1=avg_f1/(N_recalc_fitness*N_REPLICATES);
    avg_f2=avg_f2/(N_recalc_fitness*N_REPLICATES);  
    
    for(i=0;i<N_recalc_fitness;i++)
    {
        for(j=0;j<N_REPLICATES;j++)
        {
            diff_f1=f1[i][j]-avg_f1;
            diff_f2=f2[i][j]-avg_f2;
            sum_sq_diff_f1+=pow(diff_f1,2.0);
            sum_sq_diff_f2+=pow(diff_f2,2.0);
            sum_sq_diff_mean_f+=pow(diff_f1*selection->env1_weight+diff_f2*selection->env2_weight,2.0);
        }
    }
    sq_SE_f1=sum_sq_diff_f1/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));
    sq_SE_f2=sum_sq_diff_f2/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));  
    genotype->fitness1=avg_f1;
    genotype->fitness2=avg_f2;
    genotype->SE_fitness1=sqrt(sq_SE_f1);
    genotype->SE_fitness2=sqrt(sq_SE_f2); 
    genotype->avg_fitness=selection->env1_weight*avg_f1+selection->env2_weight*avg_f2;
    genotype->SE_avg_fitness=sqrt(sum_sq_diff_mean_f/(N_recalc_fitness*N_REPLICATES-1)/(N_recalc_fitness*N_REPLICATES));  
}

static int evolve_N_steps(  Genotype *resident, 
                            Genotype *mutant,
                            Mutation *mut_record, 
                            Selection *selection,
                            Output_buffer resident_info[OUTPUT_INTERVAL],
                            int *init_step,                       
                            int *N_tot_trials,        
                            int init_mRNA[MAX_GENES],   
                            float init_protein[MAX_PROTEINS],
                            RngStream RS_main,
                            RngStream RS_parallel[N_THREADS],
                            int flag_burn_in)
{
    int i,j;
    int output_counter=0;
    int flag_replaced; 
    int N_trials;
    float fitness1[HI_RESOLUTION_RECALC][N_REPLICATES],fitness2[HI_RESOLUTION_RECALC][N_REPLICATES];
    Phenotype phenotype[N_REPLICATES];
    float selection_coefficient; 
    FILE *fp;
#if OUTPUT_MUTANT_DETAILS
    Output_buffer *mutant_info;  
    int mutant_counter=0;
    int current_mutant_info_size=OUTPUT_INTERVAL*50;
    mutant_info=(Output_buffer *)malloc(current_mutant_info_size*sizeof(Output_buffer));
#endif
 
    for(i=(*init_step);i<=selection->MAX_STEPS;i++)
    {             
        flag_replaced=0;      
        N_trials=0;
        
        /*try mutations until one replaces the current genotype*/
        while(!flag_replaced) 
        {	
            N_trials++;
            (*N_tot_trials)++;
            if(N_trials>MAX_TRIALS) /*Tried too many mutation in one step.*/
            {
#if OUTPUT_MUTANT_DETAILS                
                output_mutant_info(mutant_info,mutant_counter);
#endif
                output_resident_info(resident_info,output_counter,1);
                fp=fopen(evo_summary,"a+");                              
                fprintf(fp,"Tried %d mutations, yet none could fix\n",MAX_TRIALS);               
                fclose(fp); 
                summarize_binding_sites(resident,i-1);
                return -1;
            }

            /*do mutation on a copy of the current genotype*/
            clone_genotype(resident,mutant);
            mutate(mutant,RS_main,mut_record);            
            
            /*determine if we need more space to store TFBSs*/
            calc_all_binding_sites(mutant);           
            MAX_TFBS_NUMBER=mutant->N_allocated_elements;            
         
            /*calculate the fitness of the mutant at low resolution*/
            calc_avg_fitness(mutant, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1[0], fitness2[0], RM_NONE);
            calc_fitness_stats(mutant, selection, &(fitness1[0]), &(fitness2[0]), 1); // calc fitness at low resolution

#if OUTPUT_MUTANT_DETAILS
            if(mutant_counter>=current_mutant_info_size)
            {
                current_mutant_info_size+=OUTPUT_INTERVAL*100;
                mutant_info=(Output_buffer *)realloc(mutant_info,current_mutant_info_size*sizeof(Output_buffer));
            }
            store_mutant_info(mutant,mut_record,&(mutant_info[mutant_counter]),i,*N_tot_trials);      
            mutant_counter++;
#endif
            /*Can the mutant replace the current genotype?*/
            try_replacement(resident, mutant, &flag_replaced, &selection_coefficient);
        }
        
        /*replace the current genotype by overwriting it*/
        clone_genotype(mutant,resident);        
        calc_all_binding_sites(resident); 
     
        /*calculate mutant fitness at high resolution*/ 
        /*If we are at the last step of BURN_IN, 
         * we should calculate the fitness under the post-burn-in condition, 
         * which is done in run_simulation, outside the current function*/
        if(!(i==selection->MAX_STEPS && flag_burn_in)) 
        {
            for(j=1;j<HI_RESOLUTION_RECALC;j++)  
                calc_avg_fitness(resident, selection, phenotype, init_mRNA, init_protein, RS_parallel, fitness1[j], fitness2[j], RM_NONE);              
            calc_fitness_stats(resident, selection, &(fitness1[0]), &(fitness2[0]), HI_RESOLUTION_RECALC);   
        }  
        
        /*calculate the number of motifs*/
        find_motifs(resident);
        /*save resident status to output buffer*/ 
        store_resident_info(resident, mut_record, &(resident_info[output_counter]), i, N_trials, *N_tot_trials, selection_coefficient,1); //magic number 1 means store everything
        output_counter++;       
        
        /*output network topology every OUTPUT_INTERVAL steps*/
        if(i%OUTPUT_INTERVAL==0 && i!=0)
        {  
            summarize_binding_sites(resident,i);
#if OUTPUT_MUTANT_DETAILS 
            output_mutant_info(mutant_info,mutant_counter);
            mutant_counter=0;
#endif
            if(!(i==selection->MAX_STEPS && flag_burn_in))
            {
                output_resident_info(resident_info,output_counter,1); //magic number 1 means to output everything
                output_counter=0;
            }
        }
           
        /* output rng seeds*/
#if OUTPUT_RNG_SEEDS
        unsigned long seeds[6];
        if(!(i==selection->MAX_STEPS && flag_burn_in) && i%OUTPUT_INTERVAL==0)
        {
            RngStream_GetState(RS_main,seeds);
            fp=fopen("RngSeeds.txt","a+");
            fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]); 
            for(j=0;j<N_THREADS;j++)
            {
                RngStream_GetState(RS_parallel[j],seeds);
                fprintf(fp,"%lu %lu %lu %lu %lu %lu ",seeds[0],seeds[1],seeds[2],seeds[3],seeds[4],seeds[5]);                
            }
            fprintf(fp,"\n");
            fflush(fp);
            fclose(fp);
        }
#endif       
        /*output precise hi-resolution fitness*/
        if(!(i==selection->MAX_STEPS && flag_burn_in))
        {
            if(i%OUTPUT_INTERVAL==0)
            {
                fp=fopen("saving_point.txt","w");
                fprintf(fp,"%d %d\n",i,*N_tot_trials);
                fflush(fp);
                fclose(fp);
            }    
        }
    } 
    *init_step=i;
#if OUTPUT_MUTANT_DETAILS
    free(mutant_info);
#endif
    return 0;
}

static void print_motifs(Genotype *genotype)
{
    FILE *fp; 
    int i;
    fp=fopen("N_motifs.txt","a+");
    for(i=0;i<33;i++)
        fprintf(fp,"%d ",genotype->N_motifs[i]);    
    fprintf(fp,"\n");
    fclose(fp); 
#if COUNT_NEAR_AND
    fp=fopen("N_near_AND_gate_motifs.txt","a+");
    for(i=0;i<12;i++)    
        fprintf(fp,"%d ",genotype->N_near_AND_gated_motifs[i]);
    fprintf(fp,"\n");    
    fclose(fp);    
#endif
}

static void summarize_binding_sites(Genotype *genotype,int step_i)
{
    FILE *OUTPUT1;
    int i,j,k,which_target_gene,which_TF_gene,cut_off, which_node, which_protein;
    int table[MAX_GENES][MAX_GENES],table_cp[MAX_GENES][MAX_GENES];  
    
    /* store the number of TFBSs in a matrix.*/
    for(i=0;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->ngenes;j++)
        {
            table[i][j]=0;  
            table_cp[i][j]=0;  
        }
    }
   
    /*loop through cis-reg clusters. Genes with the same cis-reg sequence must
     *have the same TFBSs
     */
    for(i=N_SIGNAL_TF;i<genotype->N_cisreg_clusters;i++)
    {
        which_target_gene=genotype->cisreg_cluster_pool[i][1][0];
       
        for(j=0;j<genotype->N_node_families;j++)    
        {           
            which_TF_gene=genotype->node_family_pool[j][1][0];            
            /*exclude weak TFBS?*/
            if(genotype->is_output[which_TF_gene]==OUTPUT_PROTEIN) 
            {
                cut_off=MAX_MISMATCH_EFFECTOR2GENES;  
            }
            else if(which_TF_gene==N_SIGNAL_TF-1)
            {
                cut_off=MAX_MISMATCH_SIGNAL2GENES;
            }
            else
            {
                cut_off=MAX_MISMATCH_NONEFFECTOR2GENES;
            }
            
            /*count number of TFBSs*/
            for(k=0;k<genotype->binding_sites_num[which_target_gene];k++)
            { 
                which_node=genotype->which_node_family[genotype->protein_pool[genotype->all_binding_sites[which_target_gene][k].tf_id][1][0]];
                if(which_node==j && genotype->all_binding_sites[which_target_gene][k].mis_match<=cut_off)
                    table_cp[which_target_gene][j]++; 
            }             
        }
        
        /*Genes with the same cis-reg sequence must have the same TFBSs*/
        for(j=0;j<genotype->cisreg_cluster_pool[i][0][0];j++)
        { 
            for(k=0;k<genotype->nproteins;k++)
                table[genotype->cisreg_cluster_pool[i][1][j]][k]=table_cp[which_target_gene][k];
        }
    }   
    
    /*Output all binding sites*/ 
    OUTPUT1=fopen("networks.txt","a+");
    fprintf(OUTPUT1,"step %d\n",step_i);
    fprintf(OUTPUT1,"Gene   ");    
    for(i=0;i<genotype->N_node_families;i++)
    {
        which_TF_gene=genotype->node_family_pool[i][1][0];  
        which_protein=genotype->which_protein[which_TF_gene];
        if(i<10)
        {
            if(genotype->protein_identity[which_protein]==ACTIVATOR)
            {
                if(genotype->is_output[which_TF_gene]==OUTPUT_PROTEIN)
                    fprintf(OUTPUT1,"A%d     ",i);
                else
                    fprintf(OUTPUT1,"a%d     ",i);
            }
            if(genotype->protein_identity[which_protein]==REPRESSOR)
            {
                if(genotype->is_output[which_TF_gene]==OUTPUT_PROTEIN)
                    fprintf(OUTPUT1,"R%d     ",i);
                else
                    fprintf(OUTPUT1,"r%d     ",i);
            }
        }
        else
        {
            if(genotype->protein_identity[which_protein]==ACTIVATOR)
            {
                if(genotype->is_output[which_TF_gene]==OUTPUT_PROTEIN)
                    fprintf(OUTPUT1,"A%d    ",i);
                else
                    fprintf(OUTPUT1,"a%d    ",i);
            }
            if(genotype->protein_identity[which_protein]==REPRESSOR)
            {
                if(genotype->is_output[which_TF_gene]==OUTPUT_PROTEIN)
                    fprintf(OUTPUT1,"R%d    ",i);
                else
                    fprintf(OUTPUT1,"r%d    ",i);
            }
        }
    }
    fprintf(OUTPUT1,"which_protein ");
    fprintf(OUTPUT1,"AND_gate_capable\n");    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        if(i<10)
            fprintf(OUTPUT1,"%d     ",i);
        else
            fprintf(OUTPUT1,"%d    ",i);
        
        for(j=0;j<genotype->N_node_families;j++)
        {
            if(table[i][j]<10)
                fprintf(OUTPUT1," %d     ",table[i][j]);
            else
                fprintf(OUTPUT1," %d    ",table[i][j]);
        }
        if(genotype->is_output[i]==OUTPUT_PROTEIN)
        {
            if(genotype->protein_identity[genotype->which_protein[i]]==ACTIVATOR)
                fprintf(OUTPUT1,"      a%d",genotype->which_node_family[i]); 
            if(genotype->protein_identity[genotype->which_protein[i]]==REPRESSOR)
                fprintf(OUTPUT1,"      r%d",genotype->which_node_family[i]);
        }
        else
        {
            if(genotype->protein_identity[genotype->which_protein[i]]==ACTIVATOR)
                fprintf(OUTPUT1,"      a%d",genotype->which_node_family[i]); 
            if(genotype->protein_identity[genotype->which_protein[i]]==REPRESSOR)
                fprintf(OUTPUT1,"      r%d",genotype->which_node_family[i]);
        }
        if(genotype->min_N_activator_to_transc[i]==1)
            fprintf(OUTPUT1,"             N\n");
        else
            fprintf(OUTPUT1,"             Y\n");
    }
    fprintf(OUTPUT1,"\n"); 
    fprintf(OUTPUT1,"\n"); 
    fclose(OUTPUT1);
  
}

static void find_motifs(Genotype *genotype)
{
    int i,j,k;  
    int N_activators, N_repressors;
    int effector_gene, effector_node, effector_protein;
    int auxiliary_tf_gene, auxiliary_tf_node, auxiliary_tf_protein;
    int repressors[MAX_PROTEINS]; //which nodes repress an effector gene
    int activators[MAX_PROTEINS];
    int regulated_by_signal[MAX_GENES];   
    int effector_self_rep,found_i1_or_NFBL;
    
    /*reset variables*/       
    for(i=0;i<MAX_GENES;i++)
    {
        genotype->gene_in_core_C1ffl[i]=0;
        for(j=0;j<MAX_PROTEINS;j++)
            genotype->TF_in_core_C1ffl[i][j]=0;
    }     
    for(i=0;i<33;i++)    
        genotype->N_motifs[i]=0;  
    
    /*find genes regulated by the signal*/    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        if(find_TFBS_of_A_on_B(genotype,N_SIGNAL_TF-1,i,MAX_MISMATCH_SIGNAL2GENES))
            regulated_by_signal[i]=YES;
        else
            regulated_by_signal[i]=NO;  
    }    
    
    /*find motifs that contain an output gene*/
    for(i=0;i<genotype->n_output_genes;i++)
    {       
        effector_gene=genotype->output_gene_ids[i];
        
        /*find out what TF node families regulate the effector genes */
        who_regulates_effector(genotype,effector_gene,activators,repressors,&N_activators,&N_repressors); 

        /***count i1ffls and NFBLs ***/  
        if(activators[0]==N_SIGNAL_TF-1) //effector is activated by the signal
        {
            effector_node=genotype->which_node_family[effector_gene]; 
            effector_protein=genotype->which_protein[effector_gene];
            
            /*Does the effector repress itself?*/
            effector_self_rep=0;               
            for(j=0;j<N_repressors;j++)
            {
                auxiliary_tf_node=repressors[j];
                if(effector_node==auxiliary_tf_node)
                {
                    effector_self_rep=1;
                    break;
                }
            }               
            
            /*set flag, which is used to determine there is stand-alone auto-repressor*/            
            found_i1_or_NFBL=0;
            
            /*check genes that repress the effector for potential motifs*/
            for(j=0;j<N_repressors;j++) 
            {
                auxiliary_tf_node=repressors[j];
                if(effector_node!=auxiliary_tf_node) // if the repressor is not a variant of the effector gene
                {
                    if(effector_self_rep==0) // if the effector does not repress itself
                    {
                        for(k=0;k<genotype->node_family_pool[auxiliary_tf_node][0][0];k++) // check all variants of the repressor
                        {
                            auxiliary_tf_gene=genotype->node_family_pool[auxiliary_tf_node][1][k]; 
                            auxiliary_tf_protein=genotype->which_protein[auxiliary_tf_gene];

                            if(regulated_by_signal[auxiliary_tf_gene]==YES)  // aux. is regulated by signal
                            {
                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][N_SIGNAL_TF-1]==ACTIVATOR) //the signal can activate aux. tf
                                {    
                                    found_i1_or_NFBL=1;
                                    if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene,MAX_MISMATCH_EFFECTOR2GENES)) //aux. tf is regulated by effector
                                    {
                                        if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) // the effector is an activator to aux. gene
                                        {                                            
                                            if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?                                           
                                                genotype->N_motifs[1]++; //I1_with_2_effector+NFBF_with_2_effector 
                                            else
                                                genotype->N_motifs[2]++; //I1+NFBL
                                        }
                                        else //overlapping I1  
                                        {
                                            if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?
                                                genotype->N_motifs[7]++;   //overlapping I1w2   
                                            else
                                                genotype->N_motifs[8]++;    //overlapping I1                       
                                        }
                                    }
                                    else // aux. tf is not regulated by effector, producing i1
                                    {
                                        if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?                                                
                                            genotype->N_motifs[3]++; //I1w2 
                                        else
                                            genotype->N_motifs[4]++; //I1
                                    }
                                } 
                            }
                            else // aux. tf is not regulated by the signal
                            {  
                                if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene,MAX_MISMATCH_EFFECTOR2GENES) &&
                                    genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) // if effector binds to and activates aux.
                                {
                                    found_i1_or_NFBL=1;
                                    if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?                                              
                                        genotype->N_motifs[5]++; //NFBFw2 
                                    else
                                        genotype->N_motifs[6]++; //NFBL     
                                }                           
                            }
                        }  
                    }
                    else //effector represses itself
                    {
                        for(k=0;k<genotype->node_family_pool[auxiliary_tf_node][0][0];k++) // check all variants of the repressor
                        {
                            auxiliary_tf_gene=genotype->node_family_pool[auxiliary_tf_node][1][k]; 
                            auxiliary_tf_protein=genotype->which_protein[auxiliary_tf_gene];

                            if(regulated_by_signal[auxiliary_tf_gene]==YES)  // aux. is regulated by signal
                            {
                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][N_SIGNAL_TF-1]==ACTIVATOR) //the signal can activate aux. tf
                                {           
                                    found_i1_or_NFBL=1;
                                    if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene,MAX_MISMATCH_EFFECTOR2GENES)) //aux. tf is regulated by effector
                                    {
                                        if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) // the effector is an activator to aux. gene
                                        {
                                            if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?                                           
                                                genotype->N_motifs[9]++; //I1w2+NFBFw2+auto_rep 
                                            else
                                                genotype->N_motifs[10]++; //I1+NFBL+auto_rep 
                                        }
                                        else //overlapping I1  
                                        {
                                            if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?
                                                genotype->N_motifs[15]++; //overlapping I1w2+auto_rep   
                                            else
                                                genotype->N_motifs[16]++;  //overlapping I1+auto_rep                            
                                        }
                                    }
                                    else // aux. tf is not regulated by effector, producing i1
                                    {
                                        if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?                                                
                                            genotype->N_motifs[11]++; //I1w2+auto_rep  
                                        else
                                            genotype->N_motifs[12]++; //I1+auto_rep 
                                    }
                                } 
                            }
                            else // aux. tf is not regulated by the signal
                            {  
                                if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene,MAX_MISMATCH_EFFECTOR2GENES) &&
                                    genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) // if effector binds to and activates aux.
                                {
                                    found_i1_or_NFBL=1;
                                    if(genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN) //auxiliary tf is effector too?                                              
                                        genotype->N_motifs[13]++; //NFBFw2+auto_rep 
                                    else
                                        genotype->N_motifs[14]++; //NFBL+auto_rep      
                                }                           
                            }
                        }  
                    }
                }   
                else //the repressor is a variant of the effector gene
                {
                    if(genotype->protein_pool[effector_protein][0][0]>1) // effector gene has more than one copy
                    {
                        for(k=0;k<genotype->node_family_pool[auxiliary_tf_node][0][0];k++)
                        {
                            auxiliary_tf_gene=genotype->node_family_pool[auxiliary_tf_node][1][k]; 
                            
                            if(effector_gene==auxiliary_tf_gene) // do not count the effector gene itself
                                continue;
                            
                            auxiliary_tf_protein=genotype->which_protein[auxiliary_tf_gene];

                            if(regulated_by_signal[auxiliary_tf_gene]==YES)  // aux. is regulated by signal
                            {
                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][N_SIGNAL_TF-1]==ACTIVATOR) //the signal can activate aux. tf
                                {                                
                                    if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene,MAX_MISMATCH_EFFECTOR2GENES)) //aux. tf is regulated by effector
                                    {
                                        if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) // the effector is an activator to aux. gene
                                        {                                       
                                            genotype->N_motifs[9]++; //I1w2+NFBFw2+auto_rep  
                                        }
                                        else //overlapping I1w1  
                                        {
                                            genotype->N_motifs[15]++; //overlapping I1w2+auto_rep   
                                        }
                                    }
                                    else // aux. tf is not regulated by effector, producing i1
                                    {                               
                                        genotype->N_motifs[11]++; //I1w2+auto_rep  
                                    }
                                } 
                            }
                            else // aux. tf is not regulated by the signal
                            {  
                                if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene,MAX_MISMATCH_EFFECTOR2GENES) &&
                                    genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) // if effector binds to and activates aux.
                                {                                                                         
                                    genotype->N_motifs[13]++; //NFBFw2+auto_rep  
                                }                           
                            }
                        }  
                    }
                    else // the variant is the effector gene itself
                    {
                        if(found_i1_or_NFBL==0)
                            genotype->N_motifs[0]++; //stand-alone auto-rep
                    }
                }
            }         
        }
    }     
}

static int find_TFBS_of_A_on_B(Genotype *genotype, int gene_A, int gene_B, int cut_off_mismatch)
{
    int site_id,protein_id;
    int found_bs=0;
    protein_id=genotype->which_protein[gene_A];
    for(site_id=0;site_id<genotype->binding_sites_num[gene_B];site_id++)
    {
        if(genotype->all_binding_sites[gene_B][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_B][site_id].mis_match<=cut_off_mismatch)
        {
            found_bs=1;
            break;
        }
    }
    return found_bs;
}

static void who_regulates_effector( Genotype *genotype, 
                                    int effector_gene, 
                                    int activators[MAX_PROTEINS], 
                                    int repressors[MAX_PROTEINS], 
                                    int *N_activators, 
                                    int *N_repressors)
{
    int i, j, protein_id, n_nodes, empty_slot;
    int protein_to_node[MAX_PROTEINS][2][MAX_PROTEINS];
    /*reset tables*/
    for(i=0;i<MAX_PROTEINS;i++)
    {
        repressors[i]=NA;    
        activators[i]=NA;  
    }  
    /*build protein_to_node_family dictionary*/
    for(i=0;i<genotype->nproteins;i++)
    {  
        /*reset*/
        for(j=0;j<genotype->N_node_families;j++)
            protein_to_node[i][1][j]=0; 
        /*which node families does the protein contain*/
        for(j=0;j<genotype->protein_pool[i][0][0];j++) 
            protein_to_node[i][1][genotype->which_node_family[genotype->protein_pool[i][1][j]]]=1;
        /*move present nodes to the front of protein_to_node, and count present nodes*/
        empty_slot=0;
        n_nodes=0; 
        for(j=0;j<genotype->N_node_families;j++)
        {    
            n_nodes+=protein_to_node[i][1][j];
            if(protein_to_node[i][1][j]==1)
            {
                protein_to_node[i][1][empty_slot]=j;
                empty_slot++;
            }
        }
        protein_to_node[i][0][0]=n_nodes;
    }
    
    /*scan binding sites for tfs that regulate effector gene*/
    for(i=0;i<genotype->binding_sites_num[effector_gene];i++)
    {
        protein_id=genotype->all_binding_sites[effector_gene][i].tf_id;
        if(protein_id==N_SIGNAL_TF-1)
        {
            if(genotype->all_binding_sites[effector_gene][i].mis_match<=MAX_MISMATCH_SIGNAL2GENES)
            {
                for(j=0;j<protein_to_node[protein_id][0][0];j++)
                {
                    /* we had planned to allow a TF to be active some genes and repress the others, i.e.
                     * the allowing TF to have locus-specific behavior. The codes below are an example 
                     * of this setting. However, we did not enable such the locus-specific behavior 
                     * during simulations. 
                     */
                    if(genotype->locus_specific_TF_behavior[effector_gene][protein_id]==ACTIVATOR) 
                        activators[protein_to_node[protein_id][1][j]]=YES;
                    else
                        repressors[protein_to_node[protein_id][1][j]]=YES;
                } 
            }
        }
        else if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==NON_OUTPUT_PROTEIN)
        {
            if(genotype->all_binding_sites[effector_gene][i].mis_match<=MAX_MISMATCH_NONEFFECTOR2GENES)
            {
                for(j=0;j<protein_to_node[protein_id][0][0];j++)
                {
                    if(genotype->locus_specific_TF_behavior[effector_gene][protein_id]==ACTIVATOR)
                        activators[protein_to_node[protein_id][1][j]]=YES;
                    else
                        repressors[protein_to_node[protein_id][1][j]]=YES;
                } 
            }
        }
        else    
        {
            if(genotype->all_binding_sites[effector_gene][i].mis_match<=MAX_MISMATCH_EFFECTOR2GENES)
            {
                for(j=0;j<protein_to_node[protein_id][0][0];j++)
                {
                    if(genotype->locus_specific_TF_behavior[effector_gene][protein_id]==ACTIVATOR)
                        activators[protein_to_node[protein_id][1][j]]=YES;
                    else
                        repressors[protein_to_node[protein_id][1][j]]=YES;
                } 
            }
        }
    }
    /* move non-zeros entries in activators and repressors to the front. */
    i=0;    
    *N_activators=0;
    *N_repressors=0;
    for(i=0;i<genotype->N_node_families;i++)
    {
        if(activators[i]==1)               
        {
            activators[*N_activators]=i;
            (*N_activators)++;
        } 
        if(repressors[i]==1)
        {
            repressors[*N_repressors]=i;
            (*N_repressors)++;
        }
    } 
}

static void tidy_output_files(char *file_genotype_summary, char *file_mutations)
{
    int i,replay_N_steps,N_tot_mutations;
    char buffer[2000];    
    FILE *fp1,*fp2;
    
    fp1=fopen("saving_point.txt","r");
    fscanf(fp1,"%d %d",&replay_N_steps,&N_tot_mutations);
    fclose(fp1);
    
    /*Basically delete the last line of the file if it is not complete*/
    fp1=fopen("RngSeeds.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<(replay_N_steps/OUTPUT_INTERVAL);i++)
    {
        fgets(buffer,2000,fp1);
        fputs(buffer,fp2);
        fflush(fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("RngSeeds.txt");
    rename("temp","RngSeeds.txt");
    
    fp1=fopen(file_genotype_summary,"r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps+2;i++)
    {
        fgets(buffer,2000,fp1);
        fputs(buffer,fp2);
        fflush(fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove(file_genotype_summary);
    rename("temp",file_genotype_summary);
    
    fp1=fopen(file_mutations,"r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,2000,fp1);
        fputs(buffer,fp2);
        fflush(fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove(file_mutations);
    rename("temp",file_mutations);
    
#if OUTPUT_MUTANT_DETAILS 
    fp1=fopen("all_mutations.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<N_tot_mutations;i++)
    {
        fgets(buffer,2000,fp1);
        fputs(buffer,fp2);
        fflush(fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("all_mutations.txt");
    rename("temp","all_mutations.txt");
    
    fp1=fopen("fitness_all_mutants.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<N_tot_mutations;i++)
    {
        fgets(buffer,2000,fp1);
        fputs(buffer,fp2);
        fflush(fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("fitness_all_mutants.txt");
    rename("temp","fitness_all_mutants.txt");
#endif
    fp1=fopen("precise_fitness.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,2000,fp1);
        fputs(buffer,fp2);
        fflush(fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("precise_fitness.txt");
    rename("temp","precise_fitness.txt");
}

void print_mutatable_parameters(Genotype *genotype,int init_or_end)
{
    int i;
    FILE *fp;    
    if(init_or_end==1)
        fp=fopen("end_mutatable_parameters.txt","w");
    else
        fp=fopen("init_mutatable_parameters.txt","w");
    
    for(i=0;i<genotype->ngenes;i++)
    {
        fprintf(fp,"%f %f %f %f %d ",genotype->active_to_intermediate_rate[i],
                                    genotype->mRNA_decay_rate[i],
                                    genotype->protein_syn_rate[i],
                                    genotype->protein_decay_rate[i],
                                    genotype->locus_length[i]);        
        fprintf(fp,"%f\n",log10(genotype->Kd[genotype->which_protein[i]]));        
    }
    fclose(fp);
}

static void store_resident_info(Genotype *resident, 
                                Mutation *mut_record, 
                                Output_buffer *resident_info, 
                                int evo_step, 
                                int N_mutations_at_current_step, 
                                int N_tot_mutations, 
                                float selection_coefficient, 
                                int flag)
{  
    int i;
 
   /*always store network motifs*/
    for(i=0;i<33;i++)
        resident_info->n_motifs[i]=resident->N_motifs[i];
//    for(i=0;i<12;i++)        
//        resident_info->n_near_AND_gated_motifs[i]=resident->N_near_AND_gated_motifs[i];
    
    /*if not called by replay_mutation, store other info*/
    if(flag!=-1) 
    {    
        resident_info->avg_f=resident->avg_fitness;
        resident_info->f1=resident->fitness1;
        resident_info->f2=resident->fitness2;
        resident_info->se_avg_f=resident->SE_avg_fitness;
        resident_info->se_f1=resident->SE_fitness1;
        resident_info->se_f2=resident->SE_fitness2;
    
        if(flag==1) //if stores everthing
        {
            resident_info->step=evo_step;
            resident_info->n_mut_at_the_step=N_mutations_at_current_step;
            resident_info->n_tot_mut=N_tot_mutations;
            resident_info->n_hit_bound=mut_record->N_hit_bound;
            resident_info->selection_coefficient=selection_coefficient;

            resident_info->n_gene=resident->ngenes;
            resident_info->n_output_genes=resident->n_output_genes;
            resident_info->n_node_families=resident->N_node_families;
            resident_info->n_act=resident->N_act;
            resident_info->n_rep=resident->N_rep;
         
            resident_info->mut_type=mut_record->mut_type;
            resident_info->which_gene=mut_record->which_gene;
            resident_info->which_protein=mut_record->which_protein;
            resident_info->which_nuc=mut_record->which_nucleotide;
            resident_info->which_kinetic=mut_record->kinetic_type;            
            resident_info->new_nuc[0]=mut_record->nuc_diff[0];
            resident_info->new_nuc[1]=mut_record->nuc_diff[1];
            resident_info->new_nuc[2]=mut_record->nuc_diff[2];
            resident_info->new_kinetic=mut_record->kinetic_diff;
        }
    }
}

static void store_mutant_info(Genotype *mutant, Mutation *mut_record, Output_buffer *mutant_info, int step, int N_tot_mutations)
{
    mutant_info->avg_f=mutant->avg_fitness;
    mutant_info->f1=mutant->fitness1;
    mutant_info->f2=mutant->fitness2;
    mutant_info->se_avg_f=mutant->SE_avg_fitness;
    mutant_info->se_f1=mutant->SE_fitness1;
    mutant_info->se_f2=mutant->SE_fitness2;
    mutant_info->step=step;
    mutant_info->n_tot_mut=N_tot_mutations;
    mutant_info->mut_type=mut_record->mut_type;
    mutant_info->which_gene=mut_record->which_gene;
    mutant_info->which_protein=mut_record->which_protein;
    mutant_info->which_nuc=mut_record->which_nucleotide;
    mutant_info->which_kinetic=mut_record->kinetic_type;
    mutant_info->new_nuc[0]=mut_record->nuc_diff[0];
    mutant_info->new_nuc[1]=mut_record->nuc_diff[1];
    mutant_info->new_nuc[2]=mut_record->nuc_diff[2];
    mutant_info->new_kinetic=mut_record->kinetic_diff;
}

static void output_mutant_info(Output_buffer *mutant_info, int N_mutant)
{
    int i;
    FILE *fp;    
    /*output mutation*/
    fp=fopen("all_mutations.txt","a+");
    for(i=0;i<N_mutant;i++)        
        fprintf(fp,"%d %d %c %d %d %d '%s' %d %a\n",
                mutant_info[i].step,
                mutant_info[i].n_tot_mut,
                mutant_info[i].mut_type,
                mutant_info[i].which_gene,                
                mutant_info[i].which_nuc,
                mutant_info[i].which_protein,
                mutant_info[i].new_nuc,
                mutant_info[i].which_kinetic,
                mutant_info[i].new_kinetic);
    fflush(fp);
    fclose(fp);
    
    /*output mutant fitness, which is low-resolution*/  
    fp=fopen("fitness_all_mutants.txt","a+");
    for(i=0;i<N_mutant;i++) 
        fprintf(fp,"%.10f %.10f %.10f %.10f %.10f %.10f\n", 
            mutant_info[i].avg_f, 
            mutant_info[i].f1,
            mutant_info[i].f2,
            mutant_info[i].se_avg_f,
            mutant_info[i].se_f1,
            mutant_info[i].se_f2);
    fflush(fp);
    fclose(fp); 
}

static void output_resident_info(Output_buffer resident_info[OUTPUT_INTERVAL], int output_counter, int flag)
{
    int i,j;
    FILE *fp;   
    
    /*output motifs*/
    fp=fopen("N_motifs.txt","a+");
    for(j=0;j<output_counter;j++)
    {
        for(i=0;i<33;i++)
            fprintf(fp,"%d ",resident_info[j].n_motifs[i]);    
        fprintf(fp,"\n");        
    }
    fflush(fp);
    fclose(fp); 
//#if COUNT_NEAR_AND
//    fp=fopen("N_near_AND_gated_motifs.txt","a+");
//    for(j=0;j<output_counter;j++)
//    {
//        for(i=0;i<12;i++)    
//            fprintf(fp,"%d ",resident_info[j].n_near_AND_gated_motifs[i]);
//        fprintf(fp,"\n");   
//    }
//    fflush(fp);
//    fclose(fp);    
//#endif
    if(flag==1) //if function is not called by replay_mutation
    {
                        
        /*output mutation info*/
        fp=fopen(mutation_file,"a+");
        for(i=0;i<output_counter;i++)
            fprintf(fp,"%c %d %d %d '%s' %d %a\n", resident_info[i].mut_type,    
                                                resident_info[i].which_gene,                                                
                                                resident_info[i].which_nuc,
                                                resident_info[i].which_protein,
                                                resident_info[i].new_nuc,
                                                resident_info[i].which_kinetic,
                                                resident_info[i].new_kinetic);
        fflush(fp);
        fclose(fp);

        /*output precise fitness */
        fp=fopen("precise_fitness.txt","a+"); 
        for(i=0;i<output_counter;i++)
            fprintf(fp,"%d %d %a %a %a %a %a %a\n",resident_info[i].n_tot_mut, 
                                            resident_info[i].n_hit_bound,
                                            resident_info[i].avg_f,                                                
                                            resident_info[i].f1,  
                                            resident_info[i].f2, 
                                            resident_info[i].se_avg_f,
                                            resident_info[i].se_f1,
                                            resident_info[i].se_f2); 
        fflush(fp);
        fclose(fp);  

        /*output a summary*/
        fp=fopen(evo_summary,"a+");
        for(i=0;i<output_counter;i++)
            fprintf(fp,"%d %d %d %d %c %f %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d\n",
                    resident_info[i].step, 
                    resident_info[i].n_tot_mut, 
                    resident_info[i].n_mut_at_the_step,
                    resident_info[i].n_hit_bound,
                    resident_info[i].mut_type,
                    resident_info[i].selection_coefficient,
                    resident_info[i].avg_f,
                    resident_info[i].f1,
                    resident_info[i].f2,
                    resident_info[i].se_avg_f,  
                    resident_info[i].se_f1,
                    resident_info[i].se_f2,
                    resident_info[i].n_gene,
                    resident_info[i].n_output_genes,                  
                    resident_info[i].n_node_families);
        fflush(fp);
        fclose(fp);
    }
}

/*use the value of flag_rm_which_TFBS to determine which 2-mismatch TFBSs are ignored during perturbation analysis
 *and determine adaptive network motifs
 */
static void ignore_TFBS(Genotype *genotype, int gene_id, int protein_id, int start_TF, int *N_min_match)
{   
    if(genotype->flag_rm_which_TFBS==RM_E2T && genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN)
    {
        if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1; // this is an effector TFBS in TF, so we require it to match at least NMIN+1 nucleotides
        else
            *N_min_match=NMIN;                    
    }
    else if(genotype->flag_rm_which_TFBS==RM_S2T && genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN)
    {
        if(protein_id==start_TF)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_ES2T && genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN)
    {
        if(protein_id==start_TF)
            *N_min_match=NMIN+1;
        else if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_E2E && genotype->is_output[gene_id]==OUTPUT_PROTEIN)
    {
        if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_T2E && genotype->is_output[gene_id]==OUTPUT_PROTEIN)
    {
        if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==NON_OUTPUT_PROTEIN && protein_id!=start_TF)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_EE2ET)
    {
        if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;                    
    }
    else if(genotype->flag_rm_which_TFBS==RM_ES2ET)
    {
        if(protein_id==start_TF && genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN && genotype->is_output[gene_id]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_EES2ETT)
    {
        if(protein_id==start_TF && genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else if(genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_TE2ET)
    {
        if(genotype->is_output[gene_id]==OUTPUT_PROTEIN && protein_id!=start_TF && genotype->is_output[genotype->protein_pool[protein_id][1][0]]==NON_OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else if(genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN && genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_TS2ET)
    {
        if(genotype->is_output[gene_id]==OUTPUT_PROTEIN && protein_id!=start_TF && genotype->is_output[genotype->protein_pool[protein_id][1][0]]==NON_OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else if(protein_id==start_TF && genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_TE2E)
    {
        if(genotype->is_output[gene_id]==OUTPUT_PROTEIN && protein_id!=start_TF)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_EST2TTE)
    {
        if(genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN && (protein_id==start_TF || genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN))
            *N_min_match=NMIN+1;
        else if(genotype->is_output[gene_id]==OUTPUT_PROTEIN && protein_id!=start_TF && genotype->is_output[genotype->protein_pool[protein_id][1][0]]==NON_OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_EET2TEE)
    {
        if(genotype->is_output[gene_id]==OUTPUT_PROTEIN && protein_id!=start_TF)
            *N_min_match=NMIN+1;
        else if(genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN && genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }
    else if(genotype->flag_rm_which_TFBS==RM_SET2TEE)
    {
        if(genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN && protein_id==start_TF)
            *N_min_match=NMIN+1;
        else if(genotype->is_output[gene_id]==OUTPUT_PROTEIN && protein_id!=start_TF)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;

    }
    else if(genotype->flag_rm_which_TFBS==RM_ESET2TTEE)
    {
        if(genotype->is_output[gene_id]==NON_OUTPUT_PROTEIN && (protein_id==start_TF || genotype->is_output[genotype->protein_pool[protein_id][1][0]]==OUTPUT_PROTEIN))
            *N_min_match=NMIN+1;
        else if(genotype->is_output[gene_id]==OUTPUT_PROTEIN && protein_id!=start_TF)
            *N_min_match=NMIN+1;
        else
            *N_min_match=NMIN;
    }            
}