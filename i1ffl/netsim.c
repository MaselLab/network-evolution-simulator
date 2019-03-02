/* -*- Mode: C; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe, Kun Xiong
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
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

/*initial conditions*/
int init_N_non_output_act=3;
int init_N_non_output_rep=3;
int init_N_output_act=1;
int init_N_output_rep=0;


/******************************************************************************
 * 
 *                     Private function prototypes
 *
 *****************************************************************************/
static void initialize_sequence(char *, int, int, RngStream);

static void initialize_genotype_fixed(Genotype *, int, int, int, int, RngStream);

static void calc_avg_fitness(Genotype *, Selection *, int [MAX_GENES], float [MAX_PROTEINS], RngStream [N_THREADS], float *, float *); 

static void clone_genotype(Genotype *, Genotype *);

static void try_replacement(Genotype *, Genotype *, int *, float*);

static void summarize_binding_sites(Genotype *,int);

static void set_signal(CellState *, Environment *, float, RngStream, int);

static void output_genotype(Genotype *);

static int evolve_N_steps(Genotype *, Genotype *,  Mutation *, Selection *, Output_buffer [OUTPUT_INTERVAL], int *, int *, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS], int);

static void run_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_PROTEINS], int, int, RngStream, RngStream [N_THREADS]);

static void continue_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS]);

static void calc_fitness_stats(Genotype *, Selection *, float (*)[N_REPLICATES], float (*)[N_REPLICATES], int);

static float calc_replicate_fitness(CellState *, Environment *, int);

static void replay_mutations(Genotype *, Mutation *, int);

static void find_motifs(Genotype *);

static int find_TFBS_of_A_on_B(Genotype *, int, int);

static void tidy_output_files(char*, char*);

static void print_motifs(Genotype *);

static void store_mutant_info(Genotype *, Mutation *, Output_buffer *, int, int);

static void store_resident_info(Genotype *, Mutation *, Output_buffer *, int , int, int, float, int);

static void output_mutant_info(Output_buffer *, int);

static void output_resident_info(Output_buffer [OUTPUT_INTERVAL], int, int);

static void who_regulates_effector(Genotype *, int , int [MAX_PROTEINS], int [MAX_PROTEINS], int *, int *);

#if PERTURB
static void remove_edges_iteratively(Genotype *);

static void modify_topology(Genotype *, Genotype *);

static void add_binding_site(Genotype *, int);

static void remove_binding_sites(Genotype *, int, int);
#endif

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
                calc_avg_fitness(resident, burn_in, init_mRNA, init_protein, RS_parallel, GR1[i], GR2[i]);
            calc_fitness_stats(resident,burn_in,&(GR1[0]),&(GR2[0]),HI_RESOLUTION_RECALC); 
        }
        else
        {
            for(i=0;i<HI_RESOLUTION_RECALC;i++)  
                calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, GR1[i], GR2[i]);  
            calc_fitness_stats(resident,selection,&(GR1[0]),&(GR2[0]),HI_RESOLUTION_RECALC); 
        }
       
        /* make title of the output file*/
        fp=fopen(evo_summary,"w");
        fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step N_hit_bound accepted_mut selection_coeff avg_fitness fitness1 fitness2 se_avg_fitness se_fitness1 se_fitness2 N_genes N_proteins N_activator N_repressor\n");
        fprintf(fp,"0 0 0 0 na na %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d \n",  
                resident->avg_fitness,               
                resident->fitness1,
                resident->fitness2,
                resident->SE_avg_fitness,
                resident->SE_fitness1,
                resident->SE_fitness2,
                resident->ngenes,
                resident->nproteins,
                resident->N_act,
                resident->N_rep);
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
    calc_all_binding_sites(resident);
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
        calc_all_binding_sites(resident);
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
        calc_all_binding_sites(resident);
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
void show_phenotype(Genotype *resident, Genotype *mutant, Mutation *mut_record, Selection *selection, int init_mRNA[MAX_GENES], float init_protein[MAX_GENES], RngStream RS_parallel[N_THREADS])
{   
    /*sampling parameters of network motifs*/
  //  if(SAMPLE_PARAMETERS)
 //       sample_motifs(resident, mut_record, selection->MAX_STEPS, RS_parallel[0]);
    
    /*replay mutations, output N_motifs.txt and networks.txt*/   
    if(REPRODUCE_GENOTYPES || SAMPLE_GENE_EXPRESSION)
    {        
        replay_mutations(resident, mut_record, selection->MAX_STEPS);    
        /*output the evolved genotype*/
        calc_all_binding_sites(resident); 
        print_mutatable_parameters(resident,1);
        summarize_binding_sites(resident,selection->MAX_STEPS);  
    }
    
    if(SAMPLE_GENE_EXPRESSION)
    {
        /*create threads*/
        omp_set_num_threads(N_THREADS);  

        int i,j;
        float fold_change[10]={10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
        float init_sig_strength[4]={10.0,40.0,160.0,640.0};
        char fold_name[32];

        for(i=0;i<10;i++)
        {
            snprintf(fold_name,sizeof(char)*32,"%i",(int)fold_change[i]);

            mkdir(fold_name,0700);
            chdir(fold_name);

            for(j=0;j<4;j++)
            {
                snprintf(fold_name,sizeof(char)*32,"%i",(int)init_sig_strength[j]);
                mkdir(fold_name,0700);
                chdir(fold_name);
                /*collection interval is 1 minute by default*/   
                selection->env1.t_development=359.9; //minutes
                selection->env2.t_development=359.9;
                selection->env1.signal_on_strength=init_sig_strength[j]*fold_change[i];    
                selection->env1.signal_off_strength=init_sig_strength[j];
                selection->env2.signal_on_strength=init_sig_strength[j]*fold_change[i]*2.0;
                selection->env2.signal_off_strength=init_sig_strength[j]*2.0;    
                selection->env1.t_signal_on=200.0; //the signal starts with "on" and last for 200 minutes, longer than the duration of developmental simulation, which means the signal is effective constant "on" 
                selection->env1.t_signal_off=180.0; 
                selection->env2.t_signal_on=200.0; //the signal is "on" for the first 10 minutes in a developmental simulation of env2     
                selection->env2.t_signal_off=180.0;

                selection->env1.t_development+=0.2;
                selection->env2.t_development+=0.2;    
                calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, NULL,NULL);   
                chdir("..");
            }    
            chdir("..");
        }

        /*collection interval is 1 minute by default*/    
    //    selection->env1.t_development+=0.2;
    //    selection->env2.t_development+=0.2;    
    //    calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, NULL,NULL);    
    }
}
#endif

#if PERTURB
void modify_network(Genotype *resident,
                    Genotype *mutant, 
                    Mutation *mut_record,  
                    Selection *selection,
                    int init_mRNA[MAX_GENES],
                    float init_protein[MAX_PROTEINS],
                    RngStream RS_parallel[N_THREADS])
{
    int i,j,k;  
    int N_motifs;
    char buffer[600],char_buffer;
    int int_buffer,step;
    float float_buffer, mean_overall_fitness, mean_fitness1, mean_fitness2, se_overall_fitness, se_fitness1, se_fitness2;
    float fitness1[HI_RESOLUTION_RECALC][N_REPLICATES],fitness2[HI_RESOLUTION_RECALC][N_REPLICATES]; 
    FILE *file_mutation,*fitness_record,*f_aft_perturbation,*f_bf_perturbation;
    
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
    
    /*skip first 2 rows of fitness_record*/
    fitness_record=fopen(evo_summary,"r");
    fgets(buffer,600,fitness_record);
    fgets(buffer,600,fitness_record);
    
    /*file to store the original fitness of the to-be-modified network*/
    f_bf_perturbation=fopen("f_bf_perturbation.txt","w");
    
    /*create threads*/
    omp_set_num_threads(N_THREADS);
    
    /*begin*/
    for(i=1;i<=selection->MAX_STEPS;i++)
    {               
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
        
        fscanf(fitness_record,"%d %d %d %d %c %f %f %f %f %f %f %f %d %d %d %d\n",
                &step,
                &int_buffer,
                &int_buffer,
                &int_buffer,
                &char_buffer,
                &float_buffer,
                &mean_overall_fitness,
                &se_overall_fitness,
                &mean_fitness1,
                &mean_fitness2,
                &se_fitness1,
                &se_fitness2,
                &int_buffer,
                &int_buffer,
                &int_buffer,
                &int_buffer);
      
        if(i>=selection->MAX_STEPS-9999)
        {            
            calc_all_binding_sites(resident);  
       //     summarize_binding_sites(resident,i); 
            find_motifs(resident); 
            N_motifs=0;
            for(j=0;j<17;j++)
                N_motifs+=resident->N_motifs[j];
#if IGNORE_BS
            if(resident->N_motifs[1]!=0 && resident->N_motifs[1]==N_motifs)
#elif MERGE_PROTEIN
            if(resident->N_motifs[13]!=0 && resident->N_motifs[11]==resident->N_motifs[13])   
#endif
            {
                for(j=0;j<HI_RESOLUTION_RECALC;j++) 
                    calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, fitness1[j], fitness2[j]);
                calc_fitness_stats( resident, selection, &(fitness1[0]), &(fitness2[0]), HI_RESOLUTION_RECALC);  
                f_aft_perturbation=fopen("f_aft_perturbation.txt","a+");
                fprintf(f_aft_perturbation,"%d %.10f %.10f %.10f %.10f %.10f %.10f\n",i,
                        resident->avg_fitness,                        
                        resident->fitness1,
                        resident->fitness2,
                        sqrt(resident->SE_avg_fitness),
                        sqrt(resident->SE_fitness1),
                        sqrt(resident->SE_fitness2));
                fclose(f_aft_perturbation);
                fprintf(f_bf_perturbation,"%d %.10f %.10f %.10f %.10f %.10f %.10f\n",
                        step,
                        mean_overall_fitness,
                        mean_fitness1,
                        mean_fitness2,
                        se_overall_fitness,
                        se_fitness1,
                        se_fitness2);                
            }    
            for(j=0;j<MAX_GENES;j++)
            {
                resident->gene_in_core_C1ffl[j]=0;
                for(k=0;k<MAX_PROTEINS;k++)
                    resident->TF_in_core_C1ffl[j][k]=0;
            }        
        }
    }  
    fclose(f_bf_perturbation);
    fclose(file_mutation);
    fclose(fitness_record);
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
    for(j=0;j<MAX_OUTPUT_PROTEINS;j++)        
        genotype->output_protein_ids[j]=NA;
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
    /* put output genes to the end*/
    for(i=0;i<genotype->n_output_genes;i++)
        genotype->output_protein_ids[i]=genotype->nproteins-genotype->n_output_genes+i;
    
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
#if EFFECTOR_NOT_TF
            if(genotype->is_output[k]==NON_OUTPUT_PROTEIN)
            {                
#endif
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
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_binding_seq_rc[j-i]) match_rc++;

                if (match_rc >= NMIN)
                {
                    /**********************************************************************/     
                    if (N_binding_sites + 1 >= genotype->N_allocated_elements) 
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
                    if(genotype->locus_specific_TF_behavior[gene_id][k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;

                }
            } 
#if EFFECTOR_NOT_TF
            }
#endif
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
            calc_all_binding_sites_copy(genotype,gene_id);          
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
        genotype->output_protein_ids[j]=i;
        j++;
    }
    for(;i<init_N_non_output_rep+init_N_non_output_act+init_N_output_act+init_N_output_rep+N_SIGNAL_TF;i++)
    {
        genotype->protein_identity[i]=REPRESSOR;
        genotype->is_output[i]=OUTPUT_PROTEIN;
        genotype->output_protein_ids[j]=i;
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
 * Set how the environmental signal should change
 */
static void set_signal(CellState *state, Environment *env, float t_burn_in, RngStream RS, int thread_ID)
{
    float t=0.0;     
    char flag;   
    
#if IRREG_SIGNAL
    int j;
    j=RngStream_RandInt(RS,0,99);
    env->external_signal=&(signal_profile_matrix[thread_ID][j][0]);
#else
    env->external_signal=NULL;             
#endif    
    
    if(env->external_signal==NULL)   
    {
        /*always start a burn-in with signal off*/       
        state->protein_number[N_SIGNAL_TF-1]=env->signal_off_strength; 
        if(t_burn_in!=0.0)
            t=t+t_burn_in; //the completion of burn_in is a fixed event, which is added in initialize_cell       
        /*after burn-in, signal should be turned "o"n*/
        flag='f'; 
#if N_SIGNAL_TF==2
        state->protein_number[0]=background_signal_strength;
#endif      
        while(t<env->t_development+t_burn_in)
        {
            if(flag=='o')
            {                
                if(env->t_signal_on!=0.0) 
                {
                    /*add a fixed event to TURN OFF signal.*/
                    add_fixed_event(-1,t+env->t_signal_on,&(state->signal_off_head),&(state->signal_off_tail));
                    t=t+env->t_signal_on; 
                }
                flag='f';                                  
            }    
            else
            {
                if(env->t_signal_off!=0.0)
                {
                    /*add when to TURN ON signal*/
                    add_fixed_event(-1,t+env->t_signal_off,&(state->signal_on_head),&(state->signal_on_tail));
                    t=t+env->t_signal_off;
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
        genotype_clone->which_cluster[i]=-1;
    /*copy which_cluster and cis-reg sequence*/
    for(i=0; i< genotype_templet->ngenes;i++)
    {
        genotype_clone->which_cluster[i]=genotype_templet->which_cluster[i];            
        memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_templet->cisreg_seq[i][0],CISREG_LEN*sizeof(char));                    
        genotype_clone->recalc_TFBS[i]=1;                
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
        genotype_clone->output_protein_ids[i]=NA;
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
        genotype_clone->output_protein_ids[i]=genotype_templet->output_protein_ids[i];           
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
}

/**
 *Calculate the fintess of a given genotype.
 *Essentially calling do_single_timestep until tdevelopment and calculate 
 *average growth rate over tdevelopment.
 */
static void calc_avg_fitness(   Genotype *genotype,
                                Selection *Selection,
                                int init_mRNA[MAX_GENES],
                                float init_protein_number[MAX_PROTEINS],
                                RngStream RS_parallel[N_THREADS],           
                                float Fitness1[N_REPLICATES],
                                float Fitness2[N_REPLICATES])         
{ 
    Phenotype timecourse1[N_REPLICATES], timecourse2[N_REPLICATES]; 
#if PHENOTYPE     
    int i,j;   
    /*alloc space and initialize values to 0.0*/
    for(i=0;i<N_REPLICATES;i++)
    {
        timecourse1[i].total_time_points=(int)Selection->env1.t_development;
        timecourse1[i].gene_specific_concentration=(float *)malloc(timecourse1[i].total_time_points*genotype->ngenes*sizeof(float));
        timecourse1[i].protein_concentration=(float *)malloc(timecourse1[i].total_time_points*genotype->nproteins*sizeof(float));
        timecourse1[i].instantaneous_fitness=(float *)malloc(timecourse1[i].total_time_points*sizeof(float));
        timecourse1[i].timepoint=0;
        for(j=0;j<timecourse1[i].total_time_points*genotype->ngenes;j++)
            timecourse1[i].gene_specific_concentration[j]=0.0;
        for(j=0;j<timecourse1[i].total_time_points*genotype->nproteins;j++)
            timecourse1[i].protein_concentration[j]=0.0;
        for(j=0;j<timecourse1[i].total_time_points;j++)
            timecourse1[i].instantaneous_fitness[j]=0.0;
        /*do the same to timecourse2*/
        timecourse2[i].total_time_points=(int)Selection->env2.t_development;
        timecourse2[i].gene_specific_concentration=(float *)malloc(timecourse2[i].total_time_points*genotype->ngenes*sizeof(float));
        timecourse2[i].protein_concentration=(float *)malloc(timecourse2[i].total_time_points*genotype->nproteins*sizeof(float));
        timecourse2[i].instantaneous_fitness=(float *)malloc(timecourse2[i].total_time_points*sizeof(float));
        timecourse2[i].timepoint=0;
        for(j=0;j<timecourse2[i].total_time_points*genotype->ngenes;j++)
            timecourse2[i].gene_specific_concentration[j]=0.0;
        for(j=0;j<timecourse2[i].total_time_points*genotype->nproteins;j++)
            timecourse2[i].protein_concentration[j]=0.0;
        for(j=0;j<timecourse2[i].total_time_points;j++)
            timecourse2[i].instantaneous_fitness[j]=0.0;
        timecourse1[i].max_change_in_probability_of_binding=0.0;
        timecourse2[i].max_change_in_probability_of_binding=0.0;
    }        
#endif
    
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
        float t_burn_in;
        float init_protein_number_clone[MAX_GENES];
        float f1[N_replicates_per_thread],f2[N_replicates_per_thread];        
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];         
        Environment Env1, Env2;

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
            Env1.t_development=Selection->env1.t_development;
            Env1.signal_on_strength=Selection->env1.signal_on_strength;
            Env1.signal_off_strength=Selection->env1.signal_off_strength;
            Env1.t_signal_on=Selection->env1.t_signal_on;
            Env1.t_signal_off=Selection->env1.t_signal_off;
            Env1.initial_effect_of_effector=Selection->env1.initial_effect_of_effector;
            Env1.effect_of_effector_aft_burn_in=Selection->env1.effect_of_effector_aft_burn_in;
            Env1.fixed_effector_effect=Selection->env1.fixed_effector_effect;
            Env1.max_duration_of_burn_in_growth_rate=Selection->env1.max_duration_of_burn_in_growth_rate;
            Env1.avg_duration_of_burn_in_growth_rate=Selection->env1.avg_duration_of_burn_in_growth_rate;
            Env1.min_peak_response=Selection->env1.min_peak_response;
            Env1.min_reduction_relative_to_peak=Selection->env1.min_reduction_relative_to_peak;
            Env1.min_response_aft_signal_change=Selection->env1.min_response_aft_signal_change;
            Env1.max_response_bf_signal_change=Selection->env1.max_response_bf_signal_change;
            Env1.window_size=Selection->env1.window_size;
            
            Env2.t_development=Selection->env2.t_development;
            Env2.signal_on_strength=Selection->env2.signal_on_strength;
            Env2.signal_off_strength=Selection->env2.signal_off_strength;
            Env2.t_signal_on=Selection->env2.t_signal_on;
            Env2.t_signal_off=Selection->env2.t_signal_off;
            Env2.initial_effect_of_effector=Selection->env2.initial_effect_of_effector;
            Env2.effect_of_effector_aft_burn_in=Selection->env2.effect_of_effector_aft_burn_in;
            Env2.fixed_effector_effect=Selection->env2.fixed_effector_effect;
            Env2.max_duration_of_burn_in_growth_rate=Selection->env2.max_duration_of_burn_in_growth_rate;
            Env2.avg_duration_of_burn_in_growth_rate=Selection->env2.avg_duration_of_burn_in_growth_rate;                   
            Env2.min_peak_response=Selection->env2.min_peak_response;
            Env2.min_reduction_relative_to_peak=Selection->env2.min_reduction_relative_to_peak;
            Env2.min_response_aft_signal_change=Selection->env2.min_response_aft_signal_change;
            Env2.max_response_bf_signal_change=Selection->env2.max_response_bf_signal_change;
            Env2.window_size=Selection->env2.window_size;
            
            for(j=0; j < MAX_GENES; j++) 
            {  
                init_mRNA_clone[j] = init_mRNA[j];
                init_protein_number_clone[j] = init_protein_number[j];
            } 
        } 
        calc_all_binding_sites(&genotype_clone); 
        
#if PERTURB 
        modify_topology(genotype, &genotype_clone);
#endif
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
         *                              TEST1 
         *
         *********************************************************************/
        for(i=0;i<N_replicates_per_thread;i++) /* env 1, usually a constant signal that matches env*/
        {  
            /*make a t_burn_in before turning on signal*/
            do
                t_burn_in=Env1.avg_duration_of_burn_in_growth_rate*expdev(RS_parallel[thread_ID]);
            while(t_burn_in>Env1.max_duration_of_burn_in_growth_rate);  
            /*initialize mRNA and protein numbers, and gene states etc.*/
            initialize_cell(&genotype_clone, &state_clone, &Env1, t_burn_in, mRNA, protein);
            
            /*set how the signal should change during simulation*/
            set_signal(&state_clone, &Env1, t_burn_in, RS_parallel[thread_ID], thread_ID);
            
            /*calcualte the rates of cellular activity based on the initial cellular state*/
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, &Env1, &(timecourse1[thread_ID*N_replicates_per_thread+i]), t_burn_in, INITIALIZATION);             
#if PHENOTYPE
            timecourse1[thread_ID*N_replicates_per_thread+i].timepoint=0;
#endif
             /*run developmental simulation until tdevelopment or encounter an error*/
            while(state_clone.t<Env1.t_development+t_burn_in) 
                do_single_timestep(&genotype_clone, &state_clone, &rate_clone, &Env1, t_burn_in, &(timecourse1[thread_ID*N_replicates_per_thread+i]), RS_parallel[thread_ID]);            

#if !PHENOTYPE                                     
            f1[i]=calc_replicate_fitness(&state_clone,&Env1,genotype_clone.n_output_genes);        
            /*free linked tables*/
#if EVOLVE_I1FFL
            free(state_clone.sampled_response);
#endif
#endif    
            free_fixedevent(&state_clone);           
        }
        /********************************************************************** 
         * 
         *                              TEST2 
         *
         *********************************************************************/
        for(i=0;i<N_replicates_per_thread;i++) 
        { 
            do
                t_burn_in=Env2.avg_duration_of_burn_in_growth_rate*expdev(RS_parallel[thread_ID]);
            while(t_burn_in>Env2.max_duration_of_burn_in_growth_rate);            
            initialize_cell(&genotype_clone, &state_clone, &Env2, t_burn_in, mRNA, protein);
            set_signal(&state_clone, &Env2, t_burn_in, RS_parallel[thread_ID], thread_ID);
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, &Env2, &(timecourse2[thread_ID*N_replicates_per_thread+i]), t_burn_in, INITIALIZATION); 
#if PHENOTYPE
            timecourse2[thread_ID*N_replicates_per_thread+i].timepoint=0;
#endif            
            while(state_clone.t<Env2.t_development+t_burn_in) 
                do_single_timestep(&genotype_clone, &state_clone, &rate_clone, &Env2, t_burn_in, &(timecourse2[thread_ID*N_replicates_per_thread+i]), RS_parallel[thread_ID]);            
        
#if !PHENOTYPE           
            f2[i]=calc_replicate_fitness(&state_clone,&Env2,genotype_clone.n_output_genes);                      
#if EVOLVE_I1FFL
            free(state_clone.sampled_response);
#endif
#endif
            free_fixedevent(&state_clone);            
        }    
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
                Fitness2[i]=f2[j];
                j++;
            }
        }
    }     
#if PHENOTYPE
    /*output timecourse*/
    int k;
    char filename[32];
    FILE *fp;   
    /*mkdir*/
    //mkdir("phenotype",0700);
	//chdir("phenotype");
    /*output output_gene ids*/
	fp=fopen("output_gene_ids.txt","w");
	for(i=0;i<genotype->n_output_genes;i++)
	fprintf(fp,"%d\n",genotype->output_protein_ids[i]);
	fclose(fp);	
    fp=fopen("output_node_ids.txt","w");
    for(i=0;i<genotype->N_node_families;i++)
    {
        if(genotype->is_output[genotype->node_family_pool[i][1][0]]==OUTPUT_PROTEIN)
            fprintf(fp,"%d\n",i);
    }
    fclose(fp);
    /*fitness: each row is a replicate*/
//    fp=fopen("fitnessA","w");
//    for(i=0;i<N_REPLICATES;i++)
//    {
//        for(j=0;j<timecourse1[i].total_time_points;j++)
//            fprintf(fp,"%f ",timecourse1[i].instantaneous_fitness[j]);
//        fprintf(fp,"\n");
//    }
//    fclose(fp);
//    fp=fopen("fitnessB","w");
//    for(i=0;i<N_REPLICATES;i++)
//    {
//        for(j=0;j<timecourse2[i].total_time_points;j++)
//            fprintf(fp,"%f ",timecourse2[i].instantaneous_fitness[j]);
//        fprintf(fp,"\n");
//    }
//    fclose(fp);
    /*proteint concentration: each protein has its own file, in which each row is a replicate*/    
    for(i=0;i<genotype->N_node_families;i++)
    {
        snprintf(filename,sizeof(char)*32,"protein%i_A",i);
        fp=fopen(filename,"w");
        for(j=0;j<N_REPLICATES;j++)
        {
            for(k=0;k<timecourse1[j].total_time_points;k++)                
                fprintf(fp,"%f ",timecourse1[j].protein_concentration[k+i*timecourse1[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    for(i=0;i<genotype->N_node_families;i++)
    {
        snprintf(filename,sizeof(char)*32,"protein%i_B",i);
        fp=fopen(filename,"w");
        for(j=0;j<N_REPLICATES;j++)
        {
            for(k=0;k<timecourse2[j].total_time_points;k++)                
                fprintf(fp,"%f ",timecourse2[j].protein_concentration[k+i*timecourse2[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    /*gene-specific concentration: each protein has its own file, in which each row is a replicate*/
//     for(i=0;i<genotype->ngenes;i++)
//    {
//        snprintf(filename,sizeof(char)*32,"gene%i_A",i);
//        fp=fopen(filename,"w");
//        for(j=0;j<N_REPLICATES;j++)
//        {
//            for(k=0;k<timecourse1[j].total_time_points;k++)                
//                fprintf(fp,"%f ",timecourse1[j].gene_specific_concentration[k+i*timecourse1[j].total_time_points]);
//            fprintf(fp,"\n");
//        }
//        fclose(fp);
//    }
//    for(i=0;i<genotype->ngenes;i++)
//    {
//        snprintf(filename,sizeof(char)*32,"gene%i_B",i);
//        fp=fopen(filename,"w");
//        for(j=0;j<N_REPLICATES;j++)
//        {
//            for(k=0;k<timecourse2[j].total_time_points;k++)                
//                fprintf(fp,"%f ",timecourse2[j].gene_specific_concentration[k+i*timecourse2[j].total_time_points]);
//            fprintf(fp,"\n");
//        }
//        fclose(fp);
//    } 
    
    /*output the maximum change in the probabilities of TF binding*/
    fp=fopen("max_change_in_binding_probability_A.txt","w");
    for(i=0;i<N_REPLICATES;i++)
        fprintf(fp,"%f\n",timecourse1[i].max_change_in_probability_of_binding);
    fclose(fp);
    
    fp=fopen("max_change_in_binding_probability_B.txt","w");
    for(i=0;i<N_REPLICATES;i++)
        fprintf(fp,"%f\n",timecourse2[i].max_change_in_probability_of_binding);
    fclose(fp);      
    
    for(i=0;i<N_THREADS;i++)
    {
        free(timecourse1[i].gene_specific_concentration);
        free(timecourse2[i].gene_specific_concentration);
        free(timecourse1[i].instantaneous_fitness);
        free(timecourse2[i].instantaneous_fitness);
        free(timecourse1[i].protein_concentration);
        free(timecourse2[i].protein_concentration);
    }
    
    for(i=0;i<N_THREADS;i++)
    {
        free(timecourse1[i].gene_specific_concentration);
        free(timecourse2[i].gene_specific_concentration);
        free(timecourse1[i].instantaneous_fitness);
        free(timecourse2[i].instantaneous_fitness);
        free(timecourse1[i].protein_concentration);
        free(timecourse2[i].protein_concentration);
    }
#endif
}

static float calc_replicate_fitness(CellState *state, Environment *env, int N_output)
{
    float fitness; 
    float steady_state_response, peak_response, min_response_bf_peak, max_response_after_peak, half_response;
    float deviation;
    int pos_peak;
    float pos_half_response;   
    int i;

#if EVOLVE_I1FFL     
    /*find peak*/
    find_max(&(state->sampled_response[0]),0,state->N_samples,&peak_response,&pos_peak); 
    /* if expression is "flat" or monotonous decrease */
    if(peak_response==0.0 || fabs(peak_response-state->sampled_response[0])/state->sampled_response[0]<=EPSILON)    
        fitness=0.0-exp_cost_factor*state->cumulative_cost; 
    else
    {   /*calculate the average signal strength at the end of simulation*/
        steady_state_response=0.0;
        for(i=state->N_samples-1;i>state->N_samples-1-env->window_size;i--)
            steady_state_response+=state->sampled_response[i];
        steady_state_response/=(float)env->window_size;        
        /* if select for response acceleration*/
        if(SELECTION==0) 
        {
            /*time to reach 50% ss level*/
            half_response=0.5*steady_state_response;
            find_x(&(state->sampled_response[0]),0,pos_peak,half_response,&pos_half_response,0); 
            /*faster is better*/
            fitness+=(env->t_development-env->t_signal_off-pos_half_response)/(env->t_development-env->t_signal_off);  

            /*response cannot be too low*/       
            fitness+=(steady_state_response<env->min_response_aft_signal_change)?steady_state_response/env->min_response_aft_signal_change:1.0;

            /*select for low ss level bf signal change*/
            fitness=(state->sampled_response[0]<env->max_response_bf_signal_change)?1.0:(env->min_peak_response-state->sampled_response[0])/(env->min_peak_response-env->max_response_bf_signal_change);

            /*reaching steady state?*/
//            deviation=0.0;
//            for(i=state->N_samples-1;i>state->N_samples-1-env->window_size;i--)
//                deviation+=fabs(state->sampled_response[i]-steady_state_response);
//            deviation/=steady_state_response;
//            fitness+=(deviation<env->max_relative_deviation_from_mean)?1.0:
            /*adds up*/
            fitness=fitness/3.0-exp_cost_factor*state->cumulative_cost;
        }
        else //select for pulse-based dynamics
        {  
            if(pos_peak==state->N_samples-1)//monotonous increase
                fitness=0.0-exp_cost_factor*state->cumulative_cost;
            else
            {
                /*select for peak level*/
                fitness+=(peak_response>env->min_peak_response)?1.0:peak_response/env->min_peak_response;

                /* look for mid point bf the peak*/
                half_response=(peak_response+state->sampled_response[0])*0.5;      
                find_x(&(state->sampled_response[0]),0,pos_peak,half_response,&pos_half_response,1); 
                /*faster is better*/
                fitness+=(env->t_development-env->t_signal_off-pos_half_response)/(env->t_development-env->t_signal_off);  

                /*select for low ss level bf signal change*/
                fitness=(state->sampled_response[0]<env->max_response_bf_signal_change)?1.0:(env->min_peak_response-state->sampled_response[0])/(env->min_peak_response-env->max_response_bf_signal_change);

                /*low ss level after signal change*/
                if(SELECTION==1) // selection for just pulse
                {
                    max_response_after_peak=(1.0-env->min_reduction_relative_to_peak)*peak_response;
                    fitness+=(steady_state_response<max_response_after_peak)?1.0:(peak_response-steady_state_response)/(peak_response-max_response_after_peak);
                }
                else                
                    fitness=(state->sampled_response[0]<env->max_response_bf_signal_change)?1.0:(env->min_peak_response-state->sampled_response[0])/(env->min_peak_response-env->max_response_bf_signal_change);
                
                /*adds up*/
                fitness=fitness/4.0-exp_cost_factor*state->cumulative_cost; 
            }        
        }  
    }
#else    
    int i;
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
    int i, output_counter;
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
    first_step=init_step;
    N_tot_trials=init_N_tot_mutations;
    
    /* first, run burn-in */
    if(burn_in->MAX_STEPS!=0)
    {
        flag_burn_in=1; 
        DUPLICATION=burn_in->temporary_DUPLICATION;                 
        SILENCING=burn_in->temporary_SILENCING;
//        N_EFFECTOR_GENES=burn_in->temporary_N_effector_genes;
//        N_TF_GENES=burn_in->temporary_N_tf_genes; 
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
                                    init_mRNA,
                                    init_protein,
                                    RS_parallel,                                        
                                    fitness1[i],
                                    fitness2[i]); 
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
//    N_EFFECTOR_GENES=selection->temporary_N_effector_genes;
//    N_TF_GENES=selection->temporary_N_tf_genes; 
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
    char buffer[200]; 
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
            fgets(buffer,200,fp);
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
            calc_avg_fitness(mutant, selection, init_mRNA, init_protein, RS_parallel, fitness1[0], fitness2[0]);
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
                calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, fitness1[j],fitness2[j]);              
            calc_fitness_stats(resident, selection, &(fitness1[0]), &(fitness2[0]), HI_RESOLUTION_RECALC);   
        }  
        
        /*calculate the number of c1-ffls*/
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
    for(i=0;i<17;i++)
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
    int i,j,k,which_protein,which_gene;
    int table[MAX_GENES][MAX_GENES];
    
    for(i=0;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->ngenes;j++)
            table[i][j]=0;     
    }
   
    for(i=N_SIGNAL_TF;i<genotype->N_cisreg_clusters;i++)
    {
        which_gene=genotype->cisreg_cluster_pool[i][1][0];
        for(j=0;j<genotype->N_node_families;j++)    
        {
            which_protein=genotype->which_protein[genotype->node_family_pool[j][1][0]];
            for(k=0;k<genotype->binding_sites_num[which_gene];k++)
            {
                if(genotype->all_binding_sites[which_gene][k].tf_id==which_protein)
                    table[which_gene][j]++; /*the numbers of the binding sites of each "node" on promoter of which_gene*/   
            }   
        }
       
        for(j=0;j<genotype->cisreg_cluster_pool[i][1][j];j++)
        {            
            for(k=0;k<genotype->N_node_families;k++)
                table[genotype->cisreg_cluster_pool[i][1][j]][k]=table[which_gene][k];
        }
    }   
    
    /*Output all binding sites*/ 
    OUTPUT1=fopen("networks.txt","a+");
    fprintf(OUTPUT1,"step %d\n",step_i);
    fprintf(OUTPUT1,"Gene   ");    
    for(i=0;i<genotype->N_node_families;i++)
    {
        if(i<10)
        {
            if(genotype->protein_identity[genotype->which_protein[genotype->node_family_pool[i][1][0]]]==ACTIVATOR)
            {
                if(genotype->is_output[genotype->node_family_pool[i][1][0]]==OUTPUT_PROTEIN)
                    fprintf(OUTPUT1,"A%d     ",i);
                else
                    fprintf(OUTPUT1,"a%d     ",i);
            }
            if(genotype->protein_identity[genotype->which_protein[genotype->node_family_pool[i][1][0]]]==REPRESSOR)
            {
                if(genotype->is_output[genotype->node_family_pool[i][1][0]]==OUTPUT_PROTEIN)
                    fprintf(OUTPUT1,"R%d     ",i);
                else
                    fprintf(OUTPUT1,"r%d     ",i);
            }
        }
        else
        {
            if(genotype->protein_identity[genotype->which_protein[genotype->node_family_pool[i][1][0]]]==ACTIVATOR)
            {
                if(genotype->is_output[genotype->node_family_pool[i][1][0]]==OUTPUT_PROTEIN)
                    fprintf(OUTPUT1,"A%d    ",i);
                else
                    fprintf(OUTPUT1,"a%d    ",i);
            }
            if(genotype->protein_identity[genotype->which_protein[genotype->node_family_pool[i][1][0]]]==REPRESSOR)
            {
                if(genotype->is_output[genotype->node_family_pool[i][1][0]]==OUTPUT_PROTEIN)
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
                fprintf(OUTPUT1," %d(%d)  ",table[i][j],genotype->locus_specific_TF_behavior[i][genotype->which_protein[genotype->node_family_pool[j][1][0]]]);
            else
                fprintf(OUTPUT1," %d(%d) ",table[i][j],genotype->locus_specific_TF_behavior[i][genotype->which_protein[genotype->node_family_pool[j][1][0]]]);
        }
        if(genotype->is_output[i]==OUTPUT_PROTEIN)
            fprintf(OUTPUT1,"      E%d",genotype->which_node_family[i]);
        else
        {
            if(genotype->protein_identity[genotype->which_protein[i]]==ACTIVATOR)
                fprintf(OUTPUT1,"      a%d",genotype->which_node_family[i]); 
            if(genotype->protein_identity[genotype->which_protein[i]]==REPRESSOR)
                fprintf(OUTPUT1,"      r%d",genotype->which_node_family[i]);
        }
        fprintf(OUTPUT1," %d \n",genotype->min_N_activator_to_transc[i]);
    }
    fprintf(OUTPUT1,"\n"); 
    fprintf(OUTPUT1,"\n"); 
    fclose(OUTPUT1);
  
}

static void find_motifs(Genotype *genotype)
{
    int i,j,k,l;
    int found_bs;
    int n_members;
    int gene_id, N_activators, N_repressors;
    int effector_gene, effector_node, effector_protein;
    int auxiliary_tf_gene, auxiliary_tf_node, auxiliary_tf_protein;
    int repressors[MAX_PROTEINS];
    int activators[MAX_PROTEINS];
    int regulated_by_signal[MAX_GENES];   
    
    /*reset variables*/       
    for(i=0;i<MAX_GENES;i++)
    {
        genotype->gene_in_core_C1ffl[i]=0;
        for(j=0;j<MAX_PROTEINS;j++)
            genotype->TF_in_core_C1ffl[i][j]=0;
    }     
    for(i=0;i<17;i++)    
        genotype->N_motifs[i]=0;  
    
    /*begin searching motifs*/    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        found_bs=0;
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {
            if(genotype->all_binding_sites[i][j].tf_id==N_SIGNAL_TF-1)                
            {
                found_bs=1;
                    break;
            }                               
        }
        if(found_bs)
            regulated_by_signal[i]=1;
        else
            regulated_by_signal[i]=0;
    } 
    
    for(i=0;i<genotype->N_cisreg_clusters;i++)
    {  
        /*does this cisreg_cluster contain output gene*/
        n_members=0;
        for(j=0;j<genotype->cisreg_cluster_pool[i][0][0];j++)
        {
            gene_id=genotype->cisreg_cluster_pool[i][1][j];
            if(genotype->is_output[gene_id]==OUTPUT_PROTEIN)
                n_members++;                
        }
        /*YES*/
        if(n_members>0)
        {
            effector_gene=gene_id;
            who_regulates_effector(genotype,effector_gene,activators,repressors,&N_activators,&N_repressors); 

            /***count ffls and NFB ***/  
            if(activators[0]==N_SIGNAL_TF-1) //effector is activated by the signal
            {
                effector_node=genotype->which_node_family[effector_gene]; //gene_id_copy encodes an effector 
                effector_protein=genotype->which_protein[effector_gene];
                /*check motifs formed with a repressor of the effector*/
                for(j=0;j<N_repressors;j++) 
                {
                    auxiliary_tf_node=repressors[j];
                    if(effector_node!=auxiliary_tf_node)
                    {
                        for(k=0;k<genotype->node_family_pool[auxiliary_tf_node][0][0];k++)
                        {
                            auxiliary_tf_gene=genotype->node_family_pool[auxiliary_tf_node][1][k]; 
                            auxiliary_tf_protein=genotype->which_protein[auxiliary_tf_gene];
                            if(regulated_by_signal[auxiliary_tf_gene]==1) // aux. TF gene is regulated by the signal
                            {
                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][N_SIGNAL_TF-1]==ACTIVATOR) //the signal can activate aux. tf
                                {                                    
                                    if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene)) // aux. tf is regulated by the effector
                                    {                                    
                                        if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) // the effector is an activator to aux. gene
                                        {  
                                            if(find_TFBS_of_A_on_B(genotype,auxiliary_tf_gene,auxiliary_tf_gene)) //aux. tf does bind to itself
                                            {
                                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][auxiliary_tf_protein]==ACTIVATOR) //aux. tf can activates itself? 
                                                {
                                                    genotype->N_motifs[1]++; //an I1-FFL+NBF+auto-activation
#if PERTURB
#if RM_PF
                                                    genotype->gene_in_core_C1ffl[auxiliary_tf_gene]=1;
                                                    genotype->TF_in_core_C1ffl[auxiliary_tf_gene][auxiliary_tf_protein]=1;
#elif RM_I1FFL
                                                    genotype->gene_in_core_C1ffl[auxiliary_tf_gene]=1;
                                                    genotype->TF_in_core_C1ffl[auxiliary_tf_gene][N_SIGNAL_TF-1]=1;                                      
#elif RM_NFBL
                                                    genotype->gene_in_core_C1ffl[auxiliary_tf_gene]=1;
                                                    genotype->TF_in_core_C1ffl[auxiliary_tf_gene][effector_protein]=1;
#elif RM_PF_NFBL
                                                    genotype->gene_in_core_C1ffl[auxiliary_tf_gene]=1;
                                                    genotype->TF_in_core_C1ffl[auxiliary_tf_gene][effector_protein]=1;
                                                    genotype->TF_in_core_C1ffl[auxiliary_tf_gene][auxiliary_tf_protein]=1;
#endif
#endif
                                                }
                                            }
                                            else
                                                genotype->N_motifs[2]++; //an I1-FFL+NBF                                            
                                        }
                                        else
                                        {
                                            genotype->N_motifs[11]++; //overlapping I1
                                            if(genotype->locus_specific_TF_behavior[effector_gene][effector_protein]==REPRESSOR &&
                                                genotype->is_output[auxiliary_tf_gene]==OUTPUT_PROTEIN)
                                            {                                                
                                                if(find_TFBS_of_A_on_B(genotype,effector_gene,effector_gene))
                                                {
                                                    genotype->N_motifs[12]++; //overlapping I1 with self-repressing effector
                                                    if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][auxiliary_tf_protein]==REPRESSOR)
                                                    {                                                        
                                                        if(find_TFBS_of_A_on_B(genotype,auxiliary_tf_gene,auxiliary_tf_gene))
                                                        {
                                                            genotype->N_motifs[13]++; //overlapping I1 with self-rep effector and self-rep aux
//#if IGNORE_BS
                                                            genotype->gene_in_core_C1ffl[effector_gene]=1;
                                                            genotype->gene_in_core_C1ffl[auxiliary_tf_gene]=1;
                                                            genotype->TF_in_core_C1ffl[effector_gene][auxiliary_tf_protein]=1;
                                                            genotype->TF_in_core_C1ffl[auxiliary_tf_gene][effector_protein]=1;
// #endif
                                                        }
                                                    }
                                                }
                                            }                                            
                                        }
                                    }
                                    else  // aux. tf is not regulated by the effector
                                    {
                                        genotype->N_motifs[0]++; // an I1FFL 
                                        if(genotype->is_output[auxiliary_tf_gene]==NON_OUTPUT_PROTEIN)
                                            genotype->N_motifs[15]++; //I1 with a non-output aux. tf
                                        else
                                            genotype->N_motifs[16]++; //I1 with an output aux. tf
                                    }
                                }                               
                            }
                            else // aux. tf is not regulated by the signal
                            {
                                /*If aux. tf is regulated by the effector?*/
                                if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene)) 
                                {                                    
                                    if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==ACTIVATOR) //if the effector is an activator to aux. gene
                                    {   
                                        if(find_TFBS_of_A_on_B(genotype,auxiliary_tf_gene,auxiliary_tf_gene)) //aux. tf does regulate itself
                                        {
                                            if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][auxiliary_tf_protein]==ACTIVATOR)//if aux. tf activates itself
                                                genotype->N_motifs[3]++; //an NBF+auto-activation
                                        }   
                                        else
                                            genotype->N_motifs[4]++; //an NBF
                                    }
                                }                                
                            }
                        }  
                    }
                }
                
                /*
                 * 
                 * check motifs formed with an activator of the effector
                 *
                 *
                 */
                for(j=N_SIGNAL_TF;j<N_activators;j++) 
                {
                    auxiliary_tf_node=activators[j];
                    if(auxiliary_tf_node!=effector_node)
                    {
                        for(k=0;k<genotype->node_family_pool[auxiliary_tf_node][0][0];k++)
                        {
                            auxiliary_tf_gene=genotype->node_family_pool[auxiliary_tf_node][1][k]; 
                            auxiliary_tf_protein=genotype->which_protein[auxiliary_tf_gene];
                            if(regulated_by_signal[auxiliary_tf_gene]==1) // aux. TF gene is regulated by the signal
                            {
                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][0]==REPRESSOR) //if the signal represses aux. tf
                                {                                   
                                    if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene)) // aux. tf is regulated by the effector
                                    {                                    
                                        if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==REPRESSOR) // the effector is an repressor to aux. gene
                                        {                                     
                                            if(find_TFBS_of_A_on_B(genotype,auxiliary_tf_gene,auxiliary_tf_gene)) //aux. tf does bind to itself
                                            {
                                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][auxiliary_tf_protein]==ACTIVATOR) //aux. tf activates itself 
                                                    genotype->N_motifs[5]++; //an I4-FFL+NBF+auto-activation
                                            }
                                            else
                                                genotype->N_motifs[6]++; //an I4-FFL+NBF
                                        }
                                    }                                   
                                }                               
                            }
                            else // aux. tf is not regulated by the signal
                            {
                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][effector_protein]==REPRESSOR) //the effector can repress aux. gene
                                {
                                    /*Is aux. tf regulated by the effector?*/
                                    if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene)) //effector has bs on aux. gene
                                    {  
                                        if(find_TFBS_of_A_on_B(genotype,auxiliary_tf_gene,auxiliary_tf_gene)) //aux. tf does regulate itself
                                        {
                                            if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][auxiliary_tf_protein]==ACTIVATOR)//if aux. tf activates itself                                            
                                                genotype->N_motifs[7]++; //an NBF+auto-activation
                                        }
                                        else
                                            genotype->N_motifs[8]++; //an NBF
                                    }    
                                }
                            }
                        }                 
                    }
                }
            } 
            
            /*
             * 
             * effector is repressed by the signal
             *
             *
             */
            if(repressors[0]==N_SIGNAL_TF-1)
            {
                effector_node=genotype->which_node_family[effector_gene]; //gene_id_copy encodes an effector                     
                for(j=N_SIGNAL_TF;j<N_activators;j++)
                {
                    auxiliary_tf_node=activators[j];
                    if(effector_node!=auxiliary_tf_node)
                    {
                        for(k=0;k<genotype->node_family_pool[auxiliary_tf_node][0][0];k++)
                        {
                            auxiliary_tf_gene=genotype->node_family_pool[auxiliary_tf_node][1][k]; 
                            auxiliary_tf_protein=genotype->which_protein[auxiliary_tf_gene];
                            if(regulated_by_signal[auxiliary_tf_gene]==1) // aux. TF gene is regulated by the signal
                            {
                                if(genotype->locus_specific_TF_behavior[auxiliary_tf_gene][0]==ACTIVATOR) //if the signal activates aux. tf
                                {                                    
                                    if(find_TFBS_of_A_on_B(genotype,effector_gene,auxiliary_tf_gene))
                                        genotype->N_motifs[9]++; //I3-FFL
                                }
                                else //the signal represses aux. tf
                                {                                   
                                    if(find_TFBS_of_A_on_B(genotype,auxiliary_tf_gene,auxiliary_tf_gene))
                                        genotype->N_motifs[10]++; //I2-FFL
                                }
                            }
                        }          
                    }
                }
            }
        }
    }     
}

static int find_TFBS_of_A_on_B(Genotype *genotype, int gene_A, int gene_B)
{
    int site_id,protein_id;
    int found_bs=0;
    protein_id=genotype->which_protein[gene_A];
    for(site_id=0;site_id<genotype->binding_sites_num[gene_B];site_id++)
    {
        if(genotype->all_binding_sites[gene_B][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_B][site_id].mis_match<=CUT_OFF_MISMATCH_TF2TF)
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
        repressors[i]=0;    
        activators[i]=0;  
    }  
    /*build protein_to_node dictionary*/
    for(i=0;i<genotype->nproteins;i++)
    {  
        /*reset*/
        for(j=0;j<genotype->N_node_families;j++)
            protein_to_node[i][1][j]=0; 
        /*which nodes does the protein contains*/
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
        for(j=0;j<protein_to_node[protein_id][0][0];j++)
        {
            if(genotype->locus_specific_TF_behavior[effector_gene][protein_id]==ACTIVATOR)
                activators[protein_to_node[protein_id][1][j]]=1;
            else
                repressors[protein_to_node[protein_id][1][j]]=1;
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
    
    fp1=fopen("N_motifs.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,2000,fp1);
        fputs(buffer,fp2);       
    }
    fclose(fp1);
    fclose(fp2);
    remove("N_motifs.txt");
    rename("temp","N_motifs.txt");
    
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

#if PERTURB
static void modify_topology(Genotype *genotype, Genotype *genotype_clone)
{
    int i, j, gene_id;   
    int likely_ancestor;
    likely_ancestor=genotype->output_protein_ids[0];
    
    for(i=0;i<17;i++)
        genotype_clone->N_motifs[i]=genotype->N_motifs[i];
    for(i=0;i<MAX_GENES;i++)
    {
        genotype_clone->gene_in_core_C1ffl[i]=genotype->gene_in_core_C1ffl[i];
        for(j=0;j<MAX_PROTEINS;j++)
            genotype_clone->TF_in_core_C1ffl[i][j]=genotype->TF_in_core_C1ffl[i][j];
    }
    
#if MERGE_PROTEIN    
    /*assume that the effector protein with the most gene copy is the ancestor of all effector proteins*/    
    for(i=1;i<genotype->n_output_proteins;i++)
        likely_ancestor=(genotype->protein_pool[likely_ancestor][0][0]>genotype->protein_pool[genotype->output_protein_ids[i]][0][0])?likely_ancestor:genotype->output_protein_ids[i];     
#endif   
    
    for(gene_id=N_SIGNAL_TF;gene_id < genotype_clone->ngenes;gene_id++)
    { 
        remove_binding_sites(genotype_clone, gene_id, likely_ancestor);
    }
}

/* ignore TF x when searching binding sites on gene y*/
/* this function is almost the same as calc_all_binding_sites_copy*/
static void remove_binding_sites(Genotype *genotype, int gene_id, int likely_ancestor)
{
    int i, j, k;
    int match,match_rc;  
    int N_hindered_BS=0;   
    int N_binding_sites=0;
    int start_TF; 
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
        
        for (k=start_TF;k<genotype->nproteins;k++) 
        { 
#if IGNORE_BS     
            if(!(genotype->gene_in_core_C1ffl[gene_id]==1 && genotype->TF_in_core_C1ffl[gene_id][k]==1))       
            {  
#endif
            tf_seq=&(genotype->tf_binding_seq[k][0]);
            tf_seq_rc=&(genotype->tf_binding_seq_rc[k][0]);    
#if MERGE_PROTEIN
            //when searching TFBS of aux tf in the cisreg of likely ancestor genes, consider aux tf has the same binding seq as the ancestor effector
            if(genotype->which_protein[gene_id]==likely_ancestor && genotype->TF_in_core_C1ffl[gene_id][k]==1)
            {
                tf_seq=&(genotype->tf_binding_seq[likely_ancestor][0]);
                tf_seq_rc=&(genotype->tf_binding_seq_rc[likely_ancestor][0]); 
            }
            //when searching TFBS of aux tf in the cisreg of likely aux tf, consider aux tf has the same binding seq as the ancestor effector
            else if(genotype->gene_in_core_C1ffl[gene_id]==1 && k==genotype->which_protein[gene_id])
            {
                tf_seq=&(genotype->tf_binding_seq[likely_ancestor][0]);
                tf_seq_rc=&(genotype->tf_binding_seq_rc[likely_ancestor][0]); 
            }
#endif
            /*find BS on the template strand*/
            match=0;
            for (j=i;j<i+TF_ELEMENT_LEN;j++) 
                if (cis_seq[j] == tf_seq[j-i]) match++; 
            if (match >= NMIN)
            {                              
                genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                      
                genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match)/(TF_ELEMENT_LEN-NMIN+1));
                genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ; 
                genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match;             
                genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                N_hindered_BS++;              
                N_binding_sites++;
                if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand.*/
            {
                match_rc=0;
                for (j=i; j < i+TF_ELEMENT_LEN; j++)                
                    if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
                if (match_rc >= NMIN)
                {                   
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                                     
                    genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(TF_ELEMENT_LEN-match_rc)/(TF_ELEMENT_LEN-NMIN+1));
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
                    genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = TF_ELEMENT_LEN-match_rc;
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                    N_hindered_BS++;                  
                    N_binding_sites++;                 
                    if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                }
            } 
#if IGNORE_BS
            }
#endif
        }/* looping through TFs ends */
    }/*end of promoter scanning*/ 
    
    genotype->binding_sites_num[gene_id]=N_binding_sites;  
    genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);
    /* calculate max_hindered_sites */    
    for(i=0;i<genotype->binding_sites_num[gene_id];i++)
    {
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id] > genotype->all_binding_sites[gene_id][i].N_hindered)?
                                      genotype->max_hindered_sites[gene_id] : genotype->all_binding_sites[gene_id][i].N_hindered;           
    }    
   
    int act_BS[MAX_TFBS_NUMBER][2],rep_BS[MAX_TFBS_NUMBER][2];
    int N_act_BS,N_rep_BS;    
    N_act_BS=1;
    N_rep_BS=1;
    for(i=0;i<genotype->binding_sites_num[gene_id];i++) /* make lists BS by their types*/    
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
    act_BS[0][0]=-1;
    act_BS[0][1]=0; 
    for(i=1;i<N_act_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][act_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][act_BS[j][0]].BS_pos<TF_ELEMENT_LEN+2*HIND_LENGTH)j--;
        act_BS[i][1]=act_BS[j][1]+1;
    } 
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

#endif

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
    for(i=0;i<17;i++)
        resident_info->n_motifs[i]=resident->N_motifs[i];
    for(i=0;i<12;i++)        
        resident_info->n_near_AND_gated_motifs[i]=resident->N_near_AND_gated_motifs[i];
    
    /*if not called by replay_mutation, store other info*/
    if(flag!=-1) 
    {    
        resident_info->avg_f=resident->avg_fitness;
        resident_info->f1=resident->fitness1;
        resident_info->f2=resident->fitness2;
        resident_info->se_avg_f=resident->SE_avg_fitness;
        resident_info->se_f1=resident->SE_avg_fitness;
        resident_info->se_f2=resident->SE_avg_fitness;
    
        if(flag==1) //if stores everthing
        {
            resident_info->step=evo_step;
            resident_info->n_mut_at_the_step=N_mutations_at_current_step;
            resident_info->n_tot_mut=N_tot_mutations;
            resident_info->n_hit_bound=mut_record->N_hit_bound;
            resident_info->selection_coefficient=selection_coefficient;

            resident_info->n_gene=resident->ngenes;
            resident_info->n_output_genes=resident->n_output_genes;
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
    mutant_info->se_f1=mutant->SE_avg_fitness;
    mutant_info->se_f2=mutant->SE_avg_fitness;
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
        fprintf(fp,"%d %d %d %c %d %d '%s' %d %a\n",
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
        for(i=0;i<17;i++)
            fprintf(fp,"%d ",resident_info[j].n_motifs[i]);    
        fprintf(fp,"\n");        
    }
    fflush(fp);
    fclose(fp); 
#if COUNT_NEAR_AND
    fp=fopen("N_near_AND_gated_motifs.txt","a+");
    for(j=0;j<output_counter;j++)
    {
        for(i=0;i<12;i++)    
            fprintf(fp,"%d ",resident_info[j].n_near_AND_gated_motifs[i]);
        fprintf(fp,"\n");   
    }
    fflush(fp);
    fclose(fp);    
#endif
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
            fprintf(fp,"%d %d %d %d %c %f %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d\n",
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
                    resident_info[i].n_act,
                    resident_info[i].n_rep);
        fflush(fp);
        fclose(fp);
    }
}