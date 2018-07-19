/* 
 * Simulator of yeast transcriptional regulatory network evolution
 * 
 * This file contains functions to initialize simulation with specified 
 * selection condition, and to output summary of genotypes and network structure.
 * 
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2018 Arizona Board of Regents (University of Arizona)
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

int MAXELEMENTS=100; 

static const float PROB_ACTIVATING=0.62;
static const float MEAN_PROTEIN_DECAY_RATE=-1.88;
static const float SD_PROTEIN_DECAY_RATE=0.561;
static const float MEAN_ACT_TO_INT_RATE=1.27;
static const float SD_ACT_TO_INT_RATE=0.226;
static const float MEAN_MRNA_DECAY_RATE=-1.49;
static const float SD_MRNA_DECAY_RATE=0.267;
static const float MEAN_PROTEIN_SYN_RATE=0.322;
static const float SD_PROTEIN_SYN_RATE=0.416;
const float MEAN_GENE_LENGTH=2.568; //log10(aa)
static const float SD_GENE_LENGTH=0.34;
static const float MIN_Kd=1.0e-9;
static const float MAX_Kd=1.0e-6;
static const float log_MIN_Kd=-9.0;
static const float log_MAX_Kd=-6.0;
static const float NS_Kd=1.0e-5;
const float KD2APP_KD=1.8e10;
static const float DEFAULT_UPDATE_INTERVAL=10.0; /*min*/
static const float MAX_TOLERABLE_CHANGE_IN_PROBABILITY_OF_BINDING=0.01;
static const float MIN_SELECTION_COEFFICIENT=1.0e-8;
/*Bounds*/
const float MAX_ACT_TO_INT_RATE=64.7;
const float MIN_ACT_TO_INT_RATE=0.59;
const float MAX_MRNA_DECAY=0.54;
const float MIN_MRNA_DECAY=7.5e-4;
const float MAX_PROTEIN_DECAY=0.69;
const float MIN_PROTEIN_DECAY=3.0e-6;
const float MAX_PROTEIN_SYN_RATE=61.4;
const float MIN_PROTEIN_SYN_RATE=4.5e-3;
const float MAX_KD=1.0e-5;
const float MIN_KD=0.0;
const int MAX_GENE_LENGTH=5000; //aa
const int MIN_GENE_LENGTH= 50; //aa

/*Number of genes*/
int N_TF_GENES=MAX_TF_GENES;
int N_EFFECTOR_GENES=MAX_EFFECTOR_GENES;

/******************************************************************************
 * 
 *                     Private function prototypes
 *
 *****************************************************************************/
static void initialize_sequence(char *, int, int, RngStream);

static void initialize_genotype_fixed(Genotype *, int, int, int, RngStream);

static void calc_avg_fitness(Genotype *, Selection *, int [MAX_GENES], float [MAX_PROTEINS], RngStream [N_THREADS], float *, float *); 

static void clone_genotype(Genotype *, Genotype *);

static float try_fixation(Genotype *, Genotype *, int, int, int *, RngStream);

static void summarize_binding_sites(Genotype *,int);

static void set_signal(CellState *, Test *, RngStream, int);

static void output_genotype(Genotype *);

static int evolve_N_steps(Genotype *, Genotype *,  Mutation *, Selection *, int *, int *, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS], int);

static void run_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_PROTEINS], int, int, RngStream, RngStream [N_THREADS]);

static void continue_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS]);

static void calc_fitness_stats(Genotype *, Selection *, float (*)[N_REPLICATES], float (*)[N_REPLICATES], int);

static void replay_mutations(Genotype *, Genotype *, Mutation *, FILE *, int);

static void find_motifs(Genotype *);

static void tidy_output_files(char*, char*);

static void print_motifs(Genotype *);

#if PERTURB
static void remove_edges_iteratively(Genotype *);

static void modify_topology(Genotype *, Genotype *);

static void add_binding_site(Genotype *, int);

static void remove_binding_sites(Genotype *, int);
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
        int replay_N_steps=0;       
        fscanf(fp,"%d",&replay_N_steps);
        fclose(fp);
        fp=fopen("sim_setup*","a+");
        fprintf(fp,"Continue simulation at step %d\n",replay_N_steps);
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
        fp=fopen(output_file,"w");
        fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step N_hit_bound accepted_mutn selection_coeff avg_fitness fitness1 fitness2 se_avg_fitness se_fitness1 se_fitness2 N_genes N_proteins N_activator N_repressor\n");
        fprintf(fp,"0 0 0 0 na na %.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d \n",  
                resident->avg_fitness,               
                resident->fitness1,
                resident->fitness2,
                sqrt(resident->sq_SE_avg_fitness),
                sqrt(resident->sq_SE_fitness1),
                sqrt(resident->sq_SE_fitness2),
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
void evolve_neutrally(Genotype *resident, Genotype *mutant, Mutation *mut_record, Selection *burn_in, Selection *selection, RngStream RS_main)
{
    int i;    
    FILE *fp;    
    /*Create title for output files*/
    fp=fopen(output_file,"a+");
    fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step N_hit_bound accepted_mutn selection_coeff avg_fitness fitness1 fitness2 se_avg_fitness se_fitness1 se_fitness2 N_genes N_proteins N_activator N_repressor\n");
    fprintf(fp,"0 0 0 na na 0.0 0.0 0.0 0.0 0.0 0.0 0 0 0 0 \n");
    fclose(fp); 
    /*BURN-IN*/              
    DUPLICATION=burn_in->temporary_DUPLICATION;                 
    SILENCING=burn_in->temporary_SILENCING;
    N_EFFECTOR_GENES=burn_in->temporary_N_effector_genes;
    N_TF_GENES=burn_in->temporary_N_tf_genes; 
    miu_ACT_TO_INT_RATE=burn_in->temporary_miu_ACT_TO_INT_RATE; 
    miu_Kd=burn_in->temporary_miu_Kd;       
    miu_protein_syn_rate=burn_in->temporary_miu_protein_syn_rate;     
    for(i=0;i<burn_in->MAX_STEPS;i++)
    {       
        clone_genotype(resident,mutant);
        mutate(mutant,RS_main,mut_record);   
        fp=fopen(mutation_file,"a+");
        fprintf(fp,"%c %d %d '%s' %d %a\n",
                mut_record->mut_type,    
                mut_record->which_gene,
                mut_record->which_nucleotide,
                mut_record->nuc_diff,
                mut_record->kinetic_type,
                mut_record->kinetic_diff);
        fclose(fp);        
        clone_genotype(mutant,resident);         
        calc_all_binding_sites(resident);
        find_motifs(resident); 
        print_motifs(resident);        
        /*output network topology every OUTPUT_INTERVAL steps*/ 
        if(i%OUTPUT_INTERVAL==0 && i!=0)
            summarize_binding_sites(resident,i);        
        /*output a summary of simulation every step*/
        output_genotype(resident);
    } 
    
    DUPLICATION=selection->temporary_DUPLICATION;                 
    SILENCING=selection->temporary_SILENCING;
    N_EFFECTOR_GENES=selection->temporary_N_effector_genes;
    N_TF_GENES=selection->temporary_N_tf_genes; 
    miu_ACT_TO_INT_RATE=selection->temporary_miu_ACT_TO_INT_RATE; 
    miu_Kd=selection->temporary_miu_Kd;       
    miu_protein_syn_rate=selection->temporary_miu_protein_syn_rate;      
    
    for(;i<selection->MAX_STEPS;i++)
    {
        clone_genotype(resident,mutant);
        mutate(mutant,RS_main,mut_record);   
        fp=fopen(mutation_file,"a+");
        fprintf(fp,"%c %d %d '%s' %d %a\n",
                mut_record->mut_type,    
                mut_record->which_gene,
                mut_record->which_nucleotide,
                mut_record->nuc_diff,
                mut_record->kinetic_type,
                mut_record->kinetic_diff);
        fclose(fp);        
        clone_genotype(mutant,resident);         
        calc_all_binding_sites(resident);
        find_motifs(resident); 
        print_motifs(resident);        
        /*output network topology every OUTPUT_INTERVAL steps*/ 
        if(i%OUTPUT_INTERVAL==0 && i!=0)
            summarize_binding_sites(resident,i);        
        /*output a summary of simulation every step*/
        output_genotype(resident);
    }
    print_mutatable_parameters(resident,1);    
}
#endif

#if PHENOTYPE
void show_phenotype(Genotype *resident, Genotype *mutant, Mutation *mut_record, Selection *selection, int init_mRNA[MAX_GENES], float init_protein[MAX_GENES], RngStream RS_parallel[N_THREADS])
{   
    FILE *fp;
    
    /*load mutation record*/
    fp=fopen("MUT.txt","r");    
    if(fp!=NULL)        
        printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
    else
    {
        printf("Loading mutation record failed! Quit program!");
#if !LOG_OFF
        LOG("Loading mutation record failed!");
#endif
        exit(-2);
    }
    
    /*replay mutations*/    
    replay_mutations(resident, mutant, mut_record, fp, selection->MAX_STEPS); 
    fclose(fp);
    
    /*output the evolved genotype*/
    calc_all_binding_sites(resident); 
    print_mutatable_parameters(resident,1);
    summarize_binding_sites(resident,selection->MAX_STEPS);   
    //exit(0);    
    
    /*create threads*/
    omp_set_num_threads(N_THREADS);  
    
    /*collection interval is 1 minute by default*/    
    selection->test1.t_development=90.1;
    selection->test2.t_development=90.1;    
    calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, NULL,NULL);    
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
    char buffer[600],char_buffer;
    int int_buffer,step;
    float float_buffer, mean_overall_fitness, mean_fitness1, mean_fitness2, se_overall_fitness, se_fitness1, se_fitness2;
    float fitness1[HI_RESOLUTION_RECALC][N_REPLICATES],fitness2[HI_RESOLUTION_RECALC][N_REPLICATES]; 
    FILE *file_mutation,*fitness_record,*f_aft_perturbation,*f_bf_perturbation;
    
    /*load mutation record*/
    file_mutation=fopen("all_mutations*.txt","r");    
    if(file_mutation!=NULL)        
        printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
    else
    {
        printf("Loading mutation record failed! Quit program!");
#if !LOG_OFF
        LOG("Loading mutation record failed!");
#endif
        exit(-2);
    } 
    
    /*skip first 2 rows of fitness_record*/
    fitness_record=fopen("evo_summary*","r");
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
        fscanf(file_mutation,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                    &(mut_record->which_gene),
                                                    &(mut_record->which_nucleotide), 
                                                    mut_record->nuc_diff,               
                                                    &(mut_record->kinetic_type),
                                                    &(mut_record->kinetic_diff));
        reproduce_mutate(resident,mut_record);        
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
            find_motifs(resident); 
#if DIRECT_REG
            if(resident->N_motifs[5]!=0 && resident->N_motifs[5]==resident->N_motifs[0])
#else          
    #if FORCE_SINGLE_FFL
            if(resident->N_motifs[23]!=0 && 
                resident->N_motifs[9]==0 && 
                resident->N_motifs[23]==resident->N_motifs[18] &&
                resident->N_motifs[27]==0 &&
                resident->N_motifs[36]==0 &&
                resident->N_motifs[38]==0) //force single ffl
    #elif FORCE_DIAMOND //networks contain only AND-gated double C1ffl
            if(resident->N_motifs[23]!=0 && 
                resident->N_motifs[9]==0 && 
                resident->N_motifs[23]==resident->N_motifs[18] && 
                resident->N_motifs[27]==0 &&
                resident->N_motifs[36]==0 &&
                resident->N_motifs[38]==0) 
    #else //DISABLE_AND_GATE. 
        #if WHICH_MOTIF==0 // disturb C1-FFL
            if(resident->N_motifs[14]!=0 && 
                resident->N_motifs[18]==0 && 
                resident->N_motifs[14]==resident->N_motifs[9] && 
                resident->N_motifs[27]==0 )
        #elif WHICH_MOTIF==1 // disturb FFL-in-diamond
            if(resident->N_motifs[23]!=0 && 
                resident->N_motifs[9]==0 && 
                resident->N_motifs[23]==resident->N_motifs[18] && 
                resident->N_motifs[27]==0 )
        #else // disturb diamond
            if(resident->N_motifs[32]!=0 && 
                resident->N_motifs[9]==0 && 
                resident->N_motifs[27]==resident->N_motifs[32] && 
                resident->N_motifs[18]==0 )
        #endif                
    #endif        
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
                        sqrt(resident->sq_SE_avg_fitness),
                        sqrt(resident->sq_SE_fitness1),
                        sqrt(resident->sq_SE_fitness2));
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
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype, int init_TF_genes, int init_N_act, int init_N_rep, int init_effector_genes, RngStream RS)
{ 
    int i,k;

    genotype->ngenes=init_effector_genes+N_SIGNAL_TF+init_TF_genes; /*including the signal genes and 1 selection gene*/
    genotype->ntfgenes=N_SIGNAL_TF+init_TF_genes; /*including the signal genes*/
    genotype->nproteins=genotype->ngenes;  /*at initialization, each protein is encoded by one copy of gene*/   
    genotype->nTF_families=genotype->nproteins-1;
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
    for(i=0;i<genotype->nTF_families;i++)
    {
        genotype->which_TF_family[i]=i;
        genotype->TF_family_pool[i][0][0]=1;
        genotype->TF_family_pool[i][1][0]=i;
    }
    initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*MAX_GENES, genotype->ngenes, RS);  // initialize cis-reg sequence
    initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*MAX_TF_GENES, genotype->ntfgenes, RS);    //initialize binding sequence of TFs    
    /* We now generate the complementary sequence of BS that are on the non-template strand.
     * The complementary sequence is used to search for BS that on the non-template strand.  
     * We also assume that all the TFs can work on both strands, but can induce expression in one direction.*/  
    for(i=0;i< genotype->ntfgenes;i++)
    {        
        for(k=0;k<TF_ELEMENT_LEN;k++)
        {
            switch (genotype->tf_seq[i][TF_ELEMENT_LEN-k-1])
            {
                case 'a': genotype->tf_seq_rc[i][k]='t'; break;
                case 't': genotype->tf_seq_rc[i][k]='a'; break;
                case 'c': genotype->tf_seq_rc[i][k]='g'; break;
                case 'g': genotype->tf_seq_rc[i][k]='c'; break;
            }
        }        
    }     
    initialize_genotype_fixed(genotype, init_N_act, init_N_rep, init_effector_genes, RS);     
    calc_all_binding_sites(genotype);
}



/*****************************************************************************
 * 
 *                           Private functions
 *
 ****************************************************************************/
static void initialize_sequence(char *Seq, int len, int num_elements, RngStream RS)
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
static void initialize_genotype_fixed(Genotype *genotype, int init_N_act, int init_N_rep, int init_effector_genes, RngStream RS)
{
    int i;
    /* the first N_SIGNAL_TF genes encode the sensor TFs. The concentration of a sensor TF
     * is determined by certain environmental signal*/
    genotype->total_loci_length=0.0;    
    for(i=N_SIGNAL_TF; i < genotype->ngenes; i++) 
    {  
        #if RANDOM_COOPERATION_LOGIC        
            genotype->min_act_to_transc[i]=RngStream_RandInt(RS,1,2); //if one activator is sufficient to induce expression, the gene is regualted by OR gate.
        #else
            genotype->min_N_activator_to_transc[i]=1; 
            genotype->min_N_activator_to_transc[genotype->ngenes-1]=2;
        #endif

        /* tf affinity */
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
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
        genotype->translation_rate[i] = pow(10.0,SD_PROTEIN_SYN_RATE*gasdev(RS)+MEAN_PROTEIN_SYN_RATE);  
        if(genotype->translation_rate[i]>MAX_PROTEIN_SYN_RATE)
            genotype->translation_rate[i]=MAX_PROTEIN_SYN_RATE;
        if(genotype->translation_rate[i]<MIN_PROTEIN_SYN_RATE)
            genotype->translation_rate[i]=MIN_PROTEIN_SYN_RATE;
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
    genotype->N_act=0;
    genotype->N_rep=0;    
    if(init_N_rep==-1 && init_N_act==-1) /*randomly generate activators and repressors*/
    {
        for(i=N_SIGNAL_TF;i<genotype->ntfgenes;i++)
        {   
            if (RngStream_RandU01(RS)<PROB_ACTIVATING) 
            {
                genotype->N_act++; 
                genotype->protein_identity[i] = ACTIVATOR;
            }
            else 
            {
                genotype->N_rep++;
                genotype->protein_identity[i]= REPRESSOR;
            }
        }
    }
    else
    {
        genotype->N_act=init_N_act;
        genotype->N_rep=init_N_rep;
        for(i=N_SIGNAL_TF;i<N_SIGNAL_TF+init_N_act;i++)            
            genotype->protein_identity[i]=ACTIVATOR;            
        for(i=N_SIGNAL_TF+init_N_act;i<genotype->ntfgenes;i++)
            genotype->protein_identity[i]=REPRESSOR;
    }
    /* parameterize sensor TF*/ 
    for(i=0;i<N_SIGNAL_TF;i++)
    {
        genotype->mRNA_decay_rate[i]=0.0; // we assume environmental signal toggles the state of sensor TF between active and inactive 
        genotype->protein_decay_rate[i]=0.0; // the concentration of sensor TF is constant.
        genotype->translation_rate[i]=0.0;
        genotype->active_to_intermediate_rate[i]=0.0; 
        genotype->protein_identity[i]=ACTIVATOR; /*make sensor TF an activator*/
        genotype->N_act++;        
        genotype->Kd[i]=pow(10.0,(log_MAX_Kd-log_MIN_Kd)*RngStream_RandU01(RS)+log_MIN_Kd); 
    }
#if RANDOMIZE_SIGNAL2
    #if N_SIGNAL_TF==2
        if(RngStream_RandU01(RS)<=0.5) // we assume there is a background "on" signal, which is sensor TF 0, in the network.
            genotype->protein_identity[1]=ACTIVATOR; // Other sensor TFs can be either activators or repressors.
        else
        {
            genotype->protein_identity[1]=REPRESSOR;
            genotype->N_act--;
            genotype->N_rep++;
        }
    #endif
#endif
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
                   
                    for(j=0;j<MAX_GENES;j++)
                    {
                        genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                        if(!genotype->all_binding_sites[j]) 
                        {  
#if !LOG_OFF
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
                if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
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

                        for(j=0;j<MAX_GENES;j++)
                        {
                            genotype->all_binding_sites[j] = realloc(genotype->all_binding_sites[j], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
                            if(!genotype->all_binding_sites[j]) 
                            {
#if !LOG_OFF
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
                    if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                }
            } 
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
    /* calculate max_unhindered_sites */
    /* max_unhindered_sites is maximum number of TFs that can bind to a cis-reg sequence at the same time*/
    /* We use it to faciliate the calculation of Pact and Prep. See calc_TF_dist_from_all_BS for its usage.*/
    int act_BS[MAXELEMENTS][2],rep_BS[MAXELEMENTS][2];
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
        for(gene_id=0;gene_id<MAX_GENES;gene_id++)
            genotype->all_binding_sites[gene_id]=realloc(genotype->all_binding_sites[gene_id], MAXELEMENTS*sizeof(AllTFBindingSites));
        genotype->N_allocated_elements=MAXELEMENTS;
    }
    for(gene_id=N_SIGNAL_TF;gene_id < genotype->ngenes;gene_id++)
    {        
        if(genotype->recalc_TFBS[gene_id]) /* do not calculate the binding sites if there's no mutation in the promoter or in TF binding seq*/
        {            
            calc_all_binding_sites_copy(genotype,gene_id);
            genotype->recalc_TFBS[gene_id]=NO;
        }
    }
}

/*
 * Set how the environmental signal should change
 */
static void set_signal(CellState *state, Test *test, RngStream RS, int thread_ID)
{
    float t=0.0;     
    char flag;   
    
#if EXTERNAL_SIGNAL
    int j;
    j=RngStream_RandInt(RS,0,99);
    test->external_signal=&(signal_profile_matrix[thread_ID][j][0]);
#else
    test->external_signal=NULL;             
#endif    
    
    if(test->external_signal==NULL)   
    {
        flag='o'; 
        state->protein_number[N_SIGNAL_TF-1]=test->signal_on_strength;    //always start with signal on
    #if N_SIGNAL_TF==2
        state->protein_number[0]=background_signal_strength;
    #endif      
        while(t<test->t_development)
        {
            if(flag=='o')
            {
                add_fixed_event(-1,t+test->t_signal_on,&(state->signal_off_head),&(state->signal_off_tail));
                flag='f';
                t=t+test->t_signal_on;                    
            }    
            else
            {
                add_fixed_event(-1,t+test->t_signal_off,&(state->signal_on_head),&(state->signal_on_tail));
                flag='o';
                t=t+test->t_signal_off;
            }
        } 
    }
    else
    {
        int time_point=1;
        state->protein_number[N_SIGNAL_TF-1]=test->external_signal[0];
        t=10.0;
        while(t<test->t_development)
        {
            add_fixed_event(time_point,t,&(state->change_signal_strength_head),&(state->change_signal_strength_tail));
            time_point++;
            t+=10.0;
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
    for(i=0; i< genotype_templet->ngenes;i++)
    {
        genotype_clone->which_cluster[i]=genotype_templet->which_cluster[i];            
        memcpy(&genotype_clone->cisreg_seq[i][0],&genotype_templet->cisreg_seq[i][0],CISREG_LEN*sizeof(char));                    
        genotype_clone->recalc_TFBS[i]=YES;                
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
    for(i=0;i<MAX_GENES;i++)
    {
        genotype_clone->which_protein[i]=NA;
        genotype_clone->min_N_activator_to_transc[i]=MAX_BINDING+1;   
    }    
    /*reset clone's tf_family_pool and which_tf_family*/
    for(i=0;i<MAX_PROTEINS;i++)
    {        
        for(j=0;j<MAX_PROTEINS;j++)
            genotype_clone->TF_family_pool[i][1][j]=NA;
        genotype_clone->TF_family_pool[i][0][0]=0;
        genotype_clone->which_TF_family[i]=NA;
    }
    /*copy from templet's tf_family_pool*/
    for(i=0;i<genotype_templet->nTF_families;i++)
    {
        genotype_clone->TF_family_pool[i][0][0]=genotype_templet->TF_family_pool[i][0][0];
        for(j=0;j<genotype_templet->TF_family_pool[i][0][0];j++)
            genotype_clone->TF_family_pool[i][1][j]=genotype_templet->TF_family_pool[i][1][j];
    }
    /*reset clone's protein_pool*/
    for(i=0;i<MAX_PROTEINS;i++)
    {            
        for(j=0;j<MAX_GENES;j++)
            genotype_clone->protein_pool[i][1][j]=NA;
        genotype_clone->protein_pool[i][0][0]=0;            
    }
    /*copy from templet's protein_pool and tf_family_pool*/
    for(i=0;i<genotype_templet->nproteins;i++)
    {            
        genotype_clone->which_TF_family[i]=genotype_templet->which_TF_family[i];
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
        genotype_clone->mRNA_decay_rate[i]=genotype_templet->mRNA_decay_rate[i];
        genotype_clone->protein_decay_rate[i]=genotype_templet->protein_decay_rate[i];
        genotype_clone->translation_rate[i]=genotype_templet->translation_rate[i];            
        genotype_clone->active_to_intermediate_rate[i]=genotype_templet->active_to_intermediate_rate[i];
        genotype_clone->which_protein[i]=genotype_templet->which_protein[i];
        genotype_clone->locus_length[i]=genotype_templet->locus_length[i];
        genotype_clone->min_N_activator_to_transc[i]=genotype_templet->min_N_activator_to_transc[i];   
    } 
    /* copy TF information*/
    for(i=0;i<MAX_PROTEINS;i++)
    {
        genotype_clone->protein_identity[i]=genotype_templet->protein_identity[i];
        genotype_clone->Kd[i]=genotype_templet->Kd[i];
    }    
    /* copy gene and protein numbers*/
    genotype_clone->ngenes=genotype_templet->ngenes;
    genotype_clone->ntfgenes=genotype_templet->ntfgenes;
    genotype_clone->nproteins=genotype_templet->nproteins;
    genotype_clone->nTF_families=genotype_templet->nTF_families;
    genotype_clone->N_act=genotype_templet->N_act;
    genotype_clone->N_rep=genotype_templet->N_rep;
    genotype_clone->total_loci_length=genotype_templet->total_loci_length;   
}

/**
 *Calculate the fintess of a given genotype.
 *Essentially calling do_single_timestep until tdevelopment and calculate 
 *average growth rate over tdevelopment.
 */
static void calc_avg_fitness(  Genotype *genotype,
                                    Selection *Selection,
                                    int init_mRNA[MAX_GENES],
                                    float init_protein_number[MAX_PROTEINS],
                                    RngStream RS_parallel[N_THREADS], 
                                    float GR1[N_REPLICATES],
                                    float GR2[N_REPLICATES])       
{   
    Phenotype timecourse1[N_REPLICATES], timecourse2[N_REPLICATES]; 
#if PHENOTYPE     
    int i,j;   
    /*alloc space and initialize values to 0.0*/
    for(i=0;i<N_REPLICATES;i++)
    {
        timecourse1[i].total_time_points=(int)Selection->test1.t_development;
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
        timecourse2[i].total_time_points=(int)Selection->test2.t_development;
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
    }        
#endif
    
    /*Making clones of a genotype, and have the clones run in parallel*/
    #pragma omp parallel num_threads(N_THREADS) 
    {
        int thread_ID=omp_get_thread_num();
//        int thread_ID=0;
        int i,j,k;
        int N_replicates_per_thread=N_REPLICATES/N_THREADS;  
        Genotype genotype_clone;
        CellState state_clone;
        GillespieRates rate_clone;
        int init_mRNA_clone[MAX_GENES]; 
        float init_protein_number_clone[MAX_GENES];
        float gr1[N_replicates_per_thread],gr2[N_replicates_per_thread];       
        int mRNA[genotype->ngenes];
        float protein[genotype->ngenes];         
        Test Test1, Test2;

        /*initialize the clone*/
        initialize_cache(&genotype_clone);
        
        /*clone genotype and initial mRNA and protein numbers*/
        #pragma omp critical
        { 
            genotype_clone.ngenes=genotype->ngenes;
            genotype_clone.ntfgenes=genotype->ntfgenes;
            genotype_clone.nproteins=genotype->nproteins;
            clone_genotype(genotype, &genotype_clone);              
            Test1.t_development=Selection->test1.t_development;
            Test1.signal_on_strength=Selection->test1.signal_on_strength;
            Test1.signal_off_strength=Selection->test1.signal_off_strength;
            Test1.t_signal_on=Selection->test1.t_signal_on;
            Test1.t_signal_off=Selection->test1.t_signal_off;
            Test1.initial_effect_of_effector=Selection->test1.initial_effect_of_effector;
            Test1.fixed_effector_effect=Selection->test1.fixed_effector_effect;
            Test2.t_development=Selection->test2.t_development;
            Test2.signal_on_strength=Selection->test2.signal_on_strength;
            Test2.signal_off_strength=Selection->test2.signal_off_strength;
            Test2.t_signal_on=Selection->test2.t_signal_on;
            Test2.t_signal_off=Selection->test2.t_signal_off;
            Test2.initial_effect_of_effector=Selection->test2.initial_effect_of_effector;
            Test2.fixed_effector_effect=Selection->test2.fixed_effector_effect;
            Test1.duration_of_burn_in_growth_rate=Selection->test1.duration_of_burn_in_growth_rate;
            Test2.duration_of_burn_in_growth_rate=Selection->test2.duration_of_burn_in_growth_rate;
            for(j=0; j < MAX_GENES; j++) 
            {  
                init_mRNA_clone[j] = init_mRNA[j];
                init_protein_number_clone[j] = init_protein_number[j];
            } 
        } 
        calc_all_binding_sites(&genotype_clone); 
        
        /*modify network */ 
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
        /* now calc fitness under the two tests*/
        /********************************************************************** 
         * 
         *                              TEST1 
         *
         *********************************************************************/
        for(i=0;i<N_replicates_per_thread;i++) /* env 1, usually a constant signal that matches env*/
        {
            /*initialize mRNA and protein numbers, and gene states etc.*/
            initialize_cell(&genotype_clone, &state_clone, &Test1, mRNA, protein);
            
            /*set how the signal should change during simulation*/
            set_signal(&state_clone, &Test1, RS_parallel[thread_ID], thread_ID);
            
            /*calcualte the rates of cellular activity based on the initial cellular state*/
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, &Test1, INITIALIZATION);             
#if PHENOTYPE
            timecourse1[thread_ID*N_replicates_per_thread+i].timepoint=0;
#endif
            /*run growth simulation until tdevelopment or encounter an error*/
            while(state_clone.t<Test1.t_development) 
                do_single_timestep(&genotype_clone, &state_clone, &rate_clone, &Test1, &(timecourse1[thread_ID*N_replicates_per_thread+i]), RS_parallel[thread_ID]);
                      
            /*calculate average growth rate*/
            gr1[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/(Test1.t_development-Test1.duration_of_burn_in_growth_rate); 
#if PHENOTYPE
            timecourse1[thread_ID*N_replicates_per_thread+i].timepoint++;
#endif          
            /*free linked tables*/
            free_fixedevent(&state_clone);           
        }   
        
        /********************************************************************** 
         * 
         *                              TEST2 
         *
         *********************************************************************/
        for(i=0;i<N_replicates_per_thread;i++) 
        {       
            initialize_cell(&genotype_clone, &state_clone, &Test2, mRNA, protein);
            set_signal(&state_clone, &Test2, RS_parallel[thread_ID], thread_ID);
            calc_all_rates(&genotype_clone, &state_clone, &rate_clone, &Test2, INITIALIZATION); 
#if PHENOTYPE
            timecourse2[thread_ID*N_replicates_per_thread+i].timepoint=0;
#endif            
            while(state_clone.t<Test2.t_development) 
                do_single_timestep(&genotype_clone, &state_clone, &rate_clone, &Test2, &(timecourse2[thread_ID*N_replicates_per_thread+i]), RS_parallel[thread_ID]);            
        
            gr2[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/(Test2.t_development-Test2.duration_of_burn_in_growth_rate);
#if PHENOTYPE
            timecourse2[thread_ID*N_replicates_per_thread+i].timepoint++;
#endif           
            free_fixedevent(&state_clone);            
        } 
        /*free linked tables*/
        for(j=0;j<MAX_GENES;j++)
            free(genotype_clone.all_binding_sites[j]);
#if !PHENOTYPE       
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
#endif
    }     
#if PHENOTYPE
    /*output timecourse*/
    int k;
    char filename[32];
    FILE *fp;   
    /*fitness: each row is a replicate*/
    fp=fopen("fitnessA","w");
    for(i=0;i<N_REPLICATES;i++)
    {
        for(j=0;j<timecourse1[i].total_time_points;j++)
            fprintf(fp,"%f ",timecourse1[i].instantaneous_fitness[j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
    fp=fopen("fitnessB","w");
    for(i=0;i<N_REPLICATES;i++)
    {
        for(j=0;j<timecourse2[i].total_time_points;j++)
            fprintf(fp,"%f ",timecourse2[i].instantaneous_fitness[j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
    /*proteint concentration: each protein has its own file, in which each row is a replicate*/    
    for(i=0;i<genotype->nproteins;i++)
    {
        snprintf(filename,sizeof(char)*32,"protein%i_A",i);
        fp=fopen(filename,"w");
        for(j=0;j<N_REPLICATES;j++)
        {
            for(k=0;k<timecourse1[j].total_time_points;k++)                
                fprintf(fp,"%f ",timecourse1[i].protein_concentration[k+i*timecourse1[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    for(i=0;i<genotype->nproteins;i++)
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
     for(i=0;i<genotype->ngenes;i++)
    {
        snprintf(filename,sizeof(char)*32,"gene%i_A",i);
        fp=fopen(filename,"w");
        for(j=0;j<N_REPLICATES;j++)
        {
            for(k=0;k<timecourse1[j].total_time_points;k++)                
                fprintf(fp,"%f ",timecourse1[j].protein_concentration[k+i*timecourse1[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    for(i=0;i<genotype->ngenes;i++)
    {
        snprintf(filename,sizeof(char)*32,"gene%i_B",i);
        fp=fopen(filename,"w");
        for(j=0;j<N_REPLICATES;j++)
        {
            for(k=0;k<timecourse2[j].total_time_points;k++)                
                fprintf(fp,"%f ",timecourse2[j].protein_concentration[k+i*timecourse2[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
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
        genotype->Kd[j]=-1.0;
        genotype->locus_length[j]=0;
        for(k=0;k<MAX_GENES;k++)        
            genotype->cisreg_cluster[j][k]=NA;
    }    
    for(j=0;j<MAX_GENES;j++)
        genotype->cisreg_cluster[MAX_GENES][j]=NA;
    /* initialize variables that applies to protein */
    for(j=0;j<MAX_PROTEINS;j++)
    {
        genotype->which_TF_family[j]=NA;
        genotype->protein_pool[j][0][0]=0;
        genotype->TF_family_pool[j][0][0]=0;
        for(k=0;k<MAX_GENES;k++)        
            genotype->protein_pool[j][1][k]=NA; 
        for(k=0;k<MAX_PROTEINS;k++) 
            genotype->TF_family_pool[j][1][k]=NA;
        genotype->protein_identity[j]=NON_TF;
    }
    /* alloc space for binding sites*/
    genotype->N_allocated_elements=MAXELEMENTS;
    for(j=0;j<MAX_GENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(MAXELEMENTS*sizeof(AllTFBindingSites));
        if (!(genotype->all_binding_sites[j])) 
        {  
#if !LOG_OFF
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

/**
 * Given the fitness of the resident and a mutant, decide whether the mutant can replace the resident
 */
static float try_fixation(Genotype *resident, Genotype *mutant, int N_measurement_resident, int N_measurement_mutant, int *fixation, RngStream RS)
{ 
    float s;
    s=(mutant->avg_fitness-resident->avg_fitness)/fabs(resident->avg_fitness);
    if(s>=MIN_SELECTION_COEFFICIENT)
        *fixation=1;
    else          
        *fixation=0;     
    return s;
}

static void replay_mutations(Genotype *resident, Genotype *mutant, Mutation *mut_record, FILE *file_mutation, int replay_N_steps)
{
    int i;
    remove("proportion_c1ffl.txt"); 
    remove("summary_BS.txt");  
    
    calc_all_binding_sites(resident);
    summarize_binding_sites(resident,0); 
    find_motifs(resident);
    print_motifs(resident);   
    
    for(i=0;i<replay_N_steps;i++)
    {        
        clone_genotype(resident,mutant);
        fscanf(file_mutation,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                    &(mut_record->which_gene),
                                                    &(mut_record->which_nucleotide), 
                                                    mut_record->nuc_diff,               
                                                    &(mut_record->kinetic_type),
                                                    &(mut_record->kinetic_diff));
        reproduce_mutate(mutant,mut_record);        
        clone_genotype(mutant,resident); 
        calc_all_binding_sites(resident);
        if(i%OUTPUT_INTERVAL==0)
            summarize_binding_sites(resident,i); 
        find_motifs(resident);
        print_motifs(resident);   
    }
    printf("Reproduce mutations successfully!\n");
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
    first_step=init_step;
    N_tot_trials=init_N_tot_mutations; 
 
    /* first, run burn-in */
    if(burn_in->MAX_STEPS!=0)
    {
        flag_burn_in=1; 
        DUPLICATION=burn_in->temporary_DUPLICATION;                 
        SILENCING=burn_in->temporary_SILENCING;
        N_EFFECTOR_GENES=burn_in->temporary_N_effector_genes;
        N_TF_GENES=burn_in->temporary_N_tf_genes; 
        miu_ACT_TO_INT_RATE=burn_in->temporary_miu_ACT_TO_INT_RATE; 
        miu_Kd=burn_in->temporary_miu_Kd;       
        miu_protein_syn_rate=burn_in->temporary_miu_protein_syn_rate; 
        float fitness1[HI_RESOLUTION_RECALC][N_REPLICATES],fitness2[HI_RESOLUTION_RECALC][N_REPLICATES];   
        
        evolve_N_steps( resident, 
                        mutant,
                        mut_record, 
                        burn_in,
                        &first_step,                     
                        &N_tot_trials, 
                        init_mRNA, 
                        init_protein,
                        RS_main,
                        RS_parallel,
                        flag_burn_in);    
        
        /*Calculate fitness of the current genotype under post burn_in condition*/
        /*The "if" is always true when the simulation is run from the beginning,
         *i.e. when init_step=0. But when continuing a simulation from a saving  
         *point after the burn-in, the "if" is always false*/
        if(init_step==burn_in->MAX_STEPS)
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
            output_genotype(resident);
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
            fprintf(fp,"%d %d %a %a %a %a %a %a\n",N_tot_trials, 
                                                mut_record->N_hit_bound,
                                                resident->avg_fitness,                                                
                                                resident->fitness1,
                                                resident->fitness2,
                                                resident->sq_SE_avg_fitness,
                                                resident->sq_SE_fitness1,
                                                resident->sq_SE_fitness2);
            fclose(fp);        
            /* marks the last step at which all state of the program has been output*/
            fp=fopen("saving_point.txt","w");
            fprintf(fp,"%d\n",burn_in->MAX_STEPS);
            fclose(fp);
        }
    }    
    
    /* post-burn-in simulations*/
    flag_burn_in=0;    
    DUPLICATION=selection->temporary_DUPLICATION;                 
    SILENCING=selection->temporary_SILENCING;
    N_EFFECTOR_GENES=selection->temporary_N_effector_genes;
    N_TF_GENES=selection->temporary_N_tf_genes; 
    miu_ACT_TO_INT_RATE=selection->temporary_miu_ACT_TO_INT_RATE; 
    miu_Kd=selection->temporary_miu_Kd;       
    miu_protein_syn_rate=selection->temporary_miu_protein_syn_rate; 
    
    evolve_N_steps(resident, 
                    mutant,
                    mut_record, 
                    selection,
                    &first_step,                   
                    &N_tot_trials,   
                    init_mRNA,
                    init_protein,
                    RS_main,
                    RS_parallel,
                    flag_burn_in); 
    calc_all_binding_sites(resident);
    summarize_binding_sites(resident,selection->MAX_STEPS); /*snapshot of the final distribution binding sites */
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
    tidy_output_files(output_file,mutation_file);
    
    /* set genotype based on previous steps*/
    fp=fopen(mutation_file,"r");
    if(fp!=NULL)
        replay_mutations(resident, mutant, mut_record, fp, replay_N_steps);
    else
    {   
#if !LOG_OFF
        LOG("cannot open mutation_file\n"); 
#endif
        exit(-2);
    }
    fclose(fp);

    /* load random number seeds*/
    fp=fopen("RngSeeds.txt","r");
    if(fp!=NULL)
    {
        for(i=0;i<replay_N_steps/OUTPUT_INTERVAL;i++)
        {
            for(j=0;j<N_THREADS;j++)        
            {
                fscanf(fp,"%lu %lu %lu %lu %lu %lu ", &(rng_seeds[j][0]),
                                                        &(rng_seeds[j][1]),
                                                        &(rng_seeds[j][2]),
                                                        &(rng_seeds[j][3]),
                                                        &(rng_seeds[j][4]),
                                                        &(rng_seeds[j][5]));
            }
            fscanf(fp,"%lu %lu %lu %lu %lu %lu \n", &(rng_seeds[N_THREADS][0]),
                                                    &(rng_seeds[N_THREADS][1]),
                                                    &(rng_seeds[N_THREADS][2]),
                                                    &(rng_seeds[N_THREADS][3]),
                                                    &(rng_seeds[N_THREADS][4]),
                                                    &(rng_seeds[N_THREADS][5]));
        }
    }
    else
    {   
#if !LOG_OFF
        LOG("cannot open RngSeeds.txt\n");     
#endif
        exit(-2);
    }
    fclose(fp);
    RngStream_SetSeed(RS_main,rng_seeds[0]);
    for(i=0;i<N_THREADS;i++)
        RngStream_SetSeed(RS_parallel[i],rng_seeds[i+1]);
    
    /* load fitness,N_tot_mutations,N_hit_boundary*/
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
                                            &(resident->sq_SE_avg_fitness),
                                            &(resident->sq_SE_fitness1),
                                            &(resident->sq_SE_fitness2));
    }
    else
    {   
#if !LOG_OFF
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

static void calc_fitness_stats( Genotype *genotype,
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
            genotype->fitness_measurement[counter]=selection->test1_weight*f1[i][j]+selection->test2_weight*f2[i][j];
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
            sum_sq_diff_mean_f+=pow(diff_f1*selection->test1_weight+diff_f2*selection->test2_weight,2.0);
        }
    }
    sq_SE_f1=sum_sq_diff_f1/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));
    sq_SE_f2=sum_sq_diff_f2/(N_recalc_fitness*N_REPLICATES*(N_recalc_fitness*N_REPLICATES-1));     
    genotype->fitness1=avg_f1;
    genotype->fitness2=avg_f2;
    genotype->sq_SE_fitness1=sq_SE_f1*N_recalc_fitness*N_REPLICATES;
    genotype->sq_SE_fitness2=sq_SE_f2*N_recalc_fitness*N_REPLICATES;     
    genotype->avg_fitness=selection->test1_weight*avg_f1+selection->test2_weight*avg_f2;
    genotype->sq_SE_avg_fitness=sum_sq_diff_mean_f/(N_recalc_fitness*N_REPLICATES-1)/(N_recalc_fitness*N_REPLICATES); 
}

static int evolve_N_steps(  Genotype *resident, 
                            Genotype *mutant,
                            Mutation *mut_record, 
                            Selection *selection,
                            int *init_step,                       
                            int *N_tot_trials,        
                            int init_mRNA[MAX_GENES],   
                            float init_protein[MAX_PROTEINS],
                            RngStream RS_main,
                            RngStream RS_parallel[N_THREADS],
                            int flag_burn_in)
{
    int i,j;
    int fixation;
    int N_trials;
    float fitness1[HI_RESOLUTION_RECALC][N_REPLICATES],fitness2[HI_RESOLUTION_RECALC][N_REPLICATES];
    float score;  
    FILE *fp;
 
    for(i=(*init_step);i<=selection->MAX_STEPS;i++)
    {             
        fixation=0;      
        N_trials=0;
      
        while(1) /*try mutations until one replaces the current genotype*/
        {	
            N_trials++;
            (*N_tot_trials)++;
            if(N_trials>MAX_TRIALS) /*Tried too many mutation in one step.*/
            {
                fp=fopen(output_file,"a+");                              
                fprintf(fp,"Tried %d mutations, yet none could fix\n",MAX_TRIALS);
                fclose(fp); 
                summarize_binding_sites(resident,i);
                return -1;
            }
            /*do mutation on a copy of the current genotype*/
            clone_genotype(resident,mutant); 
            mutate(mutant,RS_main,mut_record);
#if OUTPUT_MUTANT_DETAILS
            /*record every mutation*/
            fp=fopen("all_mutations.txt","a+");
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
            calc_all_binding_sites(mutant);           
            MAXELEMENTS=mutant->N_allocated_elements;
            /*calculate the fitness of the mutant at low resolution*/
            calc_avg_fitness(mutant, selection, init_mRNA, init_protein, RS_parallel, fitness1[0], fitness2[0]);
            calc_fitness_stats(mutant, selection, &(fitness1[0]), &(fitness2[0]), 1); // calc fitness at low resolution
#if OUTPUT_MUTANT_DETAILS
            /*record low-resolution fitness*/
            fp=fopen("fitness_all_mutants.txt","a+");
            fprintf(fp,"%.10f %.10f\n", mutant->avg_fitness, mutant->sq_SE_avg_fitness);
            fclose(fp);     
#endif
            /*Can the mutant replace the current genotype?*/
            score=try_fixation(resident, mutant, N_REPLICATES*HI_RESOLUTION_RECALC, N_REPLICATES, &fixation, RS_main); 
            /*If yes, record relevant info*/
            if(fixation==1)
            {                    
                fp=fopen(output_file,"a+");
                fprintf(fp,"%d %d %d %d %c %f ",i, *N_tot_trials, N_trials,mut_record->N_hit_bound,mut_record->mut_type,score);
                fclose(fp);
                fp=fopen(mutation_file,"a+");
                fprintf(fp,"%c %d %d '%s' %d %a\n", mut_record->mut_type,    
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
        /*calculate the number of c1-ffls every step*/
        find_motifs(resident); 
        print_motifs(resident); 
        /*output network topology every OUTPUT_INTERVAL steps*/
        if(i%OUTPUT_INTERVAL==0 && i!=0) 
            summarize_binding_sites(resident,i);        
        /*output a summary of simulation every step*/
        if(!(i==selection->MAX_STEPS && flag_burn_in))
            output_genotype(resident);        
        /* output rng seeds*/
#if OUTPUT_RNG_SEEDS
        unsigned long seeds[6];
        if(!(i==selection->MAX_STEPS && flag_burn_in) && i%SAVING_INTERVAL==0)
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
            fclose(fp);
        }
#endif
        /*output precise hi-resolution fitness*/
        if(!(i==selection->MAX_STEPS && flag_burn_in))
        {
            fp=fopen("precise_fitness.txt","a+");
            fprintf(fp,"%d %d %a %a %a %a %a %a\n",*N_tot_trials, 
                                                mut_record->N_hit_bound,
                                                resident->avg_fitness,                                                
                                                resident->fitness1,
                                                resident->fitness2,
                                                resident->sq_SE_avg_fitness,
                                                resident->sq_SE_fitness1,
                                                resident->sq_SE_fitness2);
            fclose(fp);        
            /* marks the last step at which all state of the program has been output*/
            if(i%SAVING_INTERVAL==0)
            {
                fp=fopen("saving_point.txt","w");
                fprintf(fp,"%d\n",i);
                fclose(fp);
            }
        }
    } 
    *init_step=i;
    return 0;
}

static void output_genotype(Genotype *genotype)
{
    FILE *OUTPUT;
    OUTPUT=fopen(output_file,"a+");
    fprintf(OUTPUT,"%.10f %.10f %.10f %.10f %.10f %.10f %d %d %d %d \n",  
            genotype->avg_fitness,            
            genotype->fitness1,
            genotype->fitness2,
            sqrt(genotype->sq_SE_avg_fitness),
            sqrt(genotype->sq_SE_fitness1),
            sqrt(genotype->sq_SE_fitness2),
            genotype->ngenes,
            genotype->nproteins,
            genotype->N_act,
            genotype->N_rep);
    fclose(OUTPUT);   
}

static void print_motifs(Genotype *genotype)
{
    FILE *fp; 
    int i;
    fp=fopen("N_motifs.txt","a+");
    for(i=0;i<39;i++)
        fprintf(fp,"%d ",genotype->N_motifs[i]);    
    fprintf(fp,"\n");
    fclose(fp); 
#if COUNT_NEAR_AND
    fp=fopen("N_near_AND_gate_motifs.txt","a+");
    for(i=0;i<9;i++)    
        fprintf(fp,"%d ",genotype->N_near_AND_gated_motifs[i]);
    fprintf(fp,"\n");    
    fclose(fp);
#endif
}

static void summarize_binding_sites(Genotype *genotype, int step_i)
{
    FILE *OUTPUT1;
    int i,j;
    int table[MAX_GENES][MAX_GENES];
        
    for(i=0;i<genotype->ngenes;i++)
    {
        for(j=0;j<genotype->ngenes;j++)
            table[i][j]=0; 
    }
    
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)        
    {
        for(j=0;j<genotype->binding_sites_num[i];j++)
        {            
            table[i][genotype->all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
        }    
    }
    
    /*Output all binding sites*/ 
    OUTPUT1=fopen("networks.txt","a+");
    fprintf(OUTPUT1,"step %d\n",step_i);
    fprintf(OUTPUT1,"Gene     ");    
    for(i=0;i<genotype->nproteins-1;i++)
    {
        if(genotype->protein_identity[i]==ACTIVATOR)        
            fprintf(OUTPUT1," A%d ",i);      
        if(genotype->protein_identity[i]==REPRESSOR)        
            fprintf(OUTPUT1," R%d ",i);         
    }
    fprintf(OUTPUT1,"which_protein ");
    fprintf(OUTPUT1,"AND_gate_capable\n");
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
    {
        if(i<10)
            fprintf(OUTPUT1,"%d        ",i);
        else
            fprintf(OUTPUT1,"%d       ",i);
        
        for(j=0;j<genotype->nproteins-1;j++)
        {
            if(table[i][j]<10)
                fprintf(OUTPUT1," %d  ",table[i][j]);
            else
                fprintf(OUTPUT1," %d ",table[i][j]);
        }
        if(genotype->which_protein[i]==genotype->nproteins-1)
            fprintf(OUTPUT1,"      E  ");
        else
        {
            if(genotype->protein_identity[genotype->which_protein[i]]==ACTIVATOR)
                fprintf(OUTPUT1,"      A"); 
            else
                fprintf(OUTPUT1,"      R");
            if(genotype->which_protein[i]<10)
                fprintf(OUTPUT1,"%d ",genotype->which_protein[i]);
            else
                fprintf(OUTPUT1,"%d",genotype->which_protein[i]);
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

/*find subtypes of C1-FFLs, FFL-in-diamond, and diamonds*/
static void find_motifs(Genotype *genotype)
{
    int i,j,k,l,cluster_size;
    int found_bs;
    int gene_id,gene_id_copy,site_id,protein_id,N_copies,N_activators;
    int master_TF,aux_TF;
    int activators[MAX_PROTEINS];
    int copies_reg_by_env[MAX_GENES],copies_not_reg_by_env[MAX_GENES],N_copies_reg_by_env,N_copies_not_reg_by_env;
    int hindrance[MAX_PROTEINS][MAX_PROTEINS];
    int pos_binding_sites_of_j[MAXELEMENTS],pos_binding_sites_of_k[MAXELEMENTS],N_binding_sites_of_j,N_binding_sites_of_k; 
    int N_motifs;  
    int N_strong_BS[MAX_PROTEINS][2][20];
    float r_master_TF_deg, r_aux_TF_deg;
    
    /*record which genes and TFs form the motifs of interest. We use this information to perturb motifs*/
#if PERTURB
    for(i=0;i<MAX_GENES;i++)
    {
        genotype->gene_in_core_C1ffl[i]=0;
        for(j=0;j<MAX_PROTEINS;j++)
            genotype->TF_in_core_C1ffl[i][j]=0;
    }
#endif
    
    /*look for motifs only if the genotype contains non-signal activators*/
    if(genotype->N_act>N_SIGNAL_TF)
    {        
        for(i=0;i<39;i++)
            genotype->N_motifs[i]=0; 
#if COUNT_NEAR_AND
        for(i=0;i<9;i++)
            genotype->N_near_AND_gated_motifs[i]=0;
#endif
        /*loop through each cisreg_cluster. Genes in a cluster have the same binding sites*/
        i=0;   
        while(genotype->cisreg_cluster[i][0]!=NA) //not an empty cluster
        {              
            /*get a gene in the current cis-reg cluster*/
            gene_id=genotype->cisreg_cluster[i][0];
            gene_id_copy=gene_id;
            if(genotype->which_protein[gene_id]==genotype->nproteins-1) // is a effector gene
            {                 
                /*find how many genes in this cluster*/
                cluster_size=0;
                while(genotype->cisreg_cluster[i][cluster_size]!=NA)
                    cluster_size++;        
                /*reset table for recording activators that regulate effector gene*/
                for(j=0;j<MAX_PROTEINS;j++)
                    activators[j]=0;
#if COUNT_NEAR_AND
                /*reset table for recording strength of interaction*/
                for(j=0;j<MAX_PROTEINS;j++)
                {                   
                    N_strong_BS[j][0][0]=0;
                    for(k=0;k<20;k++)
                        N_strong_BS[j][1][k]=0;
                }
#endif
                /*scan binding sites for tfs that regulate gene_id*/
                for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                {
                    protein_id=genotype->all_binding_sites[gene_id][j].tf_id;
                    if(genotype->protein_identity[protein_id]==ACTIVATOR)
                    {
                        if(protein_id==N_SIGNAL_TF-1)
                        {
                            if(genotype->all_binding_sites[gene_id][j].mis_match<=CUT_OFF_MISMATCH_SIG2EFFECTOR)
                                activators[protein_id]=1;
                        }
                        else
                        {
                            if(genotype->all_binding_sites[gene_id][j].mis_match<=CUT_OFF_MISMATCH_TF2EFFECTOR)
                                activators[protein_id]=1;
                        }                        
                    }    
                }
                /* move non-zeros entries in activators to the front. */
                k=0; //marks the entry to which a record is copied
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
                        N_activators++; //total number of activators that regulate effector
                    }    
                    j++;
                } 
                
                /*build a table to show hindrance between binding sites on effector gene*/
                /*if every binding site of j can hinder all the binding site of k, denote hindrance[j][k]=1. Otherwise 0*/
                /*hindrance[j][j]=1 means at most one TF j can bind to the effector gene */
                /*A strict AND gate between j and K should have H[j][j]=H[k][k]=1, and H[j][k]=0*/
                /*set hindrance to 1 to begin with*/
                for(j=0;j<MAX_PROTEINS;j++)
                {
                    for(l=0;l<MAX_PROTEINS;l++)
                        hindrance[j][l]=1;
                }
                for(j=0;j<N_activators;j++)
                {
                    /*reset position table*/
                    for(l=0;l<MAXELEMENTS;l++)
                        pos_binding_sites_of_j[l]=-CISREG_LEN; 
                    /*list the positions of all the binding sites of activator j*/
                    N_binding_sites_of_j=0;
                    gene_id=genotype->cisreg_cluster[i][0];
                    for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                    {
                        if(genotype->all_binding_sites[gene_id][site_id].tf_id==activators[j])
                        {
                            if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_SIG2EFFECTOR)
                                {
                                    pos_binding_sites_of_j[N_binding_sites_of_j]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                    N_binding_sites_of_j++;
                                    if(genotype->all_binding_sites[gene_id][site_id].mis_match<(TF_ELEMENT_LEN-NMIN))
                                    {                                
                                        N_strong_BS[activators[j]][1][N_strong_BS[activators[j]][0][0]]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                        N_strong_BS[activators[j]][0][0]++;
                                    }
                                }
                            }
                            else
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF2EFFECTOR)
                                {
                                    pos_binding_sites_of_j[N_binding_sites_of_j]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                    N_binding_sites_of_j++;
                                    if(genotype->all_binding_sites[gene_id][site_id].mis_match<(TF_ELEMENT_LEN-NMIN))
                                    {                                
                                        N_strong_BS[activators[j]][1][N_strong_BS[activators[j]][0][0]]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                        N_strong_BS[activators[j]][0][0]++;
                                    }
                                }
                            }                                
                        }
                    } 
                    if(N_binding_sites_of_j==1 && genotype->min_N_activator_to_transc[gene_id_copy]==2) // there is only one binding site                     
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
                    
                    /*list the positions of all the binding sites of activator j+1...*/                  
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
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_SIG2EFFECTOR)
                                    {
                                        pos_binding_sites_of_k[N_binding_sites_of_k]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                        N_binding_sites_of_k++;
                                    }
                                }
                                else
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF2EFFECTOR)   
                                    {
                                        pos_binding_sites_of_k[N_binding_sites_of_k]=genotype->all_binding_sites[gene_id][site_id].BS_pos;
                                        N_binding_sites_of_k++;
                                    }
                                }
                            } 
                        }  
                        /*determine the relation between j and j+1, ...*/
                        if(abs(pos_binding_sites_of_j[N_binding_sites_of_j-1]-pos_binding_sites_of_k[0])>= TF_ELEMENT_LEN+2*HIND_LENGTH ||
                            abs(pos_binding_sites_of_k[N_binding_sites_of_k-1]-pos_binding_sites_of_j[0])>= TF_ELEMENT_LEN+2*HIND_LENGTH)
                        {
                            hindrance[activators[j]][activators[k]]=0; 
                            hindrance[activators[k]][activators[j]]=0;
                        }
                    }
                }
                
                /*make lists of gene copies that are regulated by environmental signal and of those that are not*/
                /*reset counters of tf genes*/
                N_copies_reg_by_env=0;
                N_copies_not_reg_by_env=0;
                for(j=0;j<MAX_GENES;j++)
                {
                    copies_reg_by_env[j]=-1;
                    copies_not_reg_by_env[j]=-1;
                }
                for(j=0;j<N_activators;j++)
                {
                    if(activators[j]>=N_SIGNAL_TF)
                    {  
                        N_copies=genotype->protein_pool[activators[j]][0][0]; // the number of genes encoding j
                        for(k=0;k<N_copies;k++)
                        {
                            gene_id=genotype->protein_pool[activators[j]][1][k];
                            /*reset flag*/
                            found_bs=0;
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1 && genotype->all_binding_sites[gene_id][site_id].tf_id<=CUT_OFF_MISMATCH_SIGNAL_TO_TF)
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

                /*******************************count motifs ************************/                
#if DIRECT_REG  
                /*count c1-ffls formed by signal and one signal-regulated copy */
                if(activators[0]==N_SIGNAL_TF-1 || activators[1]==N_SIGNAL_TF-1)
                {
                    N_motifs=cluster_size*N_copies_reg_by_env;
                    genotype->N_motifs[0]+=N_motifs;           

                    for(j=0;j<N_copies_reg_by_env;j++)
                    {
                        protein_id=genotype->which_protein[copies_reg_by_env[j]];
                        if(hindrance[N_SIGNAL_TF-1][protein_id]) 
                        {
                            if(hindrance[N_SIGNAL_TF-1][N_SIGNAL_TF-1]) // signal alone cannot activate effector
                            {
                                if(hindrance[protein_id][protein_id])
                                    genotype->N_motifs[1]+=cluster_size; // effectively no regulation
                                else
                                    genotype->N_motifs[2]+=cluster_size; // I3
                            }   
                            else
                            {
                                if(hindrance[protein_id][protein_id])
                                    genotype->N_motifs[3]+=cluster_size; //I1
                                else
                                    genotype->N_motifs[4]+=cluster_size; //impossible
                            }
                        }
                        else
                        {
                            if(hindrance[N_SIGNAL_TF-1][N_SIGNAL_TF-1])
                            {
                                if(hindrance[protein_id][protein_id])
                                {
                                    genotype->N_motifs[5]+=cluster_size; // AND-gated 
#if DISABLE_AND_GATE
                                    genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if !FORCE_MASTER_CONTROLLED
                                    genotype->TF_in_core_C1ffl[gene_id_copy][protein_id]=1;
#endif
#elif FORCE_DIAMOND
                                    genotype->gene_in_core_C1ffl[copies_reg_by_env[j]]=1;                                        
#endif
                                }
                                else
                                    genotype->N_motifs[6]+=cluster_size; // aux. tf controlled
                            }   
                            else
                            {
                                if(hindrance[protein_id][protein_id])
                                    genotype->N_motifs[7]+=cluster_size; // signal controll
                                else
                                    genotype->N_motifs[8]+=cluster_size; // OR-gated
                            }
                        }
                    }
                } 
#else
                /*count c1-ffl formed by one env-regulated copy and one unregulated copy*/
                for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=0;k<N_copies_not_reg_by_env;k++)
                    {
                        /*j and k must encode TFs of different families*/
                        if(genotype->which_TF_family[genotype->which_protein[copies_reg_by_env[j]]]!=genotype->which_TF_family[genotype->which_protein[copies_not_reg_by_env[k]]])
                        {
                            /*search bs of k on j*/
                            gene_id=copies_reg_by_env[j];
                            found_bs=0;
                            protein_id=genotype->which_protein[copies_not_reg_by_env[k]];
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
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
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(found_bs) // k is regulated by j, then we found an isolated C1FFL                           
                                {    
                                    genotype->N_motifs[9]+=cluster_size; //any logic
                                    master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    aux_TF=genotype->which_protein[copies_not_reg_by_env[k]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[10]+=cluster_size; //no regulation
                                            else  
                                                genotype->N_motifs[11]+=cluster_size; //emergent I3
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[12]+=cluster_size; //emergent I1
                                            else  
                                                genotype->N_motifs[13]+=cluster_size; //no possible
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[14]+=cluster_size; //AND-gated 
#if DISABLE_AND_GATE                                   
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                                             
#endif
                                            }
                                            else
                                            {
                                                genotype->N_motifs[15]+=cluster_size; //aux. TF controlled
#if COUNT_NEAR_AND
                                                /*count near-AND-gated motifs*/
                                                if(N_strong_BS[aux_TF][0][0]!=0 && N_strong_BS[master_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                                if(N_strong_BS[master_TF][0][0]==1 || N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                    genotype->N_near_AND_gated_motifs[0]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif  
                                            }
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[16]+=cluster_size; //master TF controlled
#if COUNT_NEAR_AND
                                                /*count near-AND-gated motifs*/
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                         	if(N_strong_BS[aux_TF][0][0]==1 || N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                genotype->N_near_AND_gated_motifs[1]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif  
                                            }
                                            else 
                                            {
                                                genotype->N_motifs[17]+=cluster_size; //OR-gated
#if COUNT_NEAR_AND
                                                /*count near-AND-gated motifs*/
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                        {
                                                            if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                                abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                                )
                                                                genotype->N_near_AND_gated_motifs[2]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif                                                
                                            }
                                        }                                     
                                    } 
                                }
                            }
                            else //j is regulated by k
                            {
                                /*search bs of j on k*/
                                gene_id=copies_not_reg_by_env[k];
                                found_bs=0;
                                protein_id=genotype->which_protein[copies_reg_by_env[j]];
                                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                                {
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(found_bs)
                                    genotype->N_motifs[36]+=cluster_size; // strongly-connected
                                else
                                    genotype->N_motifs[37]+=cluster_size; // ffl w/o source
                            }
                        }
                    }
                }
                /*count motifs formed by two env-regulated copies*/
                /*also count parallel structures*/
                for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=j+1;k<N_copies_reg_by_env;k++)
                    {
                        /*k and j must encode different TFs*/
                        if(k!=j && genotype->which_TF_family[genotype->which_protein[copies_reg_by_env[j]]]!=genotype->which_TF_family[genotype->which_protein[copies_reg_by_env[k]]])
                        {                    
                            /*search bs of k on j*/
                            gene_id=copies_reg_by_env[j];
                            found_bs=0;
                            protein_id=genotype->which_protein[copies_reg_by_env[k]];
                            for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                            {
                                if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
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
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(found_bs)  // k is regulated by j, then we've found a FFL-in-diamond                          
                                {                                    
                                    genotype->N_motifs[18]+=cluster_size; //any logic
                                    master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    aux_TF=genotype->which_protein[copies_reg_by_env[k]];
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[19]+=cluster_size; //no regulation
                                            else  
                                                genotype->N_motifs[20]+=cluster_size; //emergent I3
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[21]+=cluster_size; //emergent I1
                                            else  
                                                genotype->N_motifs[22]+=cluster_size; //impossible
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[23]+=cluster_size; //AND-gated
#if DISABLE_AND_GATE                                        
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                                        
#elif FORCE_DIAMOND                                           
                                                genotype->gene_in_core_C1ffl[copies_reg_by_env[k]]=1;
                                                genotype->TF_in_core_C1ffl[copies_reg_by_env[k]][genotype->which_TF_family[master_TF]]=1;                                           
#elif FORCE_SINGLE_FFL                                           
                                                genotype->gene_in_core_C1ffl[copies_not_reg_by_env[k]]=1;  
#endif 
                                            }
                                            else 
                                            {
                                                genotype->N_motifs[24]+=cluster_size; //aux. TF controlled
#if COUNT_NEAR_AND
                                                if(N_strong_BS[aux_TF][0][0]!=0 && N_strong_BS[master_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                                if(N_strong_BS[master_TF][0][0]==1 || N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                    genotype->N_near_AND_gated_motifs[3]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif
                                            }
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[25]+=cluster_size; //master TF controlled
#if COUNT_NEAR_AND
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                         	if(N_strong_BS[aux_TF][0][0]==1 || N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                genotype->N_near_AND_gated_motifs[4]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif
                                            }
                                            else 
                                            {
                                                genotype->N_motifs[26]+=cluster_size;  //OR-gated
#if COUNT_NEAR_AND
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                        {
                                                            if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                                abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                                )
                                                                genotype->N_near_AND_gated_motifs[5]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif
                                            }
                                        }                                     
                                    }
                                }
                                else // k is not regulated by j either, then we've found a diamond:)
                                {
                                    genotype->N_motifs[27]+=cluster_size;  //any logic
                                    
                                    /*temperarily assign j and k to master and aux. TF*/
                                    master_TF=genotype->which_protein[copies_reg_by_env[j]];
                                    aux_TF=genotype->which_protein[copies_reg_by_env[k]];
                                    /*check if the assignment is correct by comparing r_protein_deg*/
                                    /*calculate geometric mean of r_protein_deg across variants*/
                                    r_master_TF_deg=1.0;
                                    r_aux_TF_deg=1.0;
                                    for(l=0;l<genotype->protein_pool[master_TF][0][0];l++)
                                        r_master_TF_deg*=genotype->protein_decay_rate[genotype->protein_pool[master_TF][1][l]];
                                    r_master_TF_deg=pow(r_master_TF_deg,1.0/(float)genotype->protein_pool[master_TF][0][0]);
                                    for(l=0;l<genotype->protein_pool[aux_TF][0][0];l++)
                                        r_aux_TF_deg*=genotype->protein_decay_rate[genotype->protein_pool[aux_TF][1][l]];
                                    r_aux_TF_deg=pow(r_aux_TF_deg,1.0/(float)genotype->protein_pool[aux_TF][0][0]);
                                    /*master TF should have larger r_protein_deg*/
                                    master_TF=(r_master_TF_deg>r_aux_TF_deg)?genotype->which_protein[copies_reg_by_env[j]]:genotype->which_protein[copies_reg_by_env[k]];
                                    aux_TF=(r_master_TF_deg>r_aux_TF_deg)?genotype->which_protein[copies_reg_by_env[k]]:genotype->which_protein[copies_reg_by_env[j]];
                                    
                                    if(hindrance[master_TF][aux_TF])
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[28]+=cluster_size; //no regulation
                                            else  
                                                genotype->N_motifs[29]+=cluster_size; //emergent I3
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                                genotype->N_motifs[30]+=cluster_size; //emergent I1
                                            else  
                                                genotype->N_motifs[31]+=cluster_size; //impossible
                                        }    
                                    }
                                    else
                                    {
                                        if(hindrance[master_TF][master_TF])
                                        {
                                            if(hindrance[aux_TF][aux_TF])                                            
                                            {
                                                genotype->N_motifs[32]+=cluster_size; //AND-gated
#if DISABLE_AND_GATE                                        
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                             
#endif 
                                            }
                                            else  
                                            {
                                                genotype->N_motifs[33]+=cluster_size; //Aux. TF controlled
#if COUNT_NEAR_AND
                                                if(N_strong_BS[aux_TF][0][0]!=0 && N_strong_BS[master_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                         	if(N_strong_BS[master_TF][0][0]==1 || N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                genotype->N_near_AND_gated_motifs[6]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif 
                                            }
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[34]+=cluster_size; //master TF controlled
#if COUNT_NEAR_AND
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                         	if(N_strong_BS[aux_TF][0][0]==1 || N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                genotype->N_near_AND_gated_motifs[7]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif
                                            }
                                            else  
                                            {
                                                genotype->N_motifs[35]+=cluster_size; //OR-gated
#if COUNT_NEAR_AND
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                        {
                                                            if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                                abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                                )
                                                                genotype->N_near_AND_gated_motifs[8]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif
                                            }
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
                                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==protein_id && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_TF_TO_TF)
                                    {
                                        found_bs=1;
                                        break;
                                    }
                                }
                                if(!found_bs) // k is not regulated by j, we've found another FFL-in-diamond                           
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
                                            {
                                                genotype->N_motifs[23]+=cluster_size; //AND-gated
#if DISABLE_AND_GATE                                      
                                                genotype->gene_in_core_C1ffl[gene_id_copy]=1;
#if FORCE_MASTER_CONTROLLED
                                                genotype->TF_in_core_C1ffl[gene_id_copy][master_TF]=1;
#else
                                                genotype->TF_in_core_C1ffl[gene_id_copy][aux_TF]=1;
#endif                                       
#elif FORCE_DIAMOND                                       
                                                genotype->gene_in_core_C1ffl[copies_reg_by_env[j]]=1;
                                                genotype->TF_in_core_C1ffl[copies_reg_by_env[j]][genotype->which_TF_family[master_TF]]=1;                                        
#elif FORCE_SINGLE_FFL                                       
                                                genotype->gene_in_core_C1ffl[copies_reg_by_env[j]]=1; 
#endif
                                            }
                                            else 
                                            {
                                                genotype->N_motifs[24]+=cluster_size; //aux. TF controlled
#if COUNT_NEAR_AND
                                                if(N_strong_BS[aux_TF][0][0]!=0 && N_strong_BS[master_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                         	if(N_strong_BS[master_TF][0][0]==1 || N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                genotype->N_near_AND_gated_motifs[3]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif
                                            }
                                        }
                                        else
                                        {
                                            if(hindrance[aux_TF][aux_TF])
                                            {
                                                genotype->N_motifs[25]+=cluster_size; //master controlled
#if COUNT_NEAR_AND
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                            abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                            )
                                                        {
                                                         	if(N_strong_BS[aux_TF][0][0]==1 || N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)   
                                                                genotype->N_near_AND_gated_motifs[4]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif    
                                            }
                                            else  
                                            {
                                                genotype->N_motifs[26]+=cluster_size; //OR-gated
#if COUNT_NEAR_AND
                                                if(N_strong_BS[master_TF][0][0]!=0 && N_strong_BS[aux_TF][0][0]!=0)
                                                {
                                                    if(N_strong_BS[master_TF][0][0]==1 ||N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1]-N_strong_BS[master_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                    {
                                                        if(N_strong_BS[aux_TF][0][0]==1 ||N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[aux_TF][1][0]<TF_ELEMENT_LEN+2*HIND_LENGTH)
                                                        {
                                                            if(abs(N_strong_BS[aux_TF][1][0]-N_strong_BS[master_TF][1][N_strong_BS[master_TF][0][0]-1])>=TF_ELEMENT_LEN+2*HIND_LENGTH ||
                                                                abs(N_strong_BS[aux_TF][1][N_strong_BS[aux_TF][0][0]-1]-N_strong_BS[master_TF][1][0])>=TF_ELEMENT_LEN+2*HIND_LENGTH 
                                                                )
                                                                genotype->N_near_AND_gated_motifs[5]+=cluster_size;
                                                        }
                                                    }
                                                }
#endif    
                                            }
                                        }                                     
                                    }                                   
                                }
                                else
                                    genotype->N_motifs[38]+=cluster_size; // strong-connected
                            }
                        }
                    }
                }  
#endif
            }
            i++;
        }        
    }    
    else // if the genotype does not contain non-signal activator
    {
        for(i=0;i<39;i++)        
            genotype->N_motifs[i]=0;
    }
}

#if PERTURB
static void modify_topology(Genotype *genotype, Genotype *genotype_clone)
{
    int i, j, gene_id;
    
    for(i=0;i<39;i++)
        genotype_clone->N_motifs[i]=genotype->N_motifs[i];
    for(i=0;i<MAX_GENES;i++)
    {
        genotype_clone->gene_in_core_C1ffl[i]=genotype->gene_in_core_C1ffl[i];
        for(j=0;j<MAX_PROTEINS;j++)
            genotype_clone->TF_in_core_C1ffl[i][j]=genotype->TF_in_core_C1ffl[i][j];
    }
    
    for(gene_id=N_SIGNAL_TF;gene_id < genotype_clone->ngenes;gene_id++)
    {        
#if DISABLE_AND_GATE
        add_binding_site(genotype_clone, gene_id);
#else
        remove_binding_sites(genotype_clone, gene_id);
#endif
    }
}

static void add_binding_site(Genotype *genotype, int gene_id)
{    
    float temp;  
    int i,j;   
    /*make sure there is enough room*/
    if (genotype->binding_sites_num[gene_id] + 1 >= genotype->N_allocated_elements) 
    {  
        while(genotype->N_allocated_elements<=genotype->binding_sites_num[gene_id]+1)
            genotype->N_allocated_elements+=100;

        for(i=0;i<MAX_GENES;i++)
        {
            genotype->all_binding_sites[i] = realloc(genotype->all_binding_sites[i], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
            if (!genotype->all_binding_sites[i]) 
            {
#if !LOG_OFF                
                LOG("error in calc_all_binding_sites_copy\n");     
#endif
                exit(-1);                                       
            }     
        }                    
    }  
    
#if DIRECT_REG
    if(genotype->N_motifs[0]!=0 && genotype->gene_in_core_C1ffl[gene_id]==1) //modify Bs only if gene_id is regulated by AND gate  
#else
    if(genotype->gene_in_core_C1ffl[gene_id]==1)
#endif
    {
#if DIRECT_REG
        temp=1.0;
    #if FORCE_MASTER_CONTROLLED /* add a binding site of the signal to gene_id.*/        
        for(i=0;i<genotype->binding_sites_num[gene_id];i++)
        {
            if(genotype->all_binding_sites[gene_id][i].tf_id==N_SIGNAL_TF-1) 
                temp=(temp>genotype->all_binding_sites[gene_id][i].Kd)?genotype->all_binding_sites[gene_id][i].Kd:temp;  // find a strong binding site                      
        }
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].tf_id = N_SIGNAL_TF-1;
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd = temp;                   
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].BS_pos = 2*CISREG_LEN; //no way to block any existing BS
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].mis_match = 0;
        genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].N_hindered = 0;
        genotype->binding_sites_num[gene_id]++;
        genotype->N_act_BS[gene_id]++;
    #endif
#endif
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(genotype->TF_in_core_C1ffl[gene_id][i]==1)
            {
#if ADD_STRONG_TFBS
                temp=1.0;                
                for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                {
                    if(genotype->all_binding_sites[gene_id][j].tf_id==i) 
                        temp=(temp>genotype->all_binding_sites[gene_id][j].Kd)?genotype->all_binding_sites[gene_id][j].Kd:temp;  // find a strong binding site                      
                }
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd=temp;
#else
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd = KD2APP_KD*genotype->Kd[i]*pow(NS_Kd/genotype->Kd[i],(float)(TF_ELEMENT_LEN-NMIN)/(TF_ELEMENT_LEN-NMIN+1));  
#endif                
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
    /*the maximum number of activator binding sites that do not exclude each other is increased*/
    genotype->max_unhindered_sites[gene_id][1]++;    
}

/* ignore TF x when searching binding sites on gene y*/
/* this function is almost the same as calc_all_binding_sites_copy*/
static void remove_binding_sites(Genotype *genotype, int gene_id)
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
        
        for (k=start_TF;k<genotype->nproteins-1;k++) 
        { 
#if FORCE_DIAMOND
    #if !DIRECT_REG
            if(!(genotype->gene_in_core_C1ffl[gene_id]==1 && genotype->TF_in_core_C1ffl[gene_id][genotype->which_TF_family[k]]==1))
    #else
            if(!(  genotype->N_motifs[0]!=0 
                && genotype->gene_in_core_C1ffl[gene_id]==1
                && k==start_TF))
    #endif
            {
#elif FORCE_SINGLE_FFL        
            if(!(genotype->gene_in_core_C1ffl[gene_id]==1 && k==start_TF))       
            {    
#endif
            tf_seq=&(genotype->tf_seq[k][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[k][0]);            
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
#if FORCE_DIAMOND
            }
#elif FORCE_SINGLE_FFL
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
   
    int act_BS[MAXELEMENTS][2],rep_BS[MAXELEMENTS][2];
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
        else if(genotype->protein_identity[genotype->all_binding_sites[gene_id][i].tf_id]==REPRESSOR)
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

static void tidy_output_files(char *file_genotype_summary, char *file_mutations)
{
    int i,replay_N_steps;
    char buffer[600];    
    FILE *fp1,*fp2;
    
    fp1=fopen("saving_point.txt","r");
    fscanf(fp1,"%d",&replay_N_steps);
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
    
    fp1=fopen("N_motifs.txt","r");
    fp2=fopen("temp","w");
    for(i=0;i<replay_N_steps;i++)
    {
        fgets(buffer,600,fp1);
        fputs(buffer,fp2);
    }
    fclose(fp1);
    fclose(fp2);
    remove("N_motifs.txt");
    rename("temp","N_motifs.txt");
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
                                    genotype->translation_rate[i],
                                    genotype->protein_decay_rate[i],
                                    genotype->locus_length[i]);
        if(genotype->protein_identity[genotype->which_protein[i]]!=-1) //is a tf gene
            fprintf(fp,"%f\n",log10(genotype->Kd[genotype->which_protein[i]]));
        else
             fprintf(fp,"na\n");
    }
    fclose(fp);
}
