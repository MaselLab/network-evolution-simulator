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

static void set_signal(CellState *, Environment *, float, RngStream, int);

static void calc_avg_fitness(Genotype *, Selection *, int [MAX_GENES], float [MAX_PROTEINS], RngStream [N_THREADS], float *, float *); 

static void calc_fitness_stats(Genotype *, Selection *, float (*)[N_REPLICATES], float (*)[N_REPLICATES], int);

static void try_replacement(Genotype *, Genotype *, int *, float*);

static void clone_genotype(Genotype *, Genotype *);

static void summarize_binding_sites(Genotype *,int);

static int evolve_N_steps(Genotype *, Genotype *,  Mutation *, Selection *, Output_buffer [OUTPUT_INTERVAL], int *, int *, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS], int);

static void run_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int [MAX_GENES], float [MAX_PROTEINS], int, int, RngStream, RngStream [N_THREADS]);

static void continue_simulation(Genotype *, Genotype *, Mutation *, Selection *, Selection *, int, int [MAX_GENES], float [MAX_PROTEINS], RngStream, RngStream [N_THREADS]);

static void replay_mutations(Genotype *, Genotype *, Mutation *, int);

static void find_motifs(Genotype *);

static int find_TFBS_of_A_on_B(Genotype *, int, int);

static void find_activators_of_effector(Genotype *, int, int *, int [MAX_PROTEINS]);

static void build_hindrance_table(Genotype *, int, int, int [MAX_PROTEINS], int [MAX_PROTEINS][MAX_PROTEINS], int [MAX_PROTEINS][2][50]);

static void find_signal_regulated_genes(Genotype *, int, int [MAX_PROTEINS], int *, int *, int [MAX_GENES], int [MAX_GENES]);

static int find_isolated_C1FFL(Genotype *, int, int, int *, int *, int *, int *);

static int find_FFL_in_diamond(Genotype *, int, int, int *, int *, int *, int *);

static void determine_motif_logic(Genotype *, int [MAX_PROTEINS][MAX_PROTEINS], int, int, int, char, int [MAX_PROTEINS][2][50], int, int, int);

static void determine_near_AND_logic(Genotype *, int [MAX_PROTEINS][2][50], int, int, int, char, int);

static void tidy_output_files(char*, char*);

static void print_motifs(Genotype *);

static void store_mutant_info(Genotype *, Mutation *, Output_buffer *, int, int);

static void store_resident_info(Genotype *, Mutation *, Output_buffer *, int , int, int, float, int);

static void output_mutant_info(Output_buffer *, int);

static void output_resident_info(Output_buffer [OUTPUT_INTERVAL], int, int);

static void sample_motifs(Genotype *, Mutation *, int, RngStream);

static void sample_parameters(Genotype *, int, RngStream);

static void mark_genes_for_sampling(Genotype *, int, int, int, char);

static void mark_genes_to_be_perturbed(Genotype *, int, int, int, int, int);

static void modify_topology(Genotype *, Genotype *);

static void add_binding_site(Genotype *, int);

static void remove_binding_sites(Genotype *, int);

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
void evolve_neutrally(Genotype *resident, Genotype *mutant, Mutation *mut_record, Selection *burn_in, Selection *selection, RngStream RS_main)
{
    int i;    
    FILE *fp;    
    /*Create title for output files*/
    fp=fopen(evo_summary,"a+");
    fprintf(fp,"step N_tot_mut_tried N_mut_tried_this_step N_hit_bound accepted_mut selection_coeff avg_fitness fitness1 fitness2 se_avg_fitness se_fitness1 se_fitness2 N_genes N_proteins N_activator N_repressor\n");
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
    /*sampling parameters of network motifs*/
    if(SAMPLE_PARAMETERS)
        sample_motifs(resident, mut_record, selection->MAX_STEPS, RS_parallel[0]);
    
    /*replay mutations, output N_motifs.txt and networks.txt*/   
    if(REPRODUCE_GENOTYPES || SAMPLE_GENE_EXPRESSION)
    {        
        replay_mutations(resident, mutant, mut_record, selection->MAX_STEPS);    
        /*output the evolved genotype*/
        calc_all_binding_sites(resident); 
        print_mutatable_parameters(resident,1);
        summarize_binding_sites(resident,selection->MAX_STEPS);  
    }
    
    if(SAMPLE_GENE_EXPRESSION)
    {
        /*create threads*/
        omp_set_num_threads(N_THREADS); 
        /*collection interval is 1 minute by default*/    
        selection->env1.t_development=91.1; //to make 90 data points
        selection->env2.t_development=91.1; 
        calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, NULL, NULL);  
    }
}

/*randomly sampling motifs*/
static void sample_motifs(Genotype *genotype_ori,                                             
                            Mutation *mut_record,
                            int max_step,                           
                            RngStream RS)
{
    int i,N_samples,which_step;    
    Genotype genotype_copy1, genotype_copy2;
    initialize_cache(&genotype_copy1); 
    initialize_cache(&genotype_copy2);
    FILE *fp;
    int buffer_int;
    float buffer_float;
    char buffer_char;
    char buffer_string[3];
    
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
    /*Don't modify genotype_ori because it will be used elsewhere*/
    clone_genotype(genotype_ori,&genotype_copy1);
    
    /*Only the last x steps will be repeated sampled. Let's keep a copy of 
     *the genotype at the last (x-1)th step, so re-sampling can alway start here.
     */
    for(i=0;i<START_STEP_OF_SAMPLING-1;i++)
    {    
        fscanf(fp,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                        &(mut_record->which_gene),
                                        &(mut_record->which_nucleotide), 
                                        mut_record->nuc_diff,               
                                        &(mut_record->kinetic_type),
                                        &(mut_record->kinetic_diff));
        reproduce_mutate(&genotype_copy1,mut_record);  
    }
    /*close mutation_file*/    
    fclose(fp);  
    /*reset N_samples*/
    N_samples=0;    
    /*sampling repeatedly*/
    while(N_samples<SAMPLE_SIZE)
    {
        /*Keep the genotype at the last (x-1)th step intact*/
        clone_genotype(&genotype_copy1,&genotype_copy2);  
        which_step=RngStream_RandInt(RS,1,max_step-START_STEP_OF_SAMPLING+1);
        /*open mutation record and skip the entries before start_step*/
        fp=fopen(mutation_file,"r");
        for(i=0;i<START_STEP_OF_SAMPLING-1;i++)
            fscanf(fp,"%c %d %d %s %d %a\n",&buffer_char,&buffer_int,&buffer_int,&(buffer_string[0]),&buffer_int, &buffer_float);
        /*reproduce the genotype at which_step*/
        for(i=0;i<which_step;i++)
        { 
            fscanf(fp,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                        &(mut_record->which_gene),
                                                        &(mut_record->which_nucleotide), 
                                                        mut_record->nuc_diff,               
                                                        &(mut_record->kinetic_type),
                                                        &(mut_record->kinetic_diff));
            reproduce_mutate(&genotype_copy2,mut_record); 
        } 
        /*score motifs*/
		calc_all_binding_sites(&genotype_copy2); 
        find_motifs(&genotype_copy2); 
        /*does the target motif exist*/
        switch(TARGET_MOTIF) 
        {
            case 0:
                    sample_parameters(&genotype_copy2,i,RS);
                    N_samples++;
                    break;
            case 1:
                    if(genotype_copy2.N_motifs[5]!=0)
                    {
                        sample_parameters(&genotype_copy2,i,RS);
                        N_samples++;
                    }
                    break;
            case 2:
                    if(genotype_copy2.N_motifs[14]!=0)
                    {
                        sample_parameters(&genotype_copy2,i,RS);
                        N_samples++;
                    }
                    break;
            case 3:
                    if(genotype_copy2.N_motifs[32]!=0)
                    {
                        sample_parameters(&genotype_copy2,i,RS);
                        N_samples++;
                    }                
        }      
        /*close mutation record*/
        fclose(fp);   
    }
    printf("Sampling parameters successfully!\n");
}

void mark_genes_for_sampling(Genotype *genotype, int effector_gene_id, int fast_TF, int slow_TF, char which_motif)
{
    int mark_genes;    
    /*set the flag to 0*/
    mark_genes=0;    
    /*raise the flag if which_motif is the target motif*/
    switch(which_motif)
    {
        case 'D':           
                if(TARGET_MOTIF==1)
                    mark_genes=1;
                    break;
        case 'C':            
                if(TARGET_MOTIF==2)
                    mark_genes=1;
                    break;        
        case 'I':            
                if(TARGET_MOTIF==3)
                    mark_genes=1;         
    }
    /*mark genes*/
    if(mark_genes)
    {
        genotype->cis_target_to_be_perturbed[effector_gene_id]=YES;
        genotype->trans_target_to_be_perturbed[effector_gene_id][fast_TF]=YES;
        genotype->slow_TF[effector_gene_id][slow_TF]=YES;            
    }
}

/*output parameters of network motifs*/
void sample_parameters(Genotype *genotype, int step, RngStream RS)
{
    int gene_id,protein_id,cluster_id,effector_gene_id;
    FILE *fp;
    
    /*make output file*/        
    fp=fopen("parameters.txt","a+");
    
    /*output the Kd of the signal*/
    fprintf(fp,"%d 0 0.0 0.0 0.0 0.0 0.0 %f\n",step,log10(genotype->Kd[0]));
    
    /*sampling an effector gene*/
    while(1)
    {
        gene_id=RngStream_RandInt(RS,1,genotype->ngenes-1);        
        /*if sample a neutral network*/
        if(TARGET_MOTIF==0)
        {
            /*then just pick an effector gene*/
            if(genotype->which_protein[gene_id]==genotype->nproteins-1)
                break;
        }
        else
        {
            cluster_id=genotype->which_cluster[gene_id];
            /*otherwise, need to pick an effector gene that's in the target motif*/
            if(genotype->cis_target_to_be_perturbed[genotype->cisreg_cluster[cluster_id][0]]==1)       
                break;
        }
    }
    fprintf(fp,"%d -1 %f %f %f %f %f 0.0\n",step,log10(genotype->active_to_intermediate_rate[gene_id]),
                                                log10(genotype->mRNA_decay_rate[gene_id]),
                                                log10(genotype->translation_rate[gene_id]),
                                                log10(genotype->protein_decay_rate[gene_id]),
                                                (float)genotype->locus_length[gene_id]); 
    
    /*sample a fast-TF gene from the motif that contain the chosen effector gene*/
    effector_gene_id=genotype->cisreg_cluster[cluster_id][0]; 
    if(TARGET_MOTIF!=1) //under direct regulation, the signal is the fast TF
    {       
        while(1)
        {
            gene_id=RngStream_RandInt(RS,1,genotype->ngenes-1);
            protein_id=genotype->which_protein[gene_id];
            /*if sample a neutral network*/
            if(TARGET_MOTIF==0)
            {        
                /*then just pick a TF protein*/
                if(protein_id!=genotype->nproteins-1)
                    break;
            }
            else
            {   /*otherwise need to pick a TF that's in the target motif*/
                if(genotype->trans_target_to_be_perturbed[effector_gene_id][protein_id]==1)
                    break;
            }
        }
    }
    else
    {
        gene_id=N_SIGNAL_TF-1;    
        protein_id=N_SIGNAL_TF-1;
    }    
    fprintf(fp,"%d 1 %f %f %f %f %f %f\n",step,log10(genotype->active_to_intermediate_rate[gene_id]),
                                            log10(genotype->mRNA_decay_rate[gene_id]),
                                            log10(genotype->translation_rate[gene_id]),
                                            log10(genotype->protein_decay_rate[gene_id]),
                                            (float)genotype->locus_length[gene_id],
                                            log10(genotype->Kd[protein_id])); 
    /*sample a slow-TF gene*/
    if(TARGET_MOTIF!=0) //no need to sample another TF in a neutral network, 
    {
        while(1)
        {
            gene_id=RngStream_RandInt(RS,1,genotype->ngenes-1);
            protein_id=genotype->which_protein[gene_id];          
            if(genotype->slow_TF[effector_gene_id][protein_id]==1)
                break;
        }
    }
    fprintf(fp,"%d 1 %f %f %f %f %f %f\n",step,log10(genotype->active_to_intermediate_rate[gene_id]),
                                            log10(genotype->mRNA_decay_rate[gene_id]),
                                            log10(genotype->translation_rate[gene_id]),
                                            log10(genotype->protein_decay_rate[gene_id]),
                                            (float)genotype->locus_length[gene_id],
                                            log10(genotype->Kd[protein_id]));
    fclose(fp);
}
#endif

#if PERTURB
void perturbation_analysis(Genotype *resident,
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
    file_mutation=fopen(mutation_file,"r");    
    if(file_mutation!=NULL)        
        printf("LOAD MUTATION RECORD SUCCESSFUL!\n");
    else
    {
        printf("Loading mutation record failed! Quit program!");
#if MAKE_LOG
        LOG("Loading mutation record failed!");
#endif
        exit(-2);
    } 
    
    /*skip first 2 rows of fitness_record*/
    fitness_record=fopen(evo_summary,"r");
    fgets(buffer,600,fitness_record);
    fgets(buffer,600,fitness_record);   
    
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
                &mean_fitness1,
                &mean_fitness2,
                &se_overall_fitness,
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
            {
                for(j=0;j<HI_RESOLUTION_RECALC;j++) 
                    calc_avg_fitness(resident, selection, init_mRNA, init_protein, RS_parallel, fitness1[j], fitness2[j]);                
                calc_fitness_stats(resident, selection, &(fitness1[0]), &(fitness2[0]), HI_RESOLUTION_RECALC);  
                f_aft_perturbation=fopen("f_aft_perturbation.txt","a+");
                fprintf(f_aft_perturbation,"%d %.10f %.10f %.10f %.10f %.10f %.10f\n",i,
                        resident->avg_fitness,                        
                        resident->fitness1,
                        resident->fitness2,
                        resident->SE_avg_fitness,
                        resident->SE_fitness1,
                        resident->SE_fitness2);
                fclose(f_aft_perturbation);
                /*file to store the original fitness of the to-be-modified network*/
                f_bf_perturbation=fopen("f_bf_perturbation.txt","a+");
                fprintf(f_bf_perturbation,"%d %.10f %.10f %.10f %.10f %.10f %.10f\n",
                        step,
                        mean_overall_fitness,
                        mean_fitness1,
                        mean_fitness2,
                        se_overall_fitness,
                        se_fitness1,
                        se_fitness2);  
                fclose(f_bf_perturbation);
            }    
            for(j=0;j<MAX_GENES;j++)
            {
                resident->cis_target_to_be_perturbed[j]=0;
                for(k=0;k<MAX_PROTEINS;k++)
                    resident->trans_target_to_be_perturbed[j][k]=0;
            }        
        }
    }      
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
    initialize_sequence((char *)genotype->tf_seq, CONSENSUS_SEQ_LEN*MAX_TF_GENES, genotype->ntfgenes, RS);    //initialize binding sequence of TFs    
    /* We now generate the complementary sequence of BS that are on the non-template strand.
     * The complementary sequence is used to search for BS that on the non-template strand.  
     * We also assume that all the TFs can work on both strands, but can induce expression in one direction.*/  
    for(i=0;i< genotype->ntfgenes;i++)
    {        
        for(k=0;k<CONSENSUS_SEQ_LEN;k++)
        {
            switch (genotype->tf_seq[i][CONSENSUS_SEQ_LEN-k-1])
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
  
    for(i=3; i < CISREG_LEN-CONSENSUS_SEQ_LEN-3; i++) /* scan promoter */
    {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(genotype->all_binding_sites[gene_id][j].BS_pos> i-CONSENSUS_SEQ_LEN-2*HIND_LENGTH)
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
            for (j=i; j < i+CONSENSUS_SEQ_LEN; j++) /*calculate the number of nucleotides that match in each [i,i+CONSENSUS_SEQ_LEN] window. The window slides by 1 each time when scanning the promoter*/
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
#if MAKE_LOG
                            LOG("error in calc_all_binding_sites_copy\n");  
#endif
                            exit(-1);                                       
                        }     
                    }                    
                }                
                genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                      
                genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(CONSENSUS_SEQ_LEN-match)/(CONSENSUS_SEQ_LEN-NMIN+1));
                genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ; 
                genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = CONSENSUS_SEQ_LEN-match;             
                genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                N_hindered_BS++;              
                N_binding_sites++;
                if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand.*/
            {
                match_rc=0;
                for (j=i; j < i+CONSENSUS_SEQ_LEN; j++)                
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
#if MAKE_LOG
                                LOG("error in calc_all_binding_sites_copy\n");
#endif
                                exit(-1);                                       
                            }     
                        }                    
                    }
                    /************************************************************************************************************/
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                                     
                    genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(CONSENSUS_SEQ_LEN-match_rc)/(CONSENSUS_SEQ_LEN-NMIN+1));
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
                    genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = CONSENSUS_SEQ_LEN-match_rc;
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
        while(j!=0 && genotype->all_binding_sites[gene_id][act_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][act_BS[j][0]].BS_pos<CONSENSUS_SEQ_LEN+2*HIND_LENGTH)j--;
        act_BS[i][1]=act_BS[j][1]+1;
    }
    /*calculate the maximum number of repressor binding sites that do not hinder each other*/
    rep_BS[0][0]=-1;
    rep_BS[0][1]=0;
    for(i=1;i<N_rep_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][rep_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][rep_BS[j][0]].BS_pos<CONSENSUS_SEQ_LEN+2*HIND_LENGTH)j--;
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
            genotype->recalc_TFBS[gene_id]=NO;
        }
    }
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
        state->protein_number[N_SIGNAL_TF-1]=0.0; 
        if(t_burn_in!=0.0)
            t=t+t_burn_in; //the completion of burn_in is a fixed event, which is added in initialize_cell       
        /*after burn-in, signal should be turned "o"n*/
        flag='o'; 
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
        for(j=0;j<CONSENSUS_SEQ_LEN;j++)
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
 *average fitness over tdevelopment.
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
        timecourse1[i].total_time_points=(int)(Selection->env1.t_development+Selection->env1.max_duration_of_burn_in_growth_rate);
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
        timecourse2[i].total_time_points=(int)(Selection->env2.t_development+Selection->env2.max_duration_of_burn_in_growth_rate);
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
              
        /*initialize the clone*/
        initialize_cache(&genotype_clone);        
        
        /*clone genotype and initial mRNA and protein numbers*/
        #pragma omp critical
        { 
            genotype_clone.ngenes=genotype->ngenes;
            genotype_clone.ntfgenes=genotype->ntfgenes;
            genotype_clone.nproteins=genotype->nproteins;
            clone_genotype(genotype, &genotype_clone);              
            Env1.t_development=Selection->env1.t_development;
            Env1.signal_on_strength=Selection->env1.signal_on_strength;
            Env1.signal_off_strength=Selection->env1.signal_off_strength;
            Env1.signal_on_aft_burn_in=Selection->env1.signal_on_aft_burn_in;
            Env1.t_signal_on=Selection->env1.t_signal_on;
            Env1.t_signal_off=Selection->env1.t_signal_off;
            Env1.initial_effect_of_effector=Selection->env1.initial_effect_of_effector;
            Env1.effect_of_effector_aft_burn_in=Selection->env1.effect_of_effector_aft_burn_in;
            Env1.fixed_effector_effect=Selection->env1.fixed_effector_effect;
            Env1.max_duration_of_burn_in_growth_rate=Selection->env1.max_duration_of_burn_in_growth_rate;
            Env1.avg_duration_of_burn_in_growth_rate=Selection->env1.avg_duration_of_burn_in_growth_rate;
            Env2.t_development=Selection->env2.t_development;
            Env2.signal_on_strength=Selection->env2.signal_on_strength;
            Env2.signal_off_strength=Selection->env2.signal_off_strength;
            Env2.signal_on_aft_burn_in=Selection->env2.signal_on_aft_burn_in;
            Env2.t_signal_on=Selection->env2.t_signal_on;
            Env2.t_signal_off=Selection->env2.t_signal_off;
            Env2.initial_effect_of_effector=Selection->env2.initial_effect_of_effector;
            Env2.effect_of_effector_aft_burn_in=Selection->env2.effect_of_effector_aft_burn_in;
            Env2.fixed_effector_effect=Selection->env2.fixed_effector_effect;
            Env2.max_duration_of_burn_in_growth_rate=Selection->env2.max_duration_of_burn_in_growth_rate;
            Env2.avg_duration_of_burn_in_growth_rate=Selection->env2.avg_duration_of_burn_in_growth_rate;
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
        /* now calc fitness under the two environments*/
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
                      
            /*calculate average instantaneous fitness of tdevelopment*/
            f1[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/Env1.t_development; 
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
        
            f2[i]=(state_clone.cumulative_fitness-state_clone.cumulative_fitness_after_burn_in)/Env2.t_development;
#if PHENOTYPE
            timecourse2[thread_ID*N_replicates_per_thread+i].timepoint++;
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
                fprintf(fp,"%f ",timecourse1[j].protein_concentration[k+i*timecourse1[j].total_time_points]);
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
                fprintf(fp,"%f ",timecourse1[j].gene_specific_concentration[k+i*timecourse1[j].total_time_points]);
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
                fprintf(fp,"%f ",timecourse2[j].gene_specific_concentration[k+i*timecourse2[j].total_time_points]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    } 
   
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
    genotype->N_allocated_elements=MAX_TFBS_NUMBER;
    for(j=0;j<MAX_GENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(MAX_TFBS_NUMBER*sizeof(AllTFBindingSites));
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


static void replay_mutations(Genotype *resident, Genotype *mutant, Mutation *mut_record, int replay_N_steps)
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
        clone_genotype(resident,mutant);
        fscanf(fp,"%c %d %d %s %d %a\n",&(mut_record->mut_type),
                                                    &(mut_record->which_gene),
                                                    &(mut_record->which_nucleotide), 
                                                    mut_record->nuc_diff,               
                                                    &(mut_record->kinetic_type),
                                                    &(mut_record->kinetic_diff));
        reproduce_mutate(mutant,mut_record);        
        clone_genotype(mutant,resident);        
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
                            int init_step,              //init_step is either 1 or loaded from a saving point
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
        N_EFFECTOR_GENES=burn_in->temporary_N_effector_genes;
        N_TF_GENES=burn_in->temporary_N_tf_genes; 
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
    N_EFFECTOR_GENES=selection->temporary_N_effector_genes;
    N_TF_GENES=selection->temporary_N_tf_genes; 
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
    replay_mutations(resident, mutant, mut_record, replay_N_steps); 

    /* load random number seeds*/
    fp=fopen("RngSeeds.txt","r");
    if(fp!=NULL)
    {
        for(i=0;i<replay_N_steps/SAVING_INTERVAL;i++)
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
#if MAKE_LOG
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
            if(i%SAVING_INTERVAL==0)
            {
                fp=fopen("saving_point.txt","w");
                fprintf(fp,"%d %d\n",i,*N_tot_trials);
                fflush(fp);
                fclose(fp);
            }    
        }
    } 
    *init_step=i;
    free(mutant_info);
    return 0;
}

static void print_motifs(Genotype *genotype)
{
    FILE *fp; 
    int i;
    fp=fopen("N_motifs.txt","a+");
    for(i=0;i<36;i++)
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
            if(genotype->which_protein[i]==genotype->nproteins-1)
            {
                if(genotype->all_binding_sites[i][j].tf_id==N_SIGNAL_TF-1)
                {
                    if(genotype->all_binding_sites[i][j].mis_match<=CUT_OFF_MISMATCH_SIG2EFFECTOR)
                        table[i][genotype->all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
                }
                else
                {
                    if(genotype->all_binding_sites[i][j].mis_match<=CUT_OFF_MISMATCH_TF2EFFECTOR)
                        table[i][genotype->all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
                }                
            }
            else
            {
                if(genotype->all_binding_sites[i][j].tf_id==N_SIGNAL_TF-1)
                {
                    if(genotype->all_binding_sites[i][j].mis_match<=CUT_OFF_MISMATCH_SIGNAL2TF)
                        table[i][genotype->all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
                }
                else
                {
                    if(genotype->all_binding_sites[i][j].mis_match<=CUT_OFF_MISMATCH_TF2TF)
                        table[i][genotype->all_binding_sites[i][j].tf_id]++; /*the numbers of the binding sites of each tf on promoter i*/
                }          
            }
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
    int i,j,k;
    int cluster_size; 
    int found_what;
    int effector_gene_id, fast_TF_gene_id, slow_TF_gene_id;  
    int fast_TF,slow_TF;
    int N_activators, activators[MAX_PROTEINS];
    int copies_reg_by_env[MAX_GENES],copies_not_reg_by_env[MAX_GENES],N_copies_reg_by_env,N_copies_not_reg_by_env;
    int hindrance[MAX_PROTEINS][MAX_PROTEINS];    
    int strong_BS_pos[MAX_PROTEINS][2][50];    // 50 slots should be enough to store TFBSs positions on a gene
    
    /*reset records*/
    for(i=0;i<36;i++)
        genotype->N_motifs[i]=0; 
#if PERTURB
    for(i=0;i<MAX_GENES;i++)
    {
        genotype->cis_target_to_be_perturbed[i]=0;
        for(j=0;j<MAX_PROTEINS;j++)
            genotype->trans_target_to_be_perturbed[i][j]=0;
    }
#elif SAMPLE_PARAMETERS
    for(i=0;i<MAX_GENES;i++)
    {
        genotype->cis_target_to_be_perturbed[i]=0;
        for(j=0;j<MAX_PROTEINS;j++)
        {
            genotype->trans_target_to_be_perturbed[i][j]=0;
            genotype->slow_TF[i][j]=0;
        }
    }
#endif
#if COUNT_NEAR_AND
        for(i=0;i<12;i++)
            genotype->N_near_AND_gated_motifs[i]=0;       
#endif  
    
    /*look for motifs only if the genotype contains non-signal activators*/
    if(genotype->N_act>N_SIGNAL_TF)
    { 
        /*loop through each cisreg_cluster. Genes in a cluster have the same binding sites*/
        i=0;   
        while(genotype->cisreg_cluster[i][0]!=NA) //not an empty cluster
        {      
            if(genotype->which_protein[genotype->cisreg_cluster[i][0]]==genotype->nproteins-1) // is a effector gene
            {       
                effector_gene_id=genotype->cisreg_cluster[i][0];
                /*find how many genes in this cluster*/
                cluster_size=0;
                while(genotype->cisreg_cluster[i][cluster_size]!=NA)
                    cluster_size++;
#if COUNT_NEAR_AND
                /*reset table for recording strength of interaction*/
                for(j=0;j<MAX_PROTEINS;j++)
                {                   
                    strong_BS_pos[j][0][0]=0;
                    for(k=0;k<50;k++)
                        strong_BS_pos[j][1][k]=0;
                }
#endif  
                /*list activators that regulate the effector gene*/
                find_activators_of_effector(genotype, effector_gene_id, &N_activators, activators);                
                /*build a table to show hindrance between binding sites on effector gene*/
                build_hindrance_table(genotype, effector_gene_id, N_activators, activators, hindrance, strong_BS_pos);                
                /*find genes regulated by the environmental signal*/
                find_signal_regulated_genes(genotype, N_activators, activators, &N_copies_reg_by_env, &N_copies_not_reg_by_env, copies_reg_by_env, copies_not_reg_by_env);  
#if DIRECT_REG  
                /*count c1-ffls formed by signal and one signal-regulated copy */
                if(activators[0]==N_SIGNAL_TF-1)
                {                    
                    for(j=0;j<N_copies_reg_by_env;j++)
                    {
                        slow_TF=genotype->which_protein[copies_reg_by_env[j]];
                        determine_motif_logic(genotype, 
                                                hindrance, 
                                                N_SIGNAL_TF-1, 
                                                slow_TF, 
                                                cluster_size, 
                                                'D', // "D"iret regulation
                                                strong_BS_pos,
                                                effector_gene_id,
                                                N_SIGNAL_TF-1,
                                                copies_reg_by_env[j]);
                    }
                }
#else
                /*count isolated c1-ffl*/
                for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=0;k<N_copies_not_reg_by_env;k++)
                    {
                        if(find_isolated_C1FFL(genotype,copies_reg_by_env[j],copies_not_reg_by_env[k], &fast_TF, &slow_TF, &fast_TF_gene_id, &slow_TF_gene_id))
                        {                            
                            determine_motif_logic(genotype, 
                                                    hindrance,
                                                    fast_TF, 
                                                    slow_TF, 
                                                    cluster_size, 
                                                    'C', //isolated C1-FFL
                                                    strong_BS_pos,
                                                    effector_gene_id,
                                                    fast_TF_gene_id,
                                                    slow_TF_gene_id);   
                        }
                    }
                }                      
               
                /*count motifs formed by two env-regulated copies*/               
                for(j=0;j<N_copies_reg_by_env;j++)
                {
                    for(k=j+1;k<N_copies_reg_by_env;k++)
                    {
                        found_what=find_FFL_in_diamond(genotype, copies_reg_by_env[j], copies_reg_by_env[k], &fast_TF, &slow_TF, &fast_TF_gene_id, &slow_TF_gene_id);
                        if(found_what==1) //found a FFL-in-diamond                      
                            determine_motif_logic(genotype, 
                                                    hindrance,
                                                    fast_TF, 
                                                    slow_TF, 
                                                    cluster_size, 
                                                    'F', //FFL-in-diamond
                                                    strong_BS_pos,
                                                    effector_gene_id,
                                                    fast_TF_gene_id,
                                                    slow_TF_gene_id);  
                        else if(found_what==2) //found a diamond:)
                            determine_motif_logic(genotype, 
                                                    hindrance,
                                                    fast_TF, 
                                                    slow_TF, 
                                                    cluster_size, 
                                                    'I', //dImond
                                                    strong_BS_pos,
                                                    effector_gene_id,
                                                    fast_TF_gene_id,
                                                    slow_TF_gene_id);   
                    }
                }  
#endif
            }
            i++;
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

static void find_activators_of_effector(Genotype *genotype, int effector_gene, int *N_activators, int activators[MAX_PROTEINS])
{
    int i,j;
    int protein_id;
    
    /*reset table for recording activators that regulate effector gene*/
    for(i=0;i<MAX_PROTEINS;i++)
        activators[i]=0;
    
    /*scan binding sites for tfs that regulate gene_id*/
    for(i=0;i<genotype->binding_sites_num[effector_gene];i++)
    {
        protein_id=genotype->all_binding_sites[effector_gene][i].tf_id;
        if(genotype->protein_identity[protein_id]==ACTIVATOR)
        {
            if(protein_id==N_SIGNAL_TF-1)
            {
                if(genotype->all_binding_sites[effector_gene][i].mis_match<=CUT_OFF_MISMATCH_SIG2EFFECTOR)
                    activators[protein_id]=1;
            }
            else
            {
                if(genotype->all_binding_sites[effector_gene][i].mis_match<=CUT_OFF_MISMATCH_TF2EFFECTOR)
                    activators[protein_id]=1;
            }                        
        }    
    }
    
    /* move non-zeros entries in activators to the front. */
    i=0; //marks the entry to which a record is copied
    j=0;
    *N_activators=0;
    while(i<genotype->nproteins)
    {
        if(activators[i]!=0)               
        {
            activators[j]=i;
            if(j!=i)
                activators[i]=0;
            j++; 
            (*N_activators)++; //total number of activators that regulate effector
        }    
        i++;
    } 
}


/*if every binding site of i can hinder all the binding site of j, denote hindrance_table[i][j]=1. Otherwise 0*/
/*hindrance_table[i][i]=1 means at most one TF i can bind to the effector gene */
/*A strict AND gate between i and j should have H[i][i]=H[j][j]=1, and H[i][j]=0*/
static void build_hindrance_table(Genotype *genotype, 
                                    int effector_gene, 
                                    int N_activators, 
                                    int activators[MAX_PROTEINS], 
                                    int hindrance_table[MAX_PROTEINS][MAX_PROTEINS],
                                    int strong_BS_pos[MAX_PROTEINS][2][50])
{
    int i,j,k;
    int site_id;
    int pos_binding_sites_of_i[MAX_TFBS_NUMBER],pos_binding_sites_of_j[MAX_TFBS_NUMBER],N_binding_sites_of_i,N_binding_sites_of_j; 

    /*reset hindrance table*/
    for(i=0;i<MAX_PROTEINS;i++)
    {
        for(j=0;j<MAX_PROTEINS;j++)
            hindrance_table[i][j]=1;
    }
    
    for(i=0;i<N_activators;i++)
    {
        /*reset position table*/
        for(j=0;j<MAX_TFBS_NUMBER;j++)
            pos_binding_sites_of_i[j]=-CISREG_LEN; 
        /*list the positions of all the binding sites of activator i*/
        N_binding_sites_of_i=0;       
        for(site_id=0;site_id<genotype->binding_sites_num[effector_gene];site_id++)
        {
            if(genotype->all_binding_sites[effector_gene][site_id].tf_id==activators[i])
            {
                if(genotype->all_binding_sites[effector_gene][site_id].tf_id==N_SIGNAL_TF-1)
                {
                    if(genotype->all_binding_sites[effector_gene][site_id].mis_match<=CUT_OFF_MISMATCH_SIG2EFFECTOR)
                    {
                        pos_binding_sites_of_i[N_binding_sites_of_i]=genotype->all_binding_sites[effector_gene][site_id].BS_pos;
                        N_binding_sites_of_i++;
#if COUNT_NEAR_AND
                        if(genotype->all_binding_sites[effector_gene][site_id].mis_match<(CONSENSUS_SEQ_LEN-NMIN))
                        {                                
                            strong_BS_pos[activators[i]][1][strong_BS_pos[activators[i]][0][0]]=genotype->all_binding_sites[effector_gene][site_id].BS_pos;
                            strong_BS_pos[activators[i]][0][0]++;
                        }
#endif
                    }
                }
                else
                {
                    if(genotype->all_binding_sites[effector_gene][site_id].mis_match<=CUT_OFF_MISMATCH_TF2EFFECTOR)
                    {
                        pos_binding_sites_of_i[N_binding_sites_of_i]=genotype->all_binding_sites[effector_gene][site_id].BS_pos;
                        N_binding_sites_of_i++;
#if COUNT_NEAR_AND
                        if(genotype->all_binding_sites[effector_gene][site_id].mis_match<(CONSENSUS_SEQ_LEN-NMIN))
                        {                                
                            strong_BS_pos[activators[i]][1][strong_BS_pos[activators[i]][0][0]]=genotype->all_binding_sites[effector_gene][site_id].BS_pos;
                            strong_BS_pos[activators[i]][0][0]++;
                        }
#endif
                    }
                }                                
            }
        }
        
        /*if the first bs hinders the last bs, we can be sure that effecto can be bound by at most one TF i */
        if(pos_binding_sites_of_i[N_binding_sites_of_i-1]-pos_binding_sites_of_i[0]>=CONSENSUS_SEQ_LEN+2*HIND_LENGTH)    
            hindrance_table[activators[i]][activators[i]]=0;
       

        /*list the positions of all the binding sites of activator i+1...*/                  
        N_binding_sites_of_j=0;
        for(j=i+1;j<N_activators;j++)
        {
            for(k=0;k<MAX_TFBS_NUMBER;k++)
                pos_binding_sites_of_j[k]=-CISREG_LEN;
            N_binding_sites_of_j=0;                        
            for(site_id=0;site_id<genotype->binding_sites_num[effector_gene];site_id++)
            {
                if(genotype->all_binding_sites[effector_gene][site_id].tf_id==activators[j])                                
                { 
                    if(genotype->all_binding_sites[effector_gene][site_id].mis_match<=CUT_OFF_MISMATCH_TF2EFFECTOR)   
                    {
                        pos_binding_sites_of_j[N_binding_sites_of_j]=genotype->all_binding_sites[effector_gene][site_id].BS_pos;
                        N_binding_sites_of_j++;
                    }                    
                } 
            }  
            /*determine the relation between i and i+1, ...*/
            if(abs(pos_binding_sites_of_i[N_binding_sites_of_i-1]-pos_binding_sites_of_j[0])>= CONSENSUS_SEQ_LEN+2*HIND_LENGTH ||
                abs(pos_binding_sites_of_j[N_binding_sites_of_j-1]-pos_binding_sites_of_i[0])>= CONSENSUS_SEQ_LEN+2*HIND_LENGTH)
            {
                hindrance_table[activators[i]][activators[j]]=0; 
                hindrance_table[activators[j]][activators[i]]=0;
            }
        }
    }
}

/*make lists of gene copies that are regulated by environmental signal and of those that are not*/
static void find_signal_regulated_genes(Genotype *genotype, 
                                        int N_activators,
                                        int activators[MAX_PROTEINS], 
                                        int *N_regulated_genes, 
                                        int *N_unregulated_genes, 
                                        int regulated_gene_ids[MAX_GENES], 
                                        int unregualted_gene_ids[MAX_GENES])
{
    int i,j;
    int N_copies;
    int site_id, gene_id;
    int found_bs;
    
    /*reset counters and tables*/
    *N_regulated_genes=0;
    *N_unregulated_genes=0;    
    for(i=0;i<MAX_GENES;i++)
    {
        regulated_gene_ids[i]=-1;
        unregualted_gene_ids[i]=-1;
    }
    
    /*loop through each gene copy of each activator*/
    for(i=0;i<N_activators;i++)
    {
        if(activators[i]>=N_SIGNAL_TF)
        {  
            N_copies=genotype->protein_pool[activators[i]][0][0]; // the number of genes encoding i
            for(j=0;j<N_copies;j++)
            {
                gene_id=genotype->protein_pool[activators[i]][1][j];
                /*reset flag*/
                found_bs=0;
                for(site_id=0;site_id<genotype->binding_sites_num[gene_id];site_id++)
                {
                    if(genotype->all_binding_sites[gene_id][site_id].tf_id==N_SIGNAL_TF-1 && genotype->all_binding_sites[gene_id][site_id].mis_match<=CUT_OFF_MISMATCH_SIGNAL2TF)
                    {
                        found_bs=1;
                        break;
                    }
                } 
                if(found_bs)
                {
                    regulated_gene_ids[*N_regulated_genes]=gene_id;
                    (*N_regulated_genes)++;                                
                }
                else
                {
                    unregualted_gene_ids[*N_unregulated_genes]=gene_id;
                    (*N_unregulated_genes)++;                                 
                }
            }                        
        }
    } 
}

static int find_isolated_C1FFL(Genotype *genotype, int signal_regulated_gene, int unregulated_gene, int *fast_TF, int *slow_TF, int *fast_TF_gene_id, int *slow_TF_gene_id)
{   
    int return_val=0;
    
    /*signal_regulated_gene and unregulated_gene must encode TFs of different families*/
    if(genotype->which_TF_family[genotype->which_protein[signal_regulated_gene]]!=genotype->which_TF_family[genotype->which_protein[unregulated_gene]])
    { 
        if(!find_TFBS_of_A_on_B(genotype, unregulated_gene, signal_regulated_gene)) //signal_regulated_gene is not regulated by regulated_gene  
        {            
            if(find_TFBS_of_A_on_B(genotype, signal_regulated_gene, unregulated_gene)) //unregulated_gene is regulated by signal_regulated_gene, then we found an isolated C1FFL                           
            {
                return_val=1;   
                *fast_TF=genotype->which_protein[signal_regulated_gene];
                *slow_TF=genotype->which_protein[unregulated_gene];
                *fast_TF_gene_id=signal_regulated_gene;
                *slow_TF_gene_id=unregulated_gene;
            }         
        }      
    }    
    return return_val;    
}

static int find_FFL_in_diamond(Genotype *genotype, int gene_A, int gene_B, int *fast_TF, int *slow_TF, int *fast_TF_gene_id, int *slow_TF_gene_id)
{   
    float fast_TF_protein_decay_rate, slow_TF_protein_decay_rate;
    int i;
    int return_val=0;
    
    /*A and B must encode TFs of different families*/
    if(genotype->which_TF_family[genotype->which_protein[gene_A]]!=genotype->which_TF_family[genotype->which_protein[gene_B]])
    { 
        if(!find_TFBS_of_A_on_B(genotype, gene_A, gene_B)) //A does not regulate B
        {            
            if(find_TFBS_of_A_on_B(genotype, gene_B, gene_A)) //B regulates A, then we found an FFL-in-diamond   
            {
                *fast_TF=genotype->which_protein[gene_B];
                *slow_TF=genotype->which_protein[gene_A];
                *fast_TF_gene_id=gene_B;
                *slow_TF_gene_id=gene_A;
                return_val=1; 
            }
            else //a diamond
            {
                /*assuming A is fast TF and B is slow TF*/
                *fast_TF=genotype->which_protein[gene_A];
                *slow_TF=genotype->which_protein[gene_B];
                /*check if the assignment is correct by comparing r_protein_deg*/
                /*calculate geometric mean of r_protein_deg across variants*/
                fast_TF_protein_decay_rate=1.0;
                slow_TF_protein_decay_rate=1.0;
                for(i=0;i<genotype->protein_pool[*fast_TF][0][0];i++)
                    fast_TF_protein_decay_rate*=genotype->protein_decay_rate[genotype->protein_pool[*fast_TF][1][i]];
                fast_TF_protein_decay_rate=pow(fast_TF_protein_decay_rate,1.0/(float)genotype->protein_pool[*fast_TF][0][0]);
                for(i=0;i<genotype->protein_pool[*slow_TF][0][0];i++)
                    slow_TF_protein_decay_rate*=genotype->protein_decay_rate[genotype->protein_pool[*slow_TF][1][i]];
                slow_TF_protein_decay_rate=pow(slow_TF_protein_decay_rate,1.0/(float)genotype->protein_pool[*slow_TF][0][0]);
                /*fast TF should have larger r_protein_deg*/
                *fast_TF=(fast_TF_protein_decay_rate>slow_TF_protein_decay_rate)?genotype->which_protein[gene_A]:genotype->which_protein[gene_B];
                *slow_TF=(fast_TF_protein_decay_rate>slow_TF_protein_decay_rate)?genotype->which_protein[gene_B]:genotype->which_protein[gene_A];
                *fast_TF_gene_id=(fast_TF_protein_decay_rate>slow_TF_protein_decay_rate)?gene_A:gene_B;
                *slow_TF_gene_id=(fast_TF_protein_decay_rate>slow_TF_protein_decay_rate)?gene_B:gene_A;
                return_val=2;
            }
        }     
        else //A regulates B
        {   
            if(!find_TFBS_of_A_on_B(genotype, gene_B, gene_A)) //B does not regulate A, then an FFL-in-diamond       
            {
                *fast_TF=genotype->which_protein[gene_A];
                *slow_TF=genotype->which_protein[gene_B];
                *fast_TF_gene_id=gene_A;
                *slow_TF_gene_id=gene_B;
                return_val=1;             
            }
        }
    }    
    return return_val;    
}

static void determine_motif_logic(Genotype *genotype,                                     
                                    int hindrance[MAX_PROTEINS][MAX_PROTEINS], 
                                    int fast_TF, 
                                    int slow_TF, 
                                    int n_motifs,
                                    char motif_name,
                                    int strong_BS_pos[MAX_PROTEINS][2][50],
                                    int effector_gene_id,
                                    int fast_TF_gene_id, 
                                    int slow_TF_gene_id)
{  
    int start_entry_in_N_motifs;    
    /*where in genotype->N_motifs and in genotype->N_near_AND_gated_motifs to store the results*/
    switch(motif_name)
    {
        case 'D': //direct regulation
            start_entry_in_N_motifs=0;           
            break;
        case 'C': //isolated C1
            start_entry_in_N_motifs=9;         
            break;
        case 'F': //FFL-in-diamond
            start_entry_in_N_motifs=18;          
            break;
        case 'I': //diamond
            start_entry_in_N_motifs=27;  
    }
    
    /*add n_motifs to AND logic*/
    genotype->N_motifs[start_entry_in_N_motifs]+=n_motifs;
    
    /*using hindrance table to determine logic*/
    if(hindrance[fast_TF][slow_TF]) 
    {
        if(hindrance[fast_TF][slow_TF]) // signal alone cannot activate effector
        {
            if(hindrance[slow_TF][slow_TF])
                genotype->N_motifs[start_entry_in_N_motifs+1]+=n_motifs; // effectively no regulation
            else
                genotype->N_motifs[start_entry_in_N_motifs+2]+=n_motifs; // I3
        }   
        else
        {
            if(hindrance[slow_TF][slow_TF])
                genotype->N_motifs[start_entry_in_N_motifs+3]+=n_motifs; //I1
            else
                genotype->N_motifs[start_entry_in_N_motifs+4]+=n_motifs; //impossible
        }
    }
    else
    {
        if(hindrance[fast_TF][fast_TF])
        {
            if(hindrance[slow_TF][slow_TF])
            {
                genotype->N_motifs[start_entry_in_N_motifs+5]+=n_motifs; // AND-gated 
#if SAMPLE_PARAMETERS
                mark_genes_for_sampling(genotype, effector_gene_id, fast_TF, slow_TF, motif_name);               
#endif
#if COUNT_NEAR_AND
                determine_near_AND_logic(genotype, strong_BS_pos, fast_TF, slow_TF, n_motifs, motif_name, 0);              
#endif   
#if PERTURB
                mark_genes_to_be_perturbed(genotype, fast_TF, slow_TF, effector_gene_id, fast_TF_gene_id, slow_TF_gene_id);     
#endif
            }
            else
            {
                genotype->N_motifs[start_entry_in_N_motifs+6]+=n_motifs; // slow tf controlled
#if COUNT_NEAR_AND
                determine_near_AND_logic(genotype, strong_BS_pos, fast_TF, slow_TF, n_motifs, motif_name, 1);              
#endif          
            }
        }   
        else
        {
            if(hindrance[slow_TF][slow_TF])    
            {
                genotype->N_motifs[start_entry_in_N_motifs+7]+=n_motifs; // fast TF control
#if COUNT_NEAR_AND
                determine_near_AND_logic(genotype, strong_BS_pos, fast_TF, slow_TF, n_motifs, motif_name, 2);              
#endif  
            }
            else
            {
                genotype->N_motifs[start_entry_in_N_motifs+8]+=n_motifs; // OR-gated
#if COUNT_NEAR_AND
                determine_near_AND_logic(genotype, strong_BS_pos, fast_TF, slow_TF, n_motifs, motif_name, 3);              
#endif  
            }
        }
    }
}

static void determine_near_AND_logic(Genotype *genotype, 
                                        int strong_BS_pos[MAX_PROTEINS][2][50], 
                                        int fast_TF, 
                                        int slow_TF, 
                                        int n_motifs, 
                                        char motif_name, 
                                        int actual_logic)
{
    int start_entry_in_N_near_AND;
    /*where in genotype->N_near_AND_gated_motifs to store the results*/
    switch(motif_name)
    {
        case 'D':           
            start_entry_in_N_near_AND=0;
            break;
        case 'C':            
            start_entry_in_N_near_AND=0;
            break;
        case 'F':           
            start_entry_in_N_near_AND=4;
            break;
        case 'I':            
            start_entry_in_N_near_AND=8;            
    }
    
    /* For a non-AND gate to be counted as near-AND gate, fast TF and slow TF each 
     * must have at least one strong TFBSs, and the two TFs must form AND gate with
     * only the strong TFBSs */
    if(strong_BS_pos[slow_TF][0][0]!=0 && strong_BS_pos[fast_TF][0][0]!=0) // both TF have strong TFBSs
    {
        //slow TF must have only one strong TFBS or strong TFBSs of slow TF overlap with each other 
        if(strong_BS_pos[slow_TF][0][0]==1 ||strong_BS_pos[slow_TF][1][strong_BS_pos[slow_TF][0][0]-1]-strong_BS_pos[slow_TF][1][0]<CONSENSUS_SEQ_LEN+2*HIND_LENGTH) 
        {
            /*fast TF must have only one strong TFBS or strong TFBSs of slow TF overlap with each other*/
            if(strong_BS_pos[fast_TF][0][0]==1 || strong_BS_pos[fast_TF][1][strong_BS_pos[fast_TF][0][0]-1]-strong_BS_pos[fast_TF][1][0]<CONSENSUS_SEQ_LEN+2*HIND_LENGTH)   
            {
                /* at least one strong TFBS of the fast TF does not overlap with at least one strong TFBSs of the slow TF
                 * i.e. there is enough distance between strong TFBSs of the fast and strong TFBSs of the slow*/
                if(abs(strong_BS_pos[slow_TF][1][0]-strong_BS_pos[fast_TF][1][strong_BS_pos[fast_TF][0][0]-1])>=CONSENSUS_SEQ_LEN+2*HIND_LENGTH ||
                    abs(strong_BS_pos[slow_TF][1][strong_BS_pos[slow_TF][0][0]-1]-strong_BS_pos[fast_TF][1][0])>=CONSENSUS_SEQ_LEN+2*HIND_LENGTH 
                )  
                    genotype->N_near_AND_gated_motifs[start_entry_in_N_near_AND+actual_logic]+=n_motifs;
            }
        }
    }
}

#if PERTURB
/*mark the cis and trans targets of perturbation*/
static void mark_genes_to_be_perturbed(Genotype *genotype, int fast_TF, int slow_TF, int effector_gene_id, int fast_TF_gene_id, int slow_TF_gene_id)
{
    int cis_target, trans_target;
    
    switch(WHICH_CIS_TARGET)
    {
        case 0:
            cis_target=effector_gene_id;
            break;            
        case 1:
            cis_target=fast_TF_gene_id;
            break;
        case 2:
            cis_target=slow_TF_gene_id;
    }
    
    switch(WHICH_TRANS_TARGET)
    {
        case 0:
            trans_target=N_SIGNAL_TF-1;
            break;
        case 1:
            trans_target=fast_TF;
            break;
        case 2:
            trans_target=slow_TF;                
    }
    genotype->cis_target_to_be_perturbed[cis_target]=YES;
    genotype->trans_target_to_be_perturbed[cis_target][trans_target]=YES;
}


static void modify_topology(Genotype *genotype, Genotype *genotype_clone)
{
    int i, j, gene_id;
    
    for(i=0;i<36;i++)
        genotype_clone->N_motifs[i]=genotype->N_motifs[i];
    for(i=0;i<MAX_GENES;i++)
    {
        genotype_clone->cis_target_to_be_perturbed[i]=genotype->cis_target_to_be_perturbed[i];
        for(j=0;j<MAX_PROTEINS;j++)
            genotype_clone->trans_target_to_be_perturbed[i][j]=genotype->trans_target_to_be_perturbed[i][j];
    }
    
    for(gene_id=N_SIGNAL_TF;gene_id < genotype_clone->ngenes;gene_id++)
    {        
#if ADD_TFBS
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
    
    /*make sure there is enough room to add one more TFBS*/
    if (genotype->binding_sites_num[gene_id] + 1 >= genotype->N_allocated_elements) 
    {  
        while(genotype->N_allocated_elements<=genotype->binding_sites_num[gene_id]+1)
            genotype->N_allocated_elements+=100;

        for(i=0;i<MAX_GENES;i++)
        {
            genotype->all_binding_sites[i] = realloc(genotype->all_binding_sites[i], genotype->N_allocated_elements*sizeof(AllTFBindingSites));
            if (!genotype->all_binding_sites[i]) 
            {
#if MAKE_LOG                
                LOG("error in calc_all_binding_sites_copy\n");     
#endif
                exit(-1);                                       
            }     
        }                    
    }  
    
    /*add one TFBSs of the trans-target to cis-target*/
    if(genotype->cis_target_to_be_perturbed[gene_id]==YES)
    {
        for(i=0;i<genotype->nproteins;i++)
        {
            if(genotype->trans_target_to_be_perturbed[gene_id][i]==YES)
            {
                if(ADD_STRONG_TFBS)
                {
                    temp=1.0;                
                    for(j=0;j<genotype->binding_sites_num[gene_id];j++)
                    {
                        if(genotype->all_binding_sites[gene_id][j].tf_id==i) 
                            temp=(temp>genotype->all_binding_sites[gene_id][j].Kd)?genotype->all_binding_sites[gene_id][j].Kd:temp;  // find a strong binding site                      
                    }
                    genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd=temp;
                }
                else
                    genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].Kd = KD2APP_KD*genotype->Kd[i]*pow(NS_Kd/genotype->Kd[i],(float)(CONSENSUS_SEQ_LEN-NMIN)/(CONSENSUS_SEQ_LEN-NMIN+1)); 
                genotype->all_binding_sites[gene_id][genotype->binding_sites_num[gene_id]].tf_id = i;                 
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
  
    for(i=NMIN/2; i < CISREG_LEN-CONSENSUS_SEQ_LEN-NMIN/2; i++) /* scan promoter */
    {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(genotype->all_binding_sites[gene_id][j].BS_pos> i-CONSENSUS_SEQ_LEN-2*HIND_LENGTH)
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
            if(!(genotype->cis_target_to_be_perturbed[gene_id]==YES && genotype->trans_target_to_be_perturbed[gene_id][k]==YES))
            {

                tf_seq=&(genotype->tf_seq[k][0]);
                tf_seq_rc=&(genotype->tf_seq_rc[k][0]);            
                /*find BS on the template strand*/
                match=0;
                for (j=i;j<i+CONSENSUS_SEQ_LEN;j++) 
                    if (cis_seq[j] == tf_seq[j-i]) match++; 

                /*if more than NMIN base pairs are matched*/
                if(match>=NMIN)
                {                              
                    genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                      
                    genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(CONSENSUS_SEQ_LEN-match)/(CONSENSUS_SEQ_LEN-NMIN+1));
                    genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i ; 
                    genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = CONSENSUS_SEQ_LEN-match;             
                    genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                    N_hindered_BS++;              
                    N_binding_sites++;
                    if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                }
                else /*find BS on the non-template strand.*/
                {
                    match_rc=0;
                    for (j=i; j < i+CONSENSUS_SEQ_LEN; j++)                
                        if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
                    if (match_rc >= NMIN)
                    {                   
                        genotype->all_binding_sites[gene_id][N_binding_sites].tf_id = k;                                     
                        genotype->all_binding_sites[gene_id][N_binding_sites].Kd=KD2APP_KD*genotype->Kd[k]*pow(NS_Kd/genotype->Kd[k],(float)(CONSENSUS_SEQ_LEN-match_rc)/(CONSENSUS_SEQ_LEN-NMIN+1));
                        genotype->all_binding_sites[gene_id][N_binding_sites].BS_pos = i;
                        genotype->all_binding_sites[gene_id][N_binding_sites].mis_match = CONSENSUS_SEQ_LEN-match_rc;
                        genotype->all_binding_sites[gene_id][N_binding_sites].N_hindered = N_hindered_BS;
                        N_hindered_BS++;                  
                        N_binding_sites++;                 
                        if(genotype->protein_identity[k]==ACTIVATOR) genotype->N_act_BS[gene_id]++;
                    }
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
        while(j!=0 && genotype->all_binding_sites[gene_id][act_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][act_BS[j][0]].BS_pos<CONSENSUS_SEQ_LEN+2*HIND_LENGTH)j--;
        act_BS[i][1]=act_BS[j][1]+1;
    } 
    rep_BS[0][0]=-1;
    rep_BS[0][1]=0;
    for(i=1;i<N_rep_BS;i++) 
    {
        j=i-1;
        while(j!=0 && genotype->all_binding_sites[gene_id][rep_BS[i][0]].BS_pos - genotype->all_binding_sites[gene_id][rep_BS[j][0]].BS_pos<CONSENSUS_SEQ_LEN+2*HIND_LENGTH)j--;
        rep_BS[i][1]=rep_BS[j][1]+1;
    }
    genotype->max_unhindered_sites[gene_id][1]=act_BS[N_act_BS-1][1];
    genotype->max_unhindered_sites[gene_id][2]=rep_BS[N_rep_BS-1][1];
}
#endif //end of PERTURB mode

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
    for(i=0;i<(replay_N_steps/SAVING_INTERVAL);i++)
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
    for(i=0;i<36;i++)
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
            resident_info->n_effector_genes=resident->ngenes-resident->ntfgenes;
            resident_info->n_act=resident->N_act;
            resident_info->n_rep=resident->N_rep;
         
            resident_info->mut_type=mut_record->mut_type;
            resident_info->which_gene=mut_record->which_gene;
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
        fprintf(fp,"%d %d %c %d %d '%s' %d %a\n",
                mutant_info[i].step,
                mutant_info[i].n_tot_mut,
                mutant_info[i].mut_type,
                mutant_info[i].which_gene,
                mutant_info[i].which_nuc,
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
        for(i=0;i<36;i++)
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
            fprintf(fp,"%c %d %d '%s' %d %a\n", resident_info[i].mut_type,    
                                                resident_info[i].which_gene,
                                                resident_info[i].which_nuc,
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
                    resident_info[i].n_effector_genes,                  
                    resident_info[i].n_act,
                    resident_info[i].n_rep);
        fflush(fp);
        fclose(fp);
    }
}