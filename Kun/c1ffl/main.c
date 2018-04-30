/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#include <stdio.h>
#include <unistd.h>
#include "netsim.h"
#include "RngStream.h"

char error_file[32];
char mutation_file[32];
char RuntimeSumm[32];
char output_file[32];

int main()
{
    /*default output directory*/
    chdir("result");
    
    /*make rng seed*/
    unsigned long int seeds[6];
    int i, seed; 
    RngStream RS_main, RS_parallel[N_THREADS];
    seed=98;
    for(i=0;i<6;i++)
        seeds[i]=seed;
    RngStream_SetPackageSeed(seeds);    
    RS_main=RngStream_CreateStream("Main");
    /*create rng streams*/    
    for(i=0;i<N_THREADS;i++)
        RS_parallel[i]=RngStream_CreateStream(""); 
    
    /*make filenames*/    
    snprintf(error_file,sizeof(char)*32,"errors_%i.txt",seed);
    snprintf(mutation_file,sizeof(char)*32,"MUT_%i.txt",seed);
    snprintf(output_file,sizeof(char)*32,"output_%i.txt",seed);
    snprintf(RuntimeSumm,sizeof(char)*32,"RuntimeSummary_%i.txt",seed);    
    
    /* initialize genotype */    
    int init_TF_genes=6;
    int init_N_act=3;
    int init_N_rep=3;
    int init_effector_gene=1;
    Genotype resident, mutant;    
    initialize_cache(&resident);    
    initialize_cache(&mutant);    
    initialize_genotype(&resident, init_TF_genes, init_N_act, init_N_rep, init_effector_gene, RS_main);
    print_mutatable_parameters(&resident,0);
    mutant.ngenes=resident.ngenes;
    mutant.ntfgenes=resident.ntfgenes;
    mutant.nproteins=resident.nproteins; 
    
    /*Log*/
    FILE *fp;
    fp=fopen(RuntimeSumm,"w");
    fprintf(fp,"Max simulation steps=%d, Burn-in steps=%d, max trials=%d\n",MAX_MUT_STEP, BURN_IN, MAX_TRIALS);        
    fprintf(fp,"max effector gene number=%d, max tf gene number=%d\n",MAX_EFFECTOR_GENES,MAX_TF_GENES);    
    fprintf(fp,"initial TF number=%d, initial ACT number=%d, initial REP number=%d\n",init_TF_genes,init_N_act,init_N_rep); 
    fclose(fp);
    
    /* initial mRNA and protein concentration are 0*/
    int init_mRNA[MAX_GENES]; 
    float init_protein[MAX_GENES]; 
    for(i=N_SIGNAL_TF; i < MAX_GENES; i++) 
    {
        init_mRNA[i]=0;
        init_protein[i]=0.0;
    }    
    
    /* initialize mut_record */
    Mutation mut_record;    
    mut_record.kinetic_diff=0.0;
    mut_record.kinetic_type=-1;
    mut_record.mut_type='\0';
    mut_record.nuc_diff[0]='\0';
    mut_record.nuc_diff[1]='\0';
    mut_record.nuc_diff[2]='\0';
    mut_record.which_gene=-1;
    mut_record.which_nucleotide=-1;
    mut_record.N_hit_bound=0;
    
    /*set selection condition and burn_in condition*/
    Selection selection, burn_in;
    
    
    selection.temporary_DUPLICATION=1.5e-7;
    selection.temporary_SILENCING=1.5e-7;
    selection.temporary_N_effector_genes=MAX_EFFECTOR_GENES;
    selection.temporary_N_tf_genes=MAX_TF_GENES;
    selection.temporary_miu_ACT_TO_INT_RATE=1.57;
    selection.temporary_miu_Kd=-5.0;
    selection.temporary_miu_protein_syn_rate=0.021;
    selection.test1.t_development=89.9;
    selection.test2.t_development=89.9;
    selection.test1.signal_on_strength=1000.0;
    selection.test2.signal_on_strength=1000.0;
    selection.test1.signal_off_strength=0.0;
    selection.test2.signal_off_strength=0.0;
    selection.test1.t_signal_on=200.0;
    selection.test2.t_signal_on=10.0;
    selection.test1.t_signal_off=0.0;
    selection.test2.t_signal_off=200.0;
    selection.test1.initial_effect_of_effector='b';
    selection.test2.initial_effect_of_effector='d';
    selection.test1.fixed_effector_effect=1;
    selection.test2.fixed_effector_effect=1;
    selection.test1_weight=0.33;
    selection.test2_weight=0.67;
    selection.test1.duration_of_burn_in_growth_rate=0.0;
    selection.test2.duration_of_burn_in_growth_rate=0.0;   
    
    burn_in.temporary_DUPLICATION=5.25e-9;
    burn_in.temporary_SILENCING=5.25e-9;
    burn_in.temporary_N_effector_genes=2;
    burn_in.temporary_N_tf_genes=9;
    burn_in.temporary_miu_ACT_TO_INT_RATE=0.762;
    burn_in.temporary_miu_Kd=-7.5;
    burn_in.temporary_miu_protein_syn_rate=0.814;        
    burn_in.test1.t_development=89.9;
    burn_in.test2.t_development=89.9;
    burn_in.test1.signal_on_strength=1000.0;
    burn_in.test2.signal_on_strength=1000.0;
    burn_in.test1.signal_off_strength=0.0;
    burn_in.test2.signal_off_strength=0.0;
    burn_in.test1.t_signal_on=200.0;
    burn_in.test2.t_signal_on=10.0;
    burn_in.test1.t_signal_off=0.0;
    burn_in.test2.t_signal_off=200.0;
    burn_in.test1.initial_effect_of_effector='b';
    burn_in.test2.initial_effect_of_effector='d';
    burn_in.test1.fixed_effector_effect=1;
    burn_in.test2.fixed_effector_effect=1;
    burn_in.test1_weight=0.67;
    burn_in.test2_weight=0.33;   
    burn_in.test1.duration_of_burn_in_growth_rate=0.0;
    burn_in.test2.duration_of_burn_in_growth_rate=0.0;
    
    /*can also use an extra file to provide an irregular signal*/
#if EXTERNAL_SIGNAL
    int j,k;
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
#endif
    
    /*evolve network under selection*/
    init_run_pop(&resident, &mutant, &mut_record, &burn_in, &selection, init_mRNA, init_protein, RS_main, RS_parallel);
    /*just throw mutations to the network*/
//    evolve_neutrally(&resident, &mutant, &mut_record, RS_main); 
    /*replay evolution and output expression dynamics*/
//    int replay_N_steps;
//    replay_N_steps=MAX_MUT_STEP; //default
//    run_plotting(&resident, &mutant, &mut_record, &selection, init_mRNA, init_protein, replay_N_steps,RS_parallel);
    /*modify network topology and calculate fitness of the modified network*/
//    plot_alternative_fitness(&resident, &mutant, &mut_record, &selection, init_mRNA, init_protein,RS_parallel);   
    
    release_memory(&resident, &mutant, &RS_main);
    
    return 0;
}