/*
 * Simulator of yeast transcriptional regulatory network evolution
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2007-2018 Arizona Board of Regents (University of Arizona)
 */
#include <stdio.h>
#include <unistd.h>
#include "netsim.h"
#include "RngStream.h"
#include "lib.h"

/*output filenames*/
char mutation_file[32];
char setup_summary[32];
char evo_summary[32];

int main()
{
    /*default output directory*/  
    chdir("result");
    
    /*make rng seed*/
    unsigned long int seeds[6];
    int i, seed; 
    RngStream RS_main, RS_parallel[N_THREADS];
    seed=5;
    for(i=0;i<6;i++)
        seeds[i]=seed;
    RngStream_SetPackageSeed(seeds);    
    RS_main=RngStream_CreateStream("Main");
    /*create rng streams*/    
    for(i=0;i<N_THREADS;i++)
        RS_parallel[i]=RngStream_CreateStream(""); 
    
    /*make filenames*/     
    snprintf(mutation_file,sizeof(char)*32,"accepted_mutation_%i.txt",seed);
    snprintf(evo_summary,sizeof(char)*32,"evo_summary_%i.txt",seed);
    snprintf(setup_summary,sizeof(char)*32,"sim_setup_%i.txt",seed);    
    
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
    FILE *fp, *saving_point;
    saving_point=fopen("saving_point.txt","r");
    if(saving_point==NULL) //saving point not found, therefore we start from the beginning
    {
        fp=fopen(setup_summary,"w");
        fprintf(fp,"Rng seed: %d\n",seed);
        fprintf(fp,"Max effector gene number=%d, max tf gene number=%d\n",MAX_EFFECTOR_GENES,MAX_TF_GENES);    
        fprintf(fp,"initial TF number=%d, initial ACT number=%d, initial REP number=%d\n\n\n",init_TF_genes,init_N_act,init_N_rep);   
        fclose(fp);
    }
    else
        fclose(saving_point);
    
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
    
    /*set selection condition */
    Selection selection;     
    
    /*set duration of developmental simulation*/
    selection.env1.t_development=89.9; //minutes
    selection.env2.t_development=89.9;
    
    /*The code below creates signals under env1 and env2. The signal is consecutive "ONs" and "OFFs" and always starts with "ON". 
     *Use t_signal_on and t_signal_off to specify the duration of "on" and "off"*/ 
    selection.env1.signal_on_strength=1000.0;    
    selection.env1.signal_off_strength=0.0;
    selection.env2.signal_on_strength=1000.0;
    selection.env2.signal_off_strength=0.0;
    selection.env1.t_signal_on=200.0; //the signal starts with "on" and last for 200 minutes, longer than the duration of developmental simulation, which means the signal is effective constant "on" 
    selection.env1.t_signal_off=0.0; 
    selection.env2.t_signal_on=10.0; //the signal is "on" for the first 10 minutes in a developmental simulation of env2     
    selection.env2.t_signal_off=200.0; //after being "on" for the initial 10 minutes, the signal is off for 200 minutes, i.e. remains off in env2.
    
    /*Alternatively, an irregular signal can be specified with an extra file*/    
#if IRREG_SIGNAL
    int j,k;  
    /*The signal will be stored in signal_profile_matrix[A][B][C] (see netsim.h). 
     *1. The signal is specified by signal strength at each of the C time points. C is 90 by default, corresponding 
     *to the 90 minutes of simulation time.
     *2. B is the number of signals that must be provided; each developmental simulation will chooses one from the B signals.
     *The signals don't have to be different. The default value of B is 100.  
     *3. A is the number of parallel threads; each thread needs a copy of the 100 signals.   
    */  
    /*This is the default external signal profile*/
    /*The file should have 90*100 rows. Every 90 rows are the strength of a signal at each time point.*/
    fp=fopen("signal.txt","r"); 
    
    /*load signal*/
    if(fp!=NULL)
    {        
        for(j=0;j<100;j++)
        {
            for(k=0;k<90;k++)
            {           
                fscanf(fp,"%f\n",&(signal_profile_matrix[0][j][k]));
                for(i=1;i<N_THREADS;i++)
                    signal_profile_matrix[i][j][k]=signal_profile_matrix[0][j][k];                    
            }
        }
        fclose(fp);
    } 
#endif    
    
    /*Set when the effector is beneficial and when it is deleterious*/
    selection.env1.initial_effect_of_effector='b'; //effector is beneficial in the beginning of env1
    selection.env2.initial_effect_of_effector='d'; //deleterious in the beginning of env2 
    selection.env1.fixed_effector_effect=1; //the effector maintains the initial fitness effect. Making it 0 allows effector to be beneficial when signal is on and deleterious when signal is off
    selection.env2.fixed_effector_effect=1;
    
    /*how to average fitness under env1 and fitness under env2*/
    selection.env1_weight=0.33;
    selection.env2_weight=0.67;
    
    /*burn-in developmental simulation (this isn't burn-in evolution) will ignore the fitness in the first couple minutes*/    
    selection.env1.duration_of_burn_in_growth_rate=0.0; //zero means do not burn-in developmental simulation
    selection.env2.duration_of_burn_in_growth_rate=0.0;  
    
    /*Maximum evolutionary steps*/
    selection.MAX_STEPS=50000; 
    
    /*Default values of */
    selection.temporary_DUPLICATION=1.5e-7;
    selection.temporary_SILENCING=1.5e-7;
    selection.temporary_N_effector_genes=MAX_EFFECTOR_GENES;
    selection.temporary_N_tf_genes=MAX_TF_GENES;
    selection.temporary_miu_ACT_TO_INT_RATE=1.57;
    selection.temporary_miu_Kd=-5.0;
    selection.temporary_miu_protein_syn_rate=0.021; 
    
    /*record selection condition*/
    saving_point=fopen("saving_point.txt","r");
    if(saving_point==NULL) //saving point not found, therefore we start from the beginning
    {
        fp=fopen(setup_summary,"a+");
        fprintf(fp,"Selection condition:\n");
        fprintf(fp,"steps dup_rate del_rate max_N_tf max_N_effector miu_A2I miu_Kd miu_protein_syn\n");
        fprintf(fp,"%d %.2e %.2e %d %d %.3f %.3f %.3f\n\n",
                selection.MAX_STEPS,
                selection.temporary_DUPLICATION,
                selection.temporary_SILENCING,
                selection.temporary_N_tf_genes,
                selection.temporary_N_effector_genes,
                selection.temporary_miu_ACT_TO_INT_RATE,
                selection.temporary_miu_Kd,
                selection.temporary_miu_protein_syn_rate);
        fprintf(fp,"env weight dev_time burn_in_growth t_sig_on t_sig_off s_sig_on s_sig_off init_effect fixed_effect\n");
        fprintf(fp,"1 %.2f %.1f %.1f %.1f %.1f %.1f %.1f %c %d\n",
                selection.env1_weight,
                selection.env1.t_development,
                selection.env1.duration_of_burn_in_growth_rate,
                selection.env1.t_signal_on,
                selection.env1.t_signal_off,
                selection.env1.signal_on_strength,
                selection.env1.signal_off_strength,
                selection.env1.initial_effect_of_effector,
                selection.env1.fixed_effector_effect);
        fprintf(fp,"2 %.2f %.1f %.1f %.1f %.1f %.1f %.1f %c %d\n\n\n",
                selection.env2_weight,
                selection.env2.t_development,
                selection.env2.duration_of_burn_in_growth_rate,
                selection.env2.t_signal_on,
                selection.env2.t_signal_off,
                selection.env2.signal_on_strength,
                selection.env2.signal_off_strength,
                selection.env2.initial_effect_of_effector,
                selection.env2.fixed_effector_effect);   
        fclose(fp);
    }
    else
        fclose(saving_point);  
    
    /*setting burn-in evolution is similar to setting selection condition*/
    Selection burn_in;         
    burn_in.env1.t_development=89.9;
    burn_in.env2.t_development=89.9;
    burn_in.env1.signal_on_strength=1000.0;  
    burn_in.env1.signal_off_strength=0.0;
    burn_in.env2.signal_on_strength=1000.0;
    burn_in.env2.signal_off_strength=0.0;
    burn_in.env1.t_signal_on=200.0;    
    burn_in.env1.t_signal_off=0.0;
    burn_in.env2.t_signal_on=10.0;
    burn_in.env2.t_signal_off=200.0;
    burn_in.env1.initial_effect_of_effector='b';
    burn_in.env2.initial_effect_of_effector='d';
    burn_in.env1.fixed_effector_effect=1;
    burn_in.env2.fixed_effector_effect=1;
    burn_in.env1_weight=0.67;
    burn_in.env2_weight=0.33;   
    burn_in.env1.duration_of_burn_in_growth_rate=0.0;
    burn_in.env2.duration_of_burn_in_growth_rate=0.0;
    burn_in.MAX_STEPS=0; // set it to 1000 if DIRECT_REG=0    
    burn_in.temporary_DUPLICATION=5.25e-9;
    burn_in.temporary_SILENCING=5.25e-9;
    burn_in.temporary_N_effector_genes=2;
    burn_in.temporary_N_tf_genes=9;
    burn_in.temporary_miu_ACT_TO_INT_RATE=0.762;
    burn_in.temporary_miu_Kd=-7.5;
    burn_in.temporary_miu_protein_syn_rate=0.814;   
    
    saving_point=fopen("saving_point.txt","r");
    if(saving_point==NULL) //saving point not found, therefore we start from the beginning
    {
        if(burn_in.MAX_STEPS!=0)
        {
            fp=fopen(setup_summary,"a+");
            fprintf(fp,"Burn-in evolution enabled!");
            fprintf(fp,"steps dup_rate del_rate max_N_tf max_N_effector miu_A2I miu_Kd miu_protein_syn\n");
            fprintf(fp,"%d %.2e %.2e %d %d %.3f %.3f %.3f\n\n",
                    burn_in.MAX_STEPS,
                    burn_in.temporary_DUPLICATION,
                    burn_in.temporary_SILENCING,
                    burn_in.temporary_N_tf_genes,
                    burn_in.temporary_N_effector_genes,
                    burn_in.temporary_miu_ACT_TO_INT_RATE,
                    burn_in.temporary_miu_Kd,
                    burn_in.temporary_miu_protein_syn_rate);
            fprintf(fp,"env weight dev_time burn_in_growth t_sig_on t_sig_off s_sig_on s_sig_off init_effect fixed_effect\n");
            fprintf(fp,"1 %.2f %.1f %.1f %.1f %.1f %.1f %.1f %c %d\n",
                    burn_in.env1_weight,
                    burn_in.env1.t_development,
                    burn_in.env1.duration_of_burn_in_growth_rate,
                    burn_in.env1.t_signal_on,
                    burn_in.env1.t_signal_off,
                    burn_in.env1.signal_on_strength,
                    burn_in.env1.signal_off_strength,
                    burn_in.env1.initial_effect_of_effector,
                    burn_in.env1.fixed_effector_effect);
            fprintf(fp,"2 %.2f %.1f %.1f %.1f %.1f %.1f %.1f %c %d\n\n\n",
                    burn_in.env2_weight,
                    burn_in.env2.t_development,
                    burn_in.env2.duration_of_burn_in_growth_rate,
                    burn_in.env2.t_signal_on,
                    burn_in.env2.t_signal_off,
                    burn_in.env2.signal_on_strength,
                    burn_in.env2.signal_off_strength,
                    burn_in.env2.initial_effect_of_effector,
                    burn_in.env2.fixed_effector_effect);   
            fclose(fp);
        }    
    }
    else
        fclose(saving_point);  
    
    /*Different simulation mode can be enabled from netsim.h*/
#if NEUTRAL //netural evolution
    evolve_neutrally(&resident, &mutant, &mut_record, &burn_in, &selection, RS_main); 
#elif PHENOTYPE //output timecourse of expression levels of genes    
    show_phenotype(&resident, &mutant, &mut_record, &selection, init_mRNA, init_protein, RS_parallel);
#elif PERTURB //pertubation analysis
    perturbation_analysis(&resident, &mutant, &mut_record, &selection, init_mRNA, init_protein, RS_parallel); 
#else //the default is to evolve a network under selection conditions specfied above   
    evolve_under_selection(&resident, &mutant, &mut_record, &burn_in, &selection, init_mRNA, init_protein, RS_main, RS_parallel);
#endif        
    release_memory(&resident, &mutant, &RS_main, RS_parallel);    
    return 0;
}