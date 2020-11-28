/*
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
	seed=81; //set random number here
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
    int init_N_output_act=1; //number of effector genes that are activator at initialization
    int init_N_output_rep=0;
    int init_N_non_output_act=3; //number of non-effector genes that are activator at initialization
    int init_N_non_output_rep=3;  
    Genotype resident, mutant;    
    initialize_cache(&resident);    
    initialize_cache(&mutant);    
    initialize_genotype(&resident, init_N_output_act, init_N_output_rep, init_N_non_output_act, init_N_non_output_rep, RS_main);
    print_mutatable_parameters(&resident,0);
    mutant.ngenes=resident.ngenes;
    mutant.n_output_genes=resident.n_output_genes;    
    mutant.nproteins=resident.nproteins;
    mutant.N_node_families=resident.N_node_families;
    
    /*Log*/
    FILE *fp, *saving_point;
    saving_point=fopen("saving_point.txt","r");
    if(saving_point==NULL) //saving point not found, therefore we start from the beginning
    {
        fp=fopen(setup_summary,"w");
        fprintf(fp,"Rng seed: %d\n",seed);
        fprintf(fp,"Max effector gene number=%d, max gene number=%d\n",MAX_OUTPUT_GENES,MAX_GENES);    
        fprintf(fp,"initial output ACT number=%d, initial output REP number=%d\n",init_N_output_act,init_N_output_rep);   
        fprintf(fp,"initial non-output ACT number=%d, initial non-output REP number=%d\n\n\n",init_N_non_output_act,init_N_non_output_rep);   
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
    mut_record.which_protein=-1;
    mut_record.N_hit_bound=0;
    
    /*set selection condition */
    Selection selection;     
    
    /*set duration of developmental simulation*/
    selection.env1.t_development=360.1; //minutes    
    
    /*burn-in developmental simulation (this isn't burn-in evolution) will ignore the fitness in the first couple minutes*/    
    selection.env1.avg_duration_of_burn_in_growth_rate=10.0;  
    selection.env1.max_duration_of_burn_in_growth_rate=30.0; 
    
    /*Signal strength under the two stages*/ 
    selection.env1.signal_strength_stage1=100.0;    
    selection.env1.signal_strength_stage2=1000.0;
    selection.env1.t_stage1=120.0; //the duration of stage1
    selection.env1.t_stage2=2000.0; //the duration of stage2
    
    /*default parameter about signal dynamics -- these are useless when select for a pulse generator, but don't change them */
    selection.env1.signal_on_aft_burn_in=0; //turn on signal after burn-in growth  
    selection.env1.initial_effect_of_effector='d'; //effector is beneficial in the beginning of env1 
    selection.env1.fixed_effector_effect=1; //the effector maintains the initial fitness effect. Making it 0 allows effector to be beneficial when signal is on and deleterious when signal is off
        
    /*selection condition of pulse generation*/
    selection.env1.opt_peak_response=10000.0; //optimal peak expression level of the effector
    selection.env1.effector_level_before_stage2=0.1;  // target effector expression level (as a fraction of the peak level) at the end of stage 1. 
    selection.env1.fitness_decay_constant=0.693; // this is the sigma squared in Eq. 1
    selection.env1.min_reduction_relative_to_peak=0.2;  // minimal reduction in the effector expression level (as a fraction of the peak level) at the end of simulation  
    selection.env1.window_size=5; // use a 5-minute window to average the expression level of the effector
    selection.env1.max_t_mid=60.0; // this is the t_saturate in Eq. 3
    selection.env1.w1=1.0; //weight of each fitness component
    selection.env1.w2=1.0;
    selection.env1.w3=1.0;
    selection.env1.w4=1.0;
    selection.env1_weight=1.0;
    selection.env2_weight=0.0;   // there is not a second selection condition, so set its weight to zero. Don't change it
    selection.env1.is_burn_in=0;  // whether env1 is burn-in selection condition. Don't change it.    
    
    /*restrict what motif can evolve*/
    selection.effector_is_TF=1; // whether effector is a TF. 1 means NFBLs can evolve
    resident.flag_effector_is_TF=1; // make sure also change the resident genotype    
    selection.signal_ctrl_repressor=1; // whether the signal can regulate repressor. 1 means I1FFLs can evolve
    resident.flag_signal_ctrl_repressor=1;
    
    /*Maximum evolutionary steps*/
    selection.MAX_STEPS=50000; 
    
    /*Default values of mutational parameters. Don't change it*/
    selection.temporary_DUPLICATION=1.5e-7;
    selection.temporary_SILENCING=2.25e-7;
    selection.temporary_miu_ACT_TO_INT_RATE=1.57;
    selection.temporary_miu_Kd=-5.0;
    selection.temporary_miu_protein_syn_rate=0.021; 
    
    /*record selection condition*/
    saving_point=fopen("saving_point.txt","r");
    if(saving_point==NULL) //saving point not found, therefore we start from the beginning
    {
        fp=fopen(setup_summary,"a+");
        fprintf(fp,"Selection condition:\n");
        fprintf(fp,"steps dup_rate del_rate miu_A2I miu_Kd miu_protein_syn\n");
        fprintf(fp,"%d %.2e %.2e %.3f %.3f %.3f\n\n",
                selection.MAX_STEPS,
                selection.temporary_DUPLICATION,
                selection.temporary_SILENCING,             
                selection.temporary_miu_ACT_TO_INT_RATE,
                selection.temporary_miu_Kd,
                selection.temporary_miu_protein_syn_rate);
        fprintf(fp,"dev_time t_stage1 t_stage2 s_sig_stage1 s_sig_stage2\n");
        fprintf(fp,"%.2f %.1f %.1f %.1f %.1f\n",          
                selection.env1.t_development,               
                selection.env1.t_stage1,
                selection.env1.t_stage2,
                selection.env1.signal_strength_stage1,
                selection.env1.signal_strength_stage2);  
        fprintf(fp,"P_opt sigma2 lvl_before_peak lvl_after_peak t_mid window_size\n");
        fprintf(fp,"%.1f %.3f %.1f %.1f %.1f %d\n",          
                selection.env1.opt_peak_response, 
                selection.env1.effector_level_before_stage2, 
                selection.env1.fitness_decay_constant,
                selection.env1.min_reduction_relative_to_peak, 
                selection.env1.max_t_mid,
                selection.env1.window_size);
        fprintf(fp,"weight_f1 weight_f2 weight_f3 weight_f4\n");
        fprintf(fp,"%.1f %.1f %.1f %.1f\n",          
                selection.env1.w1,
                selection.env1.w2,
                selection.env1.w3,
                selection.env1.w4);
        fclose(fp);      
    }
    else
        fclose(saving_point);  
    
    /*setting burn-in evolution is similar to setting selection condition*/
    Selection burn_in;   
    burn_in.MAX_STEPS=0; // Do not burn-in evolution when select for pulse generation
       
    /*Different simulation mode can be enabled from netsim.h*/
#if NEUTRAL //netural evolution -- USELESS when select for pulse generation
    evolve_neutrally(&resident, &mutant, &mut_record, &burn_in, &selection, RS_main); 
#elif PHENOTYPE //output timecourse of expression levels of genes -- USELESS when select for pulse generation    
    show_phenotype(&resident, &mut_record, &selection, init_mRNA, init_protein, 50000, RS_parallel);
#elif PERTURB //pertubation analysis
    modify_network(&resident, &mutant, &mut_record, &selection, &burn_in, init_mRNA, init_protein, RS_parallel); 
#else //the default is to evolve a network under selection conditions specfied above   
    evolve_under_selection(&resident, &mutant, &mut_record, &burn_in, &selection, init_mRNA, init_protein, RS_main, RS_parallel);
#endif        
    release_memory(&resident, &mutant, &RS_main, RS_parallel);    
    return 0;
}
