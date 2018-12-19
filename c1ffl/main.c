/*
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007-2018 Arizona Board of Regents (University of Arizona)
 */
#include <stdio.h>
#include <unistd.h>
#include "netsim.h"

int main()
{
    unsigned long int seeds[6];
    int i,seed;
    int buffer1,buffer2;
    FILE *fp1, *fp2;	
    /*make rng seed*/
	seed=54;
    for(i=0;i<6;i++)
        seeds[i]=seed;
    /*make filenames*/    
    snprintf(error_file,sizeof(char)*32,"errors_%i.txt",seed);
    snprintf(mutation_file,sizeof(char)*32,"accepted_mutations_%i.txt",seed);
    snprintf(output_file,sizeof(char)*32,"evo_summary_%i.txt",seed);
    snprintf(RuntimeSumm,sizeof(char)*32,"sim_setup_%i.txt",seed);
    /*continue simulation?*/
    chdir("result"); 
    fp1=fopen("saving_point.txt","r");
    if(fp1==NULL)
    {
        fp2=fopen(RuntimeSumm,"w");
        fprintf(fp2,"Max simulation steps=%d, max trials=%d\n",MAX_MUT_STEP,MAX_TRIALS);
        fprintf(fp2,"TF_ELEMENT_LEN=%d, NMIN=%d, HIND_LENGTH=%d, MAX_MODE=%d\n",TF_ELEMENT_LEN,NMIN,HIND_LENGTH,MAX_BINDING);
        fprintf(fp2,"max effector gene number=%d, max tf gene number=%d\n",MAX_EFFECTOR_GENES,MAX_TF_GENES);            
        fprintf(fp2,"c_transl=%f, penalty=%f\n",c_transl,penalty);
        fprintf(fp2,"initial TF number=%d, initial ACT number=%d, initial REP number=%d\n",init_TF_genes,init_N_act,init_N_rep);  
        fprintf(fp2,"minimal activators to transcribe selection gene: %d\n",min_N_activator_to_transc_selection_protein);  
        fclose(fp2);
        init_run_pop(seeds,0);
    }
    else
    {
        fscanf(fp1,"%d %d",&buffer1,&buffer2);
        fclose(fp1);
        fp2=fopen(RuntimeSumm,"a+");
        fprintf(fp2,"Continue simulation from step %d\n",buffer1);
        fclose(fp2);
        init_run_pop(seeds,1);
    }    
    return 0;
}