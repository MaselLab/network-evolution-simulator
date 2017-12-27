/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster
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
#include <sys/stat.h>
#include <errno.h>
#include "lib.h"
#include "netsim.h"

int main()
{
    unsigned long int seeds[6];
    int i,seed;
    FILE *fp;	
    /*make rng seed*/
    seed=38;
    for(i=0;i<6;i++)
        seeds[i]=seed;
    /*make filenames*/    
    snprintf(error_file,sizeof(char)*32,"errors_%i.txt",seed);
    snprintf(mutation_file,sizeof(char)*32,"MUT_%i.txt",seed);
    snprintf(output_file,sizeof(char)*32,"output_%i.txt",seed);
    snprintf(RuntimeSumm,sizeof(char)*32,"RuntimeSummary_%i.txt",seed);
    /*continue simulation?*/
    chdir("result"); 
    fp=fopen("saving_point.txt","r");
    if(fp==NULL)
    {
        fp=fopen(RuntimeSumm,"w");
        fprintf(fp,"Max simulation steps=%d, max trials=%d\n",MAX_MUT_STEP,MAX_TRIALS);
        fprintf(fp,"TF_ELEMENT_LEN=%d, NMIN=%d, HIND_LENGTH=%d, MAX_MODE=%d\n",TF_ELEMENT_LEN,NMIN,HIND_LENGTH,MAX_BINDING);
        fprintf(fp,"max effector gene number=%d, max tf gene number=%d\n",EFFECTOR_GENES,TFGENES);            
        fprintf(fp,"c_transl=%f, penalty=%f\n",c_transl,penalty);
        fprintf(fp,"initial TF number=%d, initial ACT number=%d, initial REP number=%d\n",init_TF_genes,init_N_act,init_N_rep);  
        fprintf(fp,"minimal activators to transcribe selection gene: %d\n",min_act_to_transcr_selection_protein);  
        init_run_pop(seeds,0);
    }
    else
    {
        fclose(fp);
        init_run_pop(seeds,1);
    }    
    return 0;
}