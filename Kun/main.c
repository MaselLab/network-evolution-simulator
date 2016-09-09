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

#include "random.h"
#include "lib.h"
#include "netsim.h"



int main()
{
    float kdis[NUM_K_DISASSEMBLY]; 
    char RuntimeSumm[32],filename1[32],filename2[32],filename3[32],filename4[32],filename5[32], filename6[32]; 
    unsigned long int seed[6];
    int i,j;
    FILE *fp;
    initialize_growth_rate_parameters();
    /* get the kdis.txt values */
    read_kdisassembly(kdis);
    /* now create and run the population of cells */
    for(i=4;i<=4;i++)
    {
        for(j=0;j<6;j++)
            seed[j]=i;

        snprintf(filename1,sizeof(char)*32,"output_%i.txt",i);
        snprintf(filename2,sizeof(char)*32,"MUT_%i.txt",i);
        snprintf(filename3,sizeof(char)*32,"error_%i.txt",i);
        snprintf(filename4,sizeof(char)*32,"N_BS_act_%i.txt",i);
        snprintf(filename5,sizeof(char)*32,"N_rep_BS_%i.txt",i);
        snprintf(filename6,sizeof(char)*32,"N_BS_%i.txt",i);
        snprintf(RuntimeSumm,sizeof(char)*32,"RuntimeSummary_%i.txt",i);
        
        chdir("result/test0601/c"); 
        fp=fopen(RuntimeSumm,"w");
        fprintf(fp,"BURN_IN=%d\n",BURN_IN);
        fprintf(fp,"MAX_MUT_STEP=%d\n",MAX_MUT_STEP);
        fprintf(fp,"MAX_MODE=%d\n",MAX_MODE);
        fprintf(fp,"N_replicates=%d\n",N_replicates);
        fprintf(fp,"T-development=%f\n",tdevelopment);
        fprintf(fp,"Duration of burn-in growth rate=%f\n",duration_of_burn_in_growth_rate);
        fprintf(fp,"Environment 1: T-signalA=%f min, T-signalB=%f min, signalA as noise=%d, signalA mismatches=%d\n",env1_t_signalA, env1_t_signalB, env1_signalA_as_noise, env1_signalA_mismatches);
        fprintf(fp,"Environment 2: T-signalA=%f min, T-signalB=%f min, signalA as noise=%d, signalA mismatches=%d\n",env2_t_signalA, env2_t_signalB, env1_signalA_as_noise, env2_signalA_mismatches);
        fprintf(fp,"cost_term=%f, penalty=%f\n",cost_term,penalty);
        fprintf(fp,"initial TF number=%d, initial ACT number=%d, initial REP number=%d\n",init_TF_genes,init_N_act,init_N_rep);         
        fclose(fp);
        
        init_run_pop(kdis,filename1,filename2,filename3,filename4,filename5,filename6,seed);
    }
    return 0;
}
