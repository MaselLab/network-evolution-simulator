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

    i=6;
        for(j=0;j<6;j++)
            seed[j]=i;

        snprintf(filename1,sizeof(char)*32,"output_%i.txt",i);
        snprintf(filename2,sizeof(char)*32,"MUT_%i.txt",i);
        snprintf(filename3,sizeof(char)*32,"error_%i.txt",i);
        snprintf(filename4,sizeof(char)*32,"N_BS_act_%i.txt",i);
        snprintf(filename5,sizeof(char)*32,"N_rep_BS_%i.txt",i);
        snprintf(filename6,sizeof(char)*32,"N_BS_%i.txt",i);
        snprintf(RuntimeSumm,sizeof(char)*32,"RuntimeSummary_%i.txt",i);
        
        chdir("result"); 
        fp=fopen(RuntimeSumm,"w");
        fprintf(fp,"MAX_MODE=%d\n",MAX_MODE);
        fprintf(fp,"MAX gene copies per protein=%d\n",MAX_COPIES);
        fprintf(fp,"penalty per copies=%f\n",penalty_of_extra_copies);        
        if(RdcPdup==1)
        {
            fprintf(fp,"Reducing probability of duplcation once copy number exceeds upper limit!\n");
            fprintf(fp,"fold reduction in probabiliy of duplicaton=%f\n",reduction_in_P_dup);
        }
        else
        {
            fprintf(fp,"Reducing growth rate once copy number exceeds upper limit!\n");
            fprintf(fp,"penalty per copies=%f(percentage of maximal growth rate)\n",penalty_of_extra_copies);
        }
        fprintf(fp,"max gene number=%d, max tf gene number=%d,max protein number=%d\n",NGENES,TFGENES,NPROTEINS);            
        fprintf(fp,"cost_term=%f, penalty=%f\n",cost_term,penalty);
        fprintf(fp,"initial TF number=%d, initial ACT number=%d, initial REP number=%d\n",init_TF_genes,init_N_act,init_N_rep);  
        fprintf(fp,"minimal activators to transcribe selection gene: %d\n",min_act_to_transcr_selection_protein);       
        
        init_run_pop(kdis,RuntimeSumm,filename1,filename2,filename3,filename4,filename5,filename6,seed);
    
    return 0;
}
