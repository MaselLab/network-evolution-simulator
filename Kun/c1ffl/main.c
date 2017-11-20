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
    float kdis[NUM_K_DISASSEMBLY]; 
    char RuntimeSumm[32],filename1[32],filename2[32],filename3[32],filename4[32],filename5[32], filename6[32]; 
    unsigned long int seed[6];
    int i,j;
    FILE *fp;
    initialize_growth_rate_parameters();
    /* get the kdis.txt values */
    read_kdisassembly(kdis);
    /* now create and run the population of cells */
    i=2;
    
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
    fp=fopen("saving_point","r");
    if(fp==NULL)
    {
        fp=fopen(RuntimeSumm,"w");
        fprintf(fp,"Max simulation steps=%d, max trials=%d\n",MAX_MUT_STEP,MAX_TRIALS);
        fprintf(fp,"TF_ELEMENT_LEN=%d, NMIN=%d, HIND_LENGTH=%d, MAX_MODE=%d\n",TF_ELEMENT_LEN,NMIN,HIND_LENGTH,MAX_MODE);
        fprintf(fp,"max effector gene number=%d, max tf gene number=%d\n",EFFECTOR_GENES,TFGENES);            
        fprintf(fp,"cost of expression=%f, penalty=%f\n",h,penalty);
        fprintf(fp,"initial TF number=%d, initial ACT number=%d, initial REP number=%d\n",init_TF_genes,init_N_act,init_N_rep);  
        fprintf(fp,"minimal activators to transcribe selection gene: %d\n",min_act_to_transcr_selection_protein);  
        init_run_pop(kdis,RuntimeSumm,filename1,filename2,filename3,filename4,filename5,filename6,seed,0);
    }
    else
    {
        fclose(fp);
        init_run_pop(kdis,RuntimeSumm,filename1,filename2,filename3,filename4,filename5,filename6,seed,1);
    }    
    return 0;
}



