/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 * Authors: Joanna Masel, Alex Lancaster, Jasmin Uribe
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
#include <omp.h>

/* local includes */
#include "random.h"
#include "lib.h"
#include "netsim.h"

#define UPDATE_ALL 1
#define NO_KON_UPDATE 0

const int MAXELEMENTS=500*MAX_COPIES; 
/* start by allocating maxelements when initializing a genotype, double as needed, reduce at end */
const int MAXBOUND=500*MAX_COPIES;
const int MAXALLOC=10;
const int NMIN=4;
const float KRNA=618.0;
const float TTRANSLATION=1.0;
const float TTRANSCRIPTION=1.0;
const float PROB_ACTIVATING=0.62;
const float TRANSCRIPTINIT=8.5; /* replace betaon and betaoff */
const float DEACETYLATE=0.462;
const float ACETYLATE=0.1155;
const float PICASSEMBLY=0.0277;
const float STARTNUCLEUS=0.1;
const float KR=10.0;                  /* don't put this less than 1, weird things happen to koff calculation */
const float GASCONSTANT=8.31447;
//const float COOPERATIVITY=1.0;        /* dGibbs, relative to 1 additional specific nt */
//const float COOPERATIVE_DISTANCE=11;  /* distance co-operativity operates (changed from 20) */ 
const float NUMSITESINGENOME = 1.3e+6; /* updated from 1.8e+6 */

const float mN = 3.3e-10 * 5e6;    /* Lynch et al. (2008), Tsai et al. (2008)  */
const float rate_of_event_x= 0.1;
const float Ne = 0.5e9;            /* effective population size of yeast*/

const int avg_protein_conc = 12000;
const float penalty = 2.0e-7;

const float SUBSTITUTION = 0.33e-9; /* susbstitution rate per site per cell division*/
const float INDEL = 0.12e-9;      /* indel rate per site per cell division */
const float DUPLICATION = 8.0e-9;   /* single gene duplication rate per gene per cell division (using 120min) */
const float SILENCING = 5.7e-9;      /* per gene per cell division (120min)*/
const float MUTKINETIC = 0.1;         /* assuming 10% of the subs and indel will change kinetic rates and binding seq */
const float max_inset=3.0;
const float max_delet=3.0;
/* below are default options, can be changed on command line*/

float kon=1e-7;              /* lower value is so things run faster */
                             /* actual value should be kon=0.2225 is
                                based on 1 molecule taking
                                240seconds=4 minutes and 89% of the
                                proteins being in the nucleus*/
const float Pact_scaling = 10.0;
//float kon_after_burnin=1e-4; /* lower value value after burn is so things run faster */
//float Koff[TF_ELEMENT_LEN-NMIN+1];

int burn_in = 1;             /* disable burn-in by default */

const float CELL_CYCLE=60.0; /* yeast cell cylce duration in minutes */
float tdevelopment = 30.0;/* default  development time: can be changed at runtime */
float timemax = -1.0;      /* set an upper limit to development time (default to -1.0=no limit) */
int current_ploidy = 1;    /* ploidy can be changed at run-time: 1 = haploid, 2 = diploid */
int output = 0;
long master_seed = -284671;  /* NOTE!!! parallel computation will mess up a shared seed*/
int dummyrun = 10;          /* used to change seed */

float growth_rate_scaling = 1.0; /* set default growth rate scaling factor */
float duration_env0 = 59.9; // in minutes
float duration_env1 = 59.9;
int N_replicates=150;
int max_sampling=150;
float cost_term=0.2;
/* end default options */ 


/* protein aging term: used when c=c'+g=0, set to 1e-4 < mean-3*sd of
   Belle et al. (2006) and small with respect to most growth rates  */

float protein_aging = 1e-4;

/* file output parameters */
char *output_directory = "output";   /* default output directory */
int verbose = 0;                     /* don't log verbosely by default */ 

/* initialize the growth rate parameters: 
 * do computations here so that we can easily change the scaling factor and Pp */
void initialize_growth_rate_parameters() {
  float hc, gpeak_a, gpeak_b, Ltf;
  gpeak_a = 0.005776*growth_rate_scaling;  /* in min^-1 based on doubling time of 120 min: ln(2)/(120 min)=0.005776 */
  gpeak_b = 0.003*growth_rate_scaling;
  Pp_a = 20000;
  Pp_b = 30000;              /* mean gene expression of all proteins is 12064.28 */
  Ltf= 1418;               /* mean gene expression of only TFs is 1418 */
  hc = (gpeak_a/avg_protein_conc)*(1-(log(2-2*cost_term)/log(2)));      /* in min^-1 cost of doubling gene expression, based on Wagner (2005) 
                                                   * using {s=0.2, N=500} matches {s=10^-5, N=10^7} combination (both Ns=100) */
  h = hc/0.023;            /* using c=0.023/min from mean of distribution from Belle et al (2006)*/
  gmax_a = gpeak_a + hc*(Pp_a+(TFGENES*Ltf));    /* compute the gmax coefficient based on gpeak and other parameters */
  gmax_b = gpeak_b + hc*(Pp_b+(TFGENES*Ltf));
}

char set_base_pair(float x) {
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

void initialize_sequence(char *Seq, 
                         int len,
                         int ploidy,
                         int num_elements)
{
  float x;
  int i;
  int current_element = len/(num_elements*ploidy);
  int first, second, third, fourth;

//  LOG_VERBOSE_NOCELLID("len=%d, NGENES=%d, ploidy=%d, current_element=%d\n", len, NGENES, ploidy, current_element); 
  for (i=0; i<len/ploidy; i++) {
    first = (i / current_element)*ploidy*current_element + i % current_element;
//    second = first + current_element;
//    third = second + current_element;
//    fourth = third + current_element;
//    LOG_VERBOSE_NOCELLID("first=%d\n", first); 
    
    x = ran1(&master_seed);    
    
    Seq[first] = set_base_pair(x);
    /* clone the randomly chosen sequence for all other sites */
//    Seq[second] = Seq[first];
//    Seq[third] = Seq[first];
//    Seq[fourth] = Seq[first];
  }
//  LOG_VERBOSE_NOCELLID("length: %d, sequence is %s\n", strlen(Seq), Seq);
}

void print_genotype(Genotype *genotype, int genotype_id) {
//  int i, p;
//
//  printf("[genotype %03d] hind_pos: ", genotype_id);
//  for (i=0; i < TFGENES; i++) {
//    printf("%d ", genotype->hindrance_positions[i]);
//  }
//  printf("\n");
//
//  for (i=0; i < NGENES; i++) {
//    printf("[genotype %03d gene %02d] ", genotype_id, i);
////    printf("repl=%g, ", genotype->replication_time[i]);
//    printf("mRNAdy=%g, ", genotype->mRNAdecay[i]);
//    printf("protdy=%g, ", genotype->proteindecay[i]);
//    printf("transl=%g, ", genotype->translation[i]);
//
//    for (p=0; p < MAX_COPIES; p++) {
//      printf("act[%d]=%d, ", p, genotype->activating[i][p]);
//      printf("PICd[%d]=%g ", p, genotype->pic_disassembly[i][p]);
//    }
//    printf("\n");
//  }
}


void print_all_binding_sites(int copies[NGENES],
                             AllTFBindingSites *all_binding_sites, 
                             int numElements,
                             char tf_seq[TFGENES][MAX_COPIES][TF_ELEMENT_LEN],
                             char cisreg_seq[NGENES][MAX_COPIES][CISREG_LEN]
//                             int site_id_pos[NGENES][MAX_COPIES][2]
							 )
{
//  int i, j;
//
//  // TODO: tidy up logic, a bit messy now, but works
//  /* loop through the maximum of either the number of genes or number
//     of TFs, as they can now be different */
//  int max_elements = NGENES >= TFGENES ? NGENES : TFGENES;
//
//  for (i=0; i < max_elements; i++) {
//    j=0; 
//    /* if copies for this element isn't set, still print out the TFs */
//    while (i >= NGENES || j < copies[i]) {
//      if (i < TFGENES) 
//        printf("TF sequence gene %2d (copy %d): %.*s\n", i, j, TF_ELEMENT_LEN, tf_seq[i][j]);
//      else
//        printf("            gene %2d (copy %d) does not encode a TF\n", i, j);
//      if (i < NGENES) {
//        printf("cis-reg     gene %2d (copy %d): %.*s\n", i, j, CISREG_LEN, cisreg_seq[i][j]);
////        printf("ID range    gene %2d (copy %d): [%3d, %3d]\n", i, j, site_id_pos[i][j][0], site_id_pos[i][j][1]);
//        printf("\n");
//      } else {
//        printf("            gene %2d (copy %d): no cis-reg gene here\n", i, j);
//        printf("\n");
//        break;       /* run out of cisreg genes skip to next element */
//      }
//      j++;
//    }
//  } 
//
//  printf("numElements: %3d\n", numElements);
//  
//  for (i=0; i < numElements; i++) {
//    printf("binding site %3d:\n", i);
//    printf("       cis-reg region: %3d", all_binding_sites[i].cisreg_id);
//    printf("         cis-reg copy: %3d", all_binding_sites[i].gene_copy);
//    printf(" (sequence %.*s)\n", CISREG_LEN, cisreg_seq[all_binding_sites[i].cisreg_id][all_binding_sites[i].gene_copy]);
//    printf(" transcription-factor: %3d", all_binding_sites[i].tf_id);
//    printf(" (sequence: %.*s)\n", TF_ELEMENT_LEN, tf_seq[all_binding_sites[i].tf_id][all_binding_sites[i].gene_copy]); 
//    printf("  L-edge of %2dbp hind: %3d\n", HIND_LENGTH, all_binding_sites[i].left_edge_pos);        
//    printf("  Hind offset position: %3d\n", all_binding_sites[i].hind_pos); 
//    printf("               strand: %3d\n", all_binding_sites[i].strand);
//    printf("         Hamming dist: %3d\n", all_binding_sites[i].hamming_dist); 
//  }
}

#if 0
void print_rounding(CellState *state, GillespieRates *rates, float t)
{
  fprintf(fp_rounding[state->cell_id], "%g %d %d %d %d %d %d %d\n", 
          t, rates->koff_operations, rates->transport_operations, rates->mRNAdecay_operations, 
          rates->pic_disassembly_operations, rates->salphc_operations, rates->max_salphc_operations, rates->min_salphc_operations);
}
#endif

void initialize_genotype_fixed(Genotype *genotype, 
                               float kdis[])
{
    int i, j, p;

    genotype->N_act=0;
    genotype->N_rep=0;

//    LOG_NOCELLID("[genotype %03d] activators vs repressors ", genotype_id);

    for (i=0; i < genotype->ngenes; i++) 
    {
        genotype->mRNAdecay[i] = exp(0.4909*gasdev(&master_seed)-3.20304);
        while (genotype->mRNAdecay[i]<0.0)
          genotype->mRNAdecay[i] = exp(0.4909*gasdev(&master_seed)-3.20304);
        genotype->proteindecay[i]=-1.0;
        while (genotype->proteindecay[i] < 0.0) 
        {
            if (ran1(&master_seed) < 0.08421)
              genotype->proteindecay[i] = (float)EPSILON; /* was 0.0. changed to non-zero for stability*/
            else genotype->proteindecay[i] = exp(0.7874*gasdev(&master_seed)-3.7665);
        }
        /* dilution no longer done here, because it is now variable (function of instantaneous growth rate) */
        genotype->translation[i] = exp(0.7406*gasdev(&master_seed)+4.56);
        while (genotype->translation[i] < 0.0)
          genotype->translation[i] = exp(0.7406*gasdev(&master_seed)+4.56);

        /* make the activations the same in each copy */

        if(i<genotype->ntfgenes)
        {
            if (ran1(&master_seed)<PROB_ACTIVATING) 
            {
                genotype->N_act++;  
                for (p=0; p < MAX_COPIES; p++) 
                    genotype->activating[i][p] = 1;
            } 
            else 
            {
                genotype->N_rep++;
                for (p=0; p < MAX_COPIES; p++) 
                    genotype->activating[i][p] = 0;
            }
        }

    //    for (p=0; p < MAX_COPIES; p++) 
    //      LOG_NOFUNC("%d ", genotype->activating[i][p]);

        j = trunc(NUM_K_DISASSEMBLY * ran1(&master_seed));

        for (p=0; p < MAX_COPIES; p++) 
          genotype->pic_disassembly[i][p] = kdis[j];
    }
//    LOG_NOFUNC("\n");
 
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, but the same random TF sequence,
 * hindrance positions, replication times, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype,                        
                         float kdis[])
{ 
    int i,j,k,p;

    genotype->ngenes=10;
    genotype->ntfgenes=8;
    genotype->nproteins=10;    
    
    /* initially, each protein has only one copy of gene and mRNA*/
    for(i=0;i<genotype->nproteins;i++)
    {
        genotype->protein_pool[i][0][0]=1;
        genotype->protein_pool[i][1][0]=i;
        genotype->which_protein[i]=i;
    }
    
    initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*MAX_COPIES*NGENES, MAX_COPIES, genotype->ngenes);

    initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*MAX_COPIES*TFGENES, MAX_COPIES, genotype->ntfgenes);
 
  /* We now generate the complementary sequence of BS that are on the non-template strand.
   * The complementary sequence is used to search for BS that on the non-template strand.  
   * We also assume that all the TFs have strong orientation preference.*/  
    for(i=0;i< genotype->ntfgenes;i++)
    {
        for(j=0;j<MAX_COPIES;j++)
        {
            for(k=0;k<TF_ELEMENT_LEN;k++)
            {
                switch (genotype->tf_seq[i][j][k])
                {
                    case 'a': genotype->tf_seq_rc[i][j][k]='t'; break;
                    case 't': genotype->tf_seq_rc[i][j][k]='a'; break;
                    case 'c': genotype->tf_seq_rc[i][j][k]='g'; break;
                    case 'g': genotype->tf_seq_rc[i][j][k]='c'; break;
                }
            }
        }
    }
       
    initialize_genotype_fixed(genotype, kdis);

    /* start number of copies of gene at current_ploidy */
//    for (p=0; p < genotype->ngenes; p++) 
//    {
//        genotype->copies[p] = current_ploidy;
//    }
    
    calc_all_binding_sites(genotype);
}


/* calc the numbers of TF binding configurations:
 * determine the maximal number of repressors that can bind to a promoter
 * and the numbers of activators that can bind given the number of repressors that bind
 * Use these information to accelerate the calculation of rep-to-act ratios.
 */
//void calc_configurations(Genotype *genotype, int gene_id)
//{
//    int j,k,l;
//    int hindered_rep;
//    int hindered_act; 
//    int N_rep_BS=0;
//    int N_act_BS=0;
//    
//    genotype->N_configurations[gene_id][0]=1;
//    
//    genotype->max_N_rep_bound[gene_id]=0;   
//    
//    genotype->max_N_act_bound[gene_id]=0;
//
//    for(j=0;j<genotype->binding_sites_num[gene_id];j++)
//    {
//        hindered_rep=0;
//        hindered_act=0;
//
//        for(k=0;k<j;k++)
//        {
//            /* for BS within hindrance range*/
//           if((genotype->all_binding_sites[gene_id][k].BS_pos) > (genotype->all_binding_sites[gene_id][j].BS_pos-TF_ELEMENT_LEN-2*HIND_LENGTH))
//           {
//               if(genotype->activating[genotype->all_binding_sites[gene_id][k].tf_id][0]==0) hindered_rep++;
//               else hindered_act++;
//           }
//        }
//
//        if(genotype->activating[genotype->all_binding_sites[gene_id][j].tf_id][0]==0)
//        {
//            N_rep_BS++;
//            
//            if((N_rep_BS-hindered_rep)!=1)
//            {
//                genotype->max_N_rep_bound[gene_id]++; 
//                
//                for(l=1;l<genotype->max_N_rep_bound[gene_id]+1;l++)
//                {
//                    genotype->N_configurations[gene_id][l]=genotype->N_configurations[gene_id][l-1];
//                }
//            }
//        }
//        else
//        {
//            N_act_BS++;
//            
//            if((N_act_BS-hindered_act)!=1)
//            {
//                for(l=0;l<genotype->max_N_rep_bound[gene_id]-hindered_rep+1;l++)
//                    genotype->N_configurations[gene_id][l]++;
//            }
//        }
//    } 
//    
//    for(j=0;j<genotype->max_N_rep_bound[gene_id];j++)
//    {
//        genotype->max_N_act_bound[gene_id]=(genotype->max_N_act_bound[gene_id]>genotype->N_configurations[gene_id][j])?
//                                            genotype->max_N_act_bound[gene_id]:genotype->N_configurations[gene_id][j];
//    }
//    genotype->max_N_act_bound[gene_id]--;
//}

/*
 * compute the list binding sites for specified gene and gene copy
 */
/*considering computing all the possible configurations here*/
void calc_all_binding_sites_copy(Genotype *genotype, int gene_id, int *max_binding_sites_alloc)
{
    int i, j, k, match,match_rc;
    int N_hindered_BS=0;
    int N_hindered_BS_copy;
    int FOUND_AT_THE_SAME_POS;
    int N_binding_sites=0;  


    genotype->N_act_BS[gene_id]=0;
    genotype->N_rep_BS[gene_id]=0;
    genotype->max_hindered_sites[gene_id]=0;
  
    //make the code looks cleaner
    #define p_all_BS genotype->all_binding_sites[gene_id]
    #define WHICH_TF genotype->protein_pool[k][1][0]

    //some helper pointer 
    char *tf_seq;
    char *cis_seq;
    char *tf_seq_rc; 
    cis_seq=&(genotype->cisreg_seq[gene_id][0][0]);  
  
   for (i=0; i < CISREG_LEN-TF_ELEMENT_LEN; i++) /* scan promoter */
   {  
        /*calc the number of BS within the hindrance range*/
        N_hindered_BS=0;
        
        if(N_binding_sites>0)
        {
            for(j=0;j<N_binding_sites;j++)
            {
               if(p_all_BS[j].BS_pos> i-TF_ELEMENT_LEN-2*HIND_LENGTH)
                    N_hindered_BS++;
            }
        }  
        
        /* set flags and copy*/
        N_hindered_BS_copy=N_hindered_BS; 
        
        FOUND_AT_THE_SAME_POS=0;        
        
        /* loop through TFs */
        /* NOTE: because of duplication, ntfgenes can be bigger than the number of actual tfs*/
        for (k=0; k < genotype->nproteins-2; k++) 
        {    
            tf_seq=&(genotype->tf_seq[WHICH_TF][0][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[WHICH_TF][0][0]); 
            
            /*find BS on the template strand*/
            match=0;

            for (j=i; j < i+TF_ELEMENT_LEN; j++) 
            {
                if (cis_seq[j] == tf_seq[j-i]) match++;
            }

            if (match >= NMIN)
            {  
                if (N_binding_sites + 1 >= MAXELEMENTS) 
                {                    
                    p_all_BS = realloc(p_all_BS, 2*MAXELEMENTS*sizeof(AllTFBindingSites));
                    if (!p_all_BS) 
                    {
//                      LOG_ERROR_NOCELLID("realloc of all_binding_sites to binding_sites_num=%d \n",*max_binding_sites_alloc);
                      exit(1);
                    }
//                    else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d succeeded\n", *max_binding_sites_alloc);
                }
                
                p_all_BS[N_binding_sites].tf_id = k;
                genotype->all_binding_sites[gene_id][N_binding_sites].Koff = Koff[TF_ELEMENT_LEN-match];                
                p_all_BS[N_binding_sites].BS_pos = i ;           
                
                if(!FOUND_AT_THE_SAME_POS) /* first BS found at this pos */                   
                {
                    p_all_BS[N_binding_sites].N_hindered = N_hindered_BS;
                    
                    N_hindered_BS_copy++;
                    
                    FOUND_AT_THE_SAME_POS=1;
                }
                else /* found BS of another TF at this pos again*/
                {
                    p_all_BS[N_binding_sites].N_hindered = N_hindered_BS_copy;
                    
                    N_hindered_BS_copy++;
                }
                
                N_binding_sites++;

                if(genotype->activating[WHICH_TF][0]==1) genotype->N_act_BS[gene_id]++;
            }
            else /*find BS on the non-template strand*/
            {
                match_rc=0;

                for (j=i; j < i+TF_ELEMENT_LEN; j++) 
                {
                    if (cis_seq[j] == tf_seq_rc[j-i]) match_rc++;
                }

                if (match_rc >= NMIN)
                {
                    /**********************************************************************/     
                    if (N_binding_sites + 1 >= MAXELEMENTS) 
                    {                  
                      p_all_BS = realloc(p_all_BS, 2*MAXELEMENTS*sizeof(AllTFBindingSites));

                      if (!p_all_BS) 
                        {
    //                      LOG_ERROR_NOCELLID("realloc of all_binding_sites to N_binding_sites = %d failed.\n", *max_binding_sites_alloc);
                          exit(1);
                        }
    //                  else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to N_binding_sites = %d succeeded\n", *max_binding_sites_alloc);
                    }
                    /************************************************************************************************************/
                    p_all_BS[N_binding_sites].tf_id = k;
                    p_all_BS[N_binding_sites].Koff = Koff[TF_ELEMENT_LEN-match_rc];                   
                    p_all_BS[N_binding_sites].BS_pos = i;                    
                   
                    if(!FOUND_AT_THE_SAME_POS) /* first BS found at this pos */                   
                    {
                        p_all_BS[N_binding_sites].N_hindered = N_hindered_BS;

                        N_hindered_BS_copy++;

                        FOUND_AT_THE_SAME_POS=1;
                    }
                    else /* found BS of another TF at this pos again*/
                    {
                        p_all_BS[N_binding_sites].N_hindered = N_hindered_BS_copy;

                        N_hindered_BS_copy++;
                    }
                    
                    N_binding_sites++;
                    
                    if(genotype->activating[WHICH_TF][0]==1) genotype->N_act_BS[gene_id]++;
                }
            } 
        }/* looping through TFs ends */
    }/*end of promoter scanning*/ 
    genotype->binding_sites_num[gene_id]=N_binding_sites;  
    genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);

    /* calc max_hindered_sites */
    for(i=0;i<N_binding_sites;i++)
    {
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id] > p_all_BS[i].N_hindered)?
                                      genotype->max_hindered_sites[gene_id] : p_all_BS[i].N_hindered;
    }
}

/*
 * compute the list of binding sites for the specified number of gene
 * copies
 */
void calc_all_binding_sites(Genotype *genotype)
{
    int max_binding_site_alloc= MAXELEMENTS;
    int gene_id;  

    for(gene_id=0;gene_id < genotype->ngenes;gene_id++)
    {        
        if(genotype->re_calc[gene_id][2]) /* do not calculate the binding sites if there's no mutation in the promoter*/
        {
            calc_all_binding_sites_copy(genotype,gene_id,&max_binding_site_alloc);
        
//            calc_configurations(genotype, gene_id);
            
            genotype->re_calc[gene_id][2]=0;
        }
    }
}

int mod(int a, int b) // b is the base
{
    if(a>=0)
        return (a%b);
    else
        return (abs(a+b)%b);
}

/* calc with approximation, but greatly enhance performance*/

float calc_ratio_act_to_rep_approximation(AllTFBindingSites *BS_info,
                           int ntfgenes,
                           int max_N_hindered_BS,
                           int N_BS,
                           int N_act_BS,
                           int N_rep_BS, 
                           int activating[NGENES][MAX_COPIES],
                           int approximation,
                           float protein_conc[NGENES])
{
    double ratio_matrices[max_N_hindered_BS+1][approximation][approximation];   
    double transition_matrix[approximation][approximation];
    double sum,prob_act_over_rep=0.0;    
    double product_of_freq;    
    float Kon[ntfgenes];      
    
    int pos_of_last_record;    
    int pos_next_record;
    int i,j,k,m,n;
//    int n_col1=1;
//    int n_col2=1;
//    int n_row=1;
    
    /*calc Kon based on TF concentration*/
    for(i=0;i<ntfgenes;i++)
    {
        Kon[i]=kon*protein_conc[i];
    }
    
    /* initializing matrices to all zeros */
    for(i=0;i<max_N_hindered_BS+1;i++)
    {
        for(j=0;j<approximation;j++)
        {
            for(k=0;k<approximation;k++)
            {
                ratio_matrices[i][j][k]=0.0;
            }
        }
    }
    
    for(j=0;j<approximation;j++)
    {
        for(k=0;k<approximation;k++)
        {
            transition_matrix[j][k]=0.0;            
        }
    }
    
    /* body of the forward algorithm*/    
    pos_next_record=0; //where in the ratio_matrices to put the next record
    
    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
    
    if(activating[BS_info[0].tf_id][0]==1) // if a activator binds to this BS
    {
        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col2++;
    }
    else
    {
        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
    }    
    
    for(m=1;m<N_BS;m++)
    {
        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);

        product_of_freq = Kon[BS_info[m].tf_id];

        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
        {
            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
            {
                product_of_freq*=BS_info[n].Koff;
            }
        }

        switch(activating[BS_info[m].tf_id][0])
        {
            case 1: // a BS of activators
//                n_col2++;
                
                if(m-BS_info[m].N_hindered!=0)
                {
                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              

                  for(i=0;i<approximation;i++)
                  {
                      transition_matrix[i][0]=0.0;
                      
//                      n_col1=(n_col2<N_act_bound[i])? n_col2:N_act_bound[i];

                      for(j=1;j<approximation;j++)
                      {
                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
                      }
                  }
                }
                else
                {
                    for(i=0;i<approximation;i++)
                    {
//                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
                        
                        for(j=0;j<approximation;j++)
                        {
                            transition_matrix[i][j]=0.0;
                        }
                    }

                    transition_matrix[0][1]=1.0;
                }
                    
                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              

                for(i=0;i<approximation;i++)
                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
                    
                    for(j=0;j<approximation;j++)
                    {
                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
                                                      product_of_freq*transition_matrix[i][j];                            
                    }
                }
                break;

            case 0: // a BS of repressors
                
//                n_row++;
//                
//                n_row=(n_row < max_N_rep_bound+1)?n_row:max_N_rep_bound+1;
                
                if(m-BS_info[m].N_hindered!=0)
                {
                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1);                                 

//                  n_col1=(n_col2<N_act_bound[0])?n_col2:N_act_bound[0];
                  
                  for(j=0;j<approximation;j++)
                  {
                      transition_matrix[0][j]=0.0;
                  }
                  
                  for(i=1;i<approximation;i++)
                  {
//                      n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
                      
                      for(j=0;j<approximation;j++)
                      {
                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
                      }
                  }
                }
                else
                {
                    for(i=0;i<approximation;i++)
                    {
//                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
                        
                        for(j=0;j<approximation;j++)
                        {
                            transition_matrix[i][j]=0.0;
                        }
                    }

                    transition_matrix[1][0]=1.0;
                }
                    
                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);

                for(i=0;i<approximation;i++)
                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
                     
                    for(j=0;j<approximation;j++)
                    {
                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
                                                      product_of_freq*transition_matrix[i][j];                            
                    }
                }
                break;
        }
    }

    sum=0.0;

    for(i=0;i<approximation;i++)
    {
        for(j=0;j<approximation;j++)
        {
            sum+=ratio_matrices[pos_next_record][i][j];
        }
    }   

    for(i=0;i<approximation;i++)
    {
        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
        
        for(;j<approximation;j++)
        {
            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
        }
    } 

    return (float)(prob_act_over_rep/sum);
     // end of the forward algorithm   
}


/*this use max_N_rep_bound and *N_act_bound */
//float calc_ratio_act_to_rep(AllTFBindingSites *BS_info,
//                           int ntfgenes,
//                           int max_N_hindered_BS,
//                           int N_BS,
//                           int N_act_BS,
//                           int N_rep_BS, 
//                           int activating[NGENES][MAX_COPIES],
//                           int max_N_rep_bound,
//                           int max_N_act_bound,
//                           int *N_act_bound,
//                           float protein_conc[NGENES])
//{
//    double ratio_matrices[max_N_hindered_BS+1][max_N_rep_bound+1][max_N_act_bound+1];   
//    double transition_matrix[max_N_rep_bound+1][max_N_act_bound+1];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;    
//    float Kon[ntfgenes];      
//    
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,k,m,n;
//    int n_col1=1;
//    int n_col2=1;
//    int n_row=1;
//    
//    /*calc Kon based on TF concentration*/
//    for(i=0;i<ntfgenes;i++)
//    {
//        Kon[i]=kon*protein_conc[i];
//    }
//    
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_N_hindered_BS+1;i++)
//    {
//        for(j=0;j<max_N_rep_bound+1;j++)
//        {
//            for(k=0;k<N_act_bound[j];k++)
//            {
//                ratio_matrices[i][j][k]=0.0;
//            }
//        }
//    }
//    
//    for(j=0;j<max_N_rep_bound+1;j++)
//    {
//        for(k=0;k<N_act_bound[j];k++)
//        {
//            transition_matrix[j][k]=0.0;            
//        }
//    }
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    
//    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
//    
//    if(activating[BS_info[0].tf_id][0]==1) // if a activator binds to this BS
//    {
//        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col2++;
//    }
//    else
//    {
//        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
//    }    
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
//
//        product_of_freq = Kon[BS_info[m].tf_id];
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//            {
//                product_of_freq*=BS_info[n].Koff;
//            }
//        }
//
//        switch(activating[BS_info[m].tf_id][0])
//        {
//            case 1: // a BS of activators
//                n_col2++;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
//                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//                  for(i=0;i<n_row;i++)
//                  {
//                      transition_matrix[i][0]=0.0;
//                      
//                      n_col1=(n_col2<N_act_bound[i])? n_col2:N_act_bound[i];
//
//                      for(j=1;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
//                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[0][1]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              
//
//                for(i=0;i<n_row;i++)
//                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                    
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//
//            case 0: // a BS of repressors
//                
//                n_row++;
//                
//                n_row=(n_row < max_N_rep_bound+1)?n_row:max_N_rep_bound+1;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1);                                 
//
//                  n_col1=(n_col2<N_act_bound[0])?n_col2:N_act_bound[0];
//                  
//                  for(j=0;j<n_col1;j++)
//                  {
//                      transition_matrix[0][j]=0.0;
//                  }
//                  
//                  for(i=1;i<n_row;i++)
//                  {
//                      n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                      
//                      for(j=0;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
//                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[1][0]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
//
//                for(i=0;i<n_row;i++)
//                {
//                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                     
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//        }
//    }
//
//    sum=0.0;
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        for(j=0;j<N_act_bound[i];j++)
//        {
//            sum+=ratio_matrices[pos_next_record][i][j];
//        }
//    }   
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
//        
//        for(;j<N_act_bound[i];j++)
//        {
//            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
//        }
//    } 
//
//    return (float)(prob_act_over_rep/sum);
//     // end of the forward algorithm   
//}


/* this one does not use N_act_bound*/
//float calc_ratio_act_to_rep(AllTFBindingSites *BS_info,
//                           int ntfgenes,
//                           int max_N_hindered_BS,
//                           int N_BS,
//                           int N_act_BS,
//                           int N_rep_BS, 
//                           int activating[NGENES][MAX_COPIES],
//                           int max_N_rep_bound,
//                           int max_N_act_bound, 
//                           int *N_act_bound,
//                           float protein_conc[NGENES])
//{
//    double ratio_matrices[max_N_hindered_BS+1][max_N_rep_bound+1][max_N_act_bound+1];   
//    double transition_matrix[max_N_rep_bound+1][max_N_act_bound+1];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;    
//    float Kon[ntfgenes];      
//    
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,k,m,n;
//    int n_col1=1;
//    int n_col2=1;
//    int n_row=1;
//    
//    /*calc Kon based on TF concentration*/
//    for(i=0;i<ntfgenes;i++)
//    {
//        Kon[i]=kon*protein_conc[i];
//    }
//    
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_N_hindered_BS+1;i++)
//    {
//        for(j=0;j<max_N_rep_bound+1;j++)
//        {
//            for(k=0;k<max_N_act_bound+1;k++)
//            {
//                ratio_matrices[i][j][k]=0.0;
//            }
//        }
//    }
//    
//    for(j=0;j<max_N_rep_bound+1;j++)
//    {
//        for(k=0;k<max_N_act_bound+1;k++)
//        {
//            transition_matrix[j][k]=0.0;            
//        }
//    }
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    
//    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
//    
//    if(activating[BS_info[0].tf_id][0]==1) // if a activator binds to this BS
//    {
//        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col1++;
//    }
//    else
//    {
//        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
//    }    
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
//
//        product_of_freq = Kon[BS_info[m].tf_id];
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//            {
//                product_of_freq*=BS_info[n].Koff;
//            }
//        }
//
//        switch(activating[BS_info[m].tf_id][0])
//        {
//            case 1: // a BS of activators
//                n_col1++;
//                
//                n_col1=(n_col1<max_N_act_bound+1)? n_col1:max_N_act_bound+1;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
//                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//                  for(i=0;i<n_row;i++)
//                  {
//                      transition_matrix[i][0]=0.0;
//                      
////                      n_col1=(n_col2<N_act_bound[i])? n_col2:N_act_bound[i];
//
//                      for(j=1;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
////                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[0][1]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              
//
//                for(i=0;i<n_row;i++)
//                {
////                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                    
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//
//            case 0: // a BS of repressors
//                
//                n_row++;
//                
//                n_row=(n_row < max_N_rep_bound+1)?n_row:max_N_rep_bound+1;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1);                                 
//
////                  n_col1=(n_col2<N_act_bound[0])?n_col2:N_act_bound[0];
//                  
//                  for(j=0;j<n_col1;j++)
//                  {
//                      transition_matrix[0][j]=0.0;
//                  }
//                  
//                  for(i=1;i<n_row;i++)
//                  {
////                      n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                      
//                      for(j=0;j<n_col1;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
////                        n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                        
//                        for(j=0;j<n_col1;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[1][0]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
//
//                for(i=0;i<n_row;i++)
//                {
////                    n_col1=(n_col2<N_act_bound[i])?n_col2:N_act_bound[i];
//                     
//                    for(j=0;j<n_col1;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//        }
//    }
//
//    sum=0.0;
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        for(j=0;j<max_N_act_bound+1;j++)
//        {
//            sum+=ratio_matrices[pos_next_record][i][j];
//        }
//    }   
//
//    for(i=0;i<max_N_rep_bound+1;i++)
//    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
//        
//        for(;j<max_N_act_bound+1;j++)
//        {
//            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
//        }
//    } 
//
//    return (float)(prob_act_over_rep/sum);
//     // end of the forward algorithm   
//}


/* this one use neither max_N_rep_bound or *N_act_bound*/
//float calc_ratio_act_to_rep(AllTFBindingSites *BS_info,
//                           int ntfgenes,
//                           int max_N_hindered_BS,
//                           int N_BS,
//                           int N_act_BS,
//                           int N_rep_BS, 
//                           int activating[NGENES][MAX_COPIES],
//                           int max_N_rep_bound,
//                           int *N_act_bound,
//                           float protein_conc[NGENES])
//{
//    double ratio_matrices[max_N_hindered_BS+1][N_rep_BS+1][N_act_BS+1];   
//    double transition_matrix[N_rep_BS+1][N_act_BS+1];
//    double sum,prob_act_over_rep=0.0;    
//    double product_of_freq;    
//    float Kon[ntfgenes];      
//    
//    int pos_of_last_record;    
//    int pos_next_record;
//    int i,j,k,m,n;
//    int n_col=1;    
//    int n_row=1;
//    
//    /*calc Kon based on TF concentration*/
//    for(i=0;i<ntfgenes;i++)
//    {
//        Kon[i]=kon*protein_conc[i];
//    }
//    
//    /* initializing matrices to all zeros */
//    for(i=0;i<max_N_hindered_BS+1;i++)
//    {
//        for(j=0;j<N_rep_BS+1;j++)
//        {
//            for(k=0;k<N_act_BS+1;k++)
//            {
//                ratio_matrices[i][j][k]=0.0;
//            }
//        }
//    }
//    
//    for(j=0;j<max_N_rep_bound+1;j++)
//    {
//        for(k=0;k<N_act_BS+1;k++)
//        {
//            transition_matrix[j][k]=0.0;            
//        }
//    }
//    
//    /* body of the forward algorithm*/    
//    pos_next_record=0; //where in the ratio_matrices to put the next record
//    
//    ratio_matrices[pos_next_record][0][0]=BS_info[0].Koff;   
//    
//    if(activating[BS_info[0].tf_id][0]==1) // if a activator binds to this BS
//    {
//        ratio_matrices[pos_next_record][0][1]=Kon[BS_info[0].tf_id]; 
//        n_col++;
//    }
//    else
//    {
//        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];  
//        n_row++;
//    }    
//    
//    for(m=1;m<N_BS;m++)
//    {
//        pos_next_record=mod(pos_next_record+1,max_N_hindered_BS+1);
//
//        product_of_freq = Kon[BS_info[m].tf_id];
//
//        if(BS_info[m].N_hindered) // if binding to the current BS hinders other BS
//        {
//            for(n=m-BS_info[m].N_hindered;n<=m-1;n++)
//            {
//                product_of_freq*=BS_info[n].Koff;
//            }
//        }
//
//        switch(activating[BS_info[m].tf_id][0])
//        {
//            case 1: // a BS of activators
//                n_col++;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
//                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              
//
//                  for(i=0;i<n_row;i++)
//                  {
//                      transition_matrix[i][0]=0.0;
//
//                      for(j=1;j<n_col;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {
//                        for(j=0;j<n_col;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[0][1]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              
//
//                for(i=0;i<n_row;i++)
//                {
//                    for(j=0;j<n_col;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//
//            case 0: // a BS of repressors
//                
//                n_row++;
//                
//                if(m-BS_info[m].N_hindered!=0)
//                {
//                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); 
//                  
//                  for(j=0;j<n_col;j++)
//                  {
//                      transition_matrix[0][j]=0.0;
//                  }
//                  
//                  for(i=1;i<n_row;i++)
//                  {                      
//                      for(j=0;j<n_col;j++)
//                      {
//                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
//                      }
//                  }
//                }
//                else
//                {
//                    for(i=0;i<n_row;i++)
//                    {                        
//                        for(j=0;j<n_col;j++)
//                        {
//                            transition_matrix[i][j]=0.0;
//                        }
//                    }
//
//                    transition_matrix[1][0]=1.0;
//                }
//                    
//                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);
//
//                for(i=0;i<n_row;i++)
//                {                     
//                    for(j=0;j<n_col;j++)
//                    {
//                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
//                                                      product_of_freq*transition_matrix[i][j];                            
//                    }
//                }
//                break;
//        }
//    }
//
//    sum=0.0;
//
//    for(i=0;i<N_rep_BS+1;i++)
//    {
//        for(j=0;j<n_col;j++)
//        {
//            sum+=ratio_matrices[pos_next_record][i][j];
//        }
//    }   
//
//    for(i=0;i<N_rep_BS+1;i++)
//    {
//        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe    
//        
//        for(;j<n_col;j++)
//        {
//            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
//        }
//    } 
//
//    return (float)(prob_act_over_rep/sum);
//     // end of the forward algorithm   
//}

int add_fixed_event(int i,
                    int p,
                    float t,
                    FixedEvent **start,
                    FixedEvent **last)
{
  FixedEvent *newtime;
  int pos;

  newtime = (FixedEvent *)malloc(sizeof(FixedEvent));
  if (!newtime) {
//    printf("Out of memory\n");
    exit(1);
  }
  newtime->gene_id = i;
  newtime->copy = p;
  newtime->time = t;
//  LOG_VERBOSE_NOCELLID("adding event at time=%f for gene=%d (copy=%d)\n", t, i, p);
  pos = sls_store(newtime, start, last);
  return pos;
}

void add_time_point(float time,
                    float conc,
                    TimeCourse **start,
                    TimeCourse **last)
{
  TimeCourse *newtime;
  
  newtime = (TimeCourse *)malloc(sizeof(TimeCourse));
  if (!newtime) {
//    printf("Out of memory\n");
    exit(1);
  }
  newtime->time = time;
  newtime->concentration = conc;
  sls_store_end2(newtime, start, last);
}

void add_fixed_event_end(int gene_id,
                         int gene_copy,
                         float t,
                         FixedEvent **start,
                         FixedEvent **last)
{
  FixedEvent *newtime;
  
  newtime = (FixedEvent *)malloc(sizeof(FixedEvent));
  if (!newtime) {
//    printf("Out of memory\n");
    exit(1);
  }
  newtime->gene_id = gene_id;
  newtime->copy = gene_copy;
  newtime->time = t;
//  LOG_VERBOSE_NOCELLID("adding event end at time=%f for gene=%d (copy=%d)\n", t, gene_id, gene_copy);
  sls_store_end(newtime, start, last);
}

void delete_fixed_event(int gene_id,
                        int gene_copy,
                        int i,
                        FixedEvent **start,
                        FixedEvent **last)
{
  FixedEvent *info, *lastinfo;
  int j, done;
  
  j = -1;
  done = 0;
  info = *start;
  while (info) {
    if ((info->gene_id==gene_id && info->copy==gene_copy)) {
      j++;
      if (j == i) {
        if (info == *start) {
          *start = info->next;
          if (info == *last) *last = NULL;
        } else {
          lastinfo->next = info->next;
          if (info == *last) *last = lastinfo;
        }
        done = 1;
      } else {
        lastinfo = info;
        info = info->next;
      }
    } else {
      lastinfo = info;
      info = info->next;
    }
  }
  if (done == 0) {
//    LOG_ERROR_NOCELLID("In %d elements, couldn't find element %d to delete in gene %d (copy %d)\n",
//                       j+1, i, gene_id, gene_copy);
  }
  free(info);
}

void delete_fixed_event_start(FixedEvent **start,
                              FixedEvent **last)
{
  FixedEvent *info;
  
  info = *start;
  *start = info->next;
  if (*last == info) *last = NULL;
  free(info);
}

/*
 * initialize the cell state with the specified initial protein
 * concentration, mean mRNA number and mRNA decay and whether to do
 * burn-in of high kon rate or not
 */
void initialize_cell(CellState *state,
                     int ngenes,
                     int nproteins,
                     int *protein_pool[NPROTEINS][2], 
                     float mRNAdecay[NGENES],
                     float meanmRNA[NGENES],
                     float init_protein_conc[NPROTEINS],
                     long int *seed)
{
    int i, j, k, totalmRNA;
    float t;

    /* start cell size at 1.0 */
    state->cell_size = 1.0; 
    state->cell_size_copy = 1.0;

    /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
    state->growth_rate = 0.0;
    state->cumulative_growth_rate =1.0;
    
    state->env_change=0;

    state->mRNA_transcr_time_end = NULL;
    state->mRNA_transcr_time_end_last = NULL;
    state->mRNA_transl_time_end = NULL;
    state->mRNA_transl_time_end_last = NULL;
    state->env0_time_end = NULL;
    state->env0_time_end_last = NULL;
    state->env1_time_end = NULL;
    state->env1_time_end_last = NULL;

    for (i=0; i < ngenes; i++) 
    {
        for (j=0; j < MAX_COPIES; j++) 
        {
            state->active[i][j] = NUC_NO_PIC;// was ON_WITH_NUCLEOSOME
        }

        totalmRNA = (int) poidev(meanmRNA[i],seed);
        state->mRNA_nuclear_num[i] = (int) bnldev(STARTNUCLEUS, totalmRNA, seed);
        state->mRNA_cyto_num[i] = totalmRNA - state->mRNA_nuclear_num[i];
        state->mRNA_transl_cyto_num[i] = 0;

        for (k=0; k<state->mRNA_cyto_num[i]; k++) 
        {
            t = expdev(seed) / mRNAdecay[i];
            if (t < TTRANSLATION) 
            {
                (state->mRNA_cyto_num[i])--;
                (state->mRNA_transl_cyto_num[i])++;
        //        LOG_VERBOSE("add translation event time=%g for gene=%d\n", (TTRANSLATION-t), i);
                add_fixed_event(i, -1, TTRANSLATION-t, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
            }
        } 

        int total_mRNA_transcribing = (int) poidev(meanmRNA[i]*TTRANSCRIPTION*mRNAdecay[i], seed);

        /* split it up evenly between the copies */
        int mRNA_copy1 = trunc(total_mRNA_transcribing/current_ploidy);
        int mRNA_copy2 = total_mRNA_transcribing - mRNA_copy1;

        for (j=0; j < MAX_COPIES; j++) 
        {
            if (j < current_ploidy)  
            {
                state->mRNA_transcr_num[i][j] = (j==0) ? mRNA_copy1 : mRNA_copy2;
        //        LOG_VERBOSE_NOCELLID("initializing state->mRNA_transcr_num[%2d][%2d]=%d\n", i, j, state->mRNA_transcr_num[i][j]);
                for (k=0; k < state->mRNA_transcr_num[i][j]; k++)
                  add_fixed_event(i, j, ran1(seed)*TTRANSCRIPTION, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));
            } 
            else 
            {
                state->mRNA_transcr_num[i][j] = 0;
            }
        }
    }
    
    /* initiate protein concentration*/
    for (i=0; i < ngenes; i++) 
            state->gene_specific_protein_conc[i] = init_protein_conc[i];
    
    for(i=0;i<nproteins;i++)
    {
        state->protein_conc[i]=0.0;
        
        for(j=0;j<protein_pool[i][0][0];j++)
            state->protein_conc[i]+=state->gene_specific_protein_conc[protein_pool[i][1][j]];
    }
    
//  for (j=0; j < MAX_COPIES; j++) {  
//    int pos = 0;
//    for (i=0; i < NGENES; i++) {
//      if (genotype->copies[i] > j) {  
//        state->state_change_ids[ACETYLATION_STATE][j][pos] = i;
//        LOG_VERBOSE("Initializing statechange gene=%d, ploidy=%d state_change_ids[%d][%d]=%d\n", i, j, j, 
//                    pos, state->state_change_ids[ACETYLATION_STATE][j][pos]);
//        pos++;
//      }
//    }
//  }
  // TODO: add more env changing 
    add_fixed_event(0,0,duration_env0,&(state->env0_time_end),&(state->env0_time_end_last));
}

/* 
 * change in the number of mRNAs in the cytoplasm affects the kon
 * rates, note we must have already updated protein_conc first
 */
void change_mRNA_cytoplasm(int i,
                           Genotype *genotype,
                           CellState *state,
                           GillespieRates *rates)
{
  float salphc; 
  
  /* number of mRNAs in cytoplasm affects kon rates */
  salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / state->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  
//  LOG_VERBOSE("change_mRNA_cytoplasm[%d]: mRNA=%d, transl rate=%g, protein decay=%g, salphc=%g\n", 
//              i, state->mRNA_cyto_num[i], genotype->translation[i], state->konvalues[i][KON_PROTEIN_DECAY_INDEX], salphc); 

//  state->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / state->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  state->konvalues[i][KON_SALPHC_INDEX] = salphc;
}

void calc_all_rates(Genotype *genotype,
                    CellState *state,
                    GillespieRates *rates,                     
                    int UPDATE_KONVALUES)
{
    int i;
    float salphc; 
    float protein_decay;

    if(UPDATE_KONVALUES) /*Since these values are also updated in update_protein_conc_cell_size, this is only used in initilization*/
    {
        for (i=0; i < genotype->ngenes; i++) 
        {
            /* if protein decay is otherwise going to fall below cut-off, use aging term */
            protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;

            salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);

            state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = protein_decay;

            state->konvalues[i][KON_SALPHC_INDEX] = salphc;

//          LOG_VERBOSE("protein decay[%d]=%g\n", i, state->konvalues[i][KON_PROTEIN_DECAY_INDEX]);
        }
    }
    
    /* reset rates and operations */
    rates->transport=0.0;
    rates->mRNAdecay=0.0;
    rates->pic_disassembly=0.0;
    rates->acetylation=0.0;
    rates->deacetylation=0.0;
    rates->transcript_init=0;
    rates->pic_assembly=0.0;
    rates->subtotal=0.0;
    
    for(i=0;i<genotype->ngenes;i++)
    {
        rates->acetylation_rate[i]=0.0;
        rates->deacetylation_rate[i]=0.0;
        rates->pic_assembly_rate[i]=0.0;
        rates->pic_disassembly_rate[i]=0.0;
        rates->mRNAdecay_rate[i]=0.0;
        rates->transport_rate[i]=0.0; 
    } 

    /* update rates*/
    for(i=0;i<genotype->ngenes;i++)
    {
         /* update RNA transport rate*/
        rates->transport_rate[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
        rates->transport += rates->transport_rate[i];
//        rates->transport_operations++;
//        LOG_VERBOSE("mRNA_nuclear_num=%d, initializing transport[%d]=%g\n", state->mRNA_nuclear_num[i], i, rates->transport_rate[i]);
        
        /* update mRNAdecay rate based on the total number of mRNAs in both
           cytoplasm (mRNA_cyto_num) and ones that have only just recently arrived
           (mRNA_transl_cyto_num) */
        rates->mRNAdecay_rate[i] = genotype->mRNAdecay[i] * ((float) state->mRNA_cyto_num[i] + (float) state->mRNA_transl_cyto_num[i]);
        rates->mRNAdecay += rates->mRNAdecay_rate[i];
//        rates->mRNAdecay_operations++;          
    }
    
    for (i=0; i < genotype->ngenes; i++) 
    {  
        /* calc TF binding prob distribution*/
        /* if 1) there is a place to copy the ratio from,
         * and 2) we can copy from that place*/
        if((genotype->re_calc[i][0]!=-1)&&(genotype->re_calc[genotype->re_calc[i][0]][1])) 
            state->Pact[i]=state->Pact[genotype->re_calc[i][0]];
            
        else /* otherwise, we need to calc the ratio*/
//            state->Pact[i]=calc_ratio_act_to_rep(genotype->all_binding_sites[i],
//                                                genotype->ntfgenes,
//                                                genotype->max_hindered_sites[i],
//                                                genotype->binding_sites_num[i],
//                                                genotype->N_act_BS[i],
//                                                genotype->N_rep_BS[i],
//                                                genotype->activating, 
//                                                genotype->max_N_rep_bound[i],
//                                                genotype->max_N_act_bound[i],
//                                                genotype->N_configurations[i],
//                                                state->protein_conc);
             state->Pact[i]=calc_ratio_act_to_rep_approximation(genotype->all_binding_sites[i],
                                                genotype->ntfgenes,
                                                genotype->max_hindered_sites[i],
                                                genotype->binding_sites_num[i],
                                                genotype->N_act_BS[i],
                                                genotype->N_rep_BS[i],
                                                genotype->activating, 
                                                4,
                                                state->protein_conc);
        
        /* calc other rates*/
        switch (state->active[i][0])
        {
            case NUC_NO_PIC:
                rates->acetylation_rate[i]=state->Pact[i]*ACETYLATE*Pact_scaling;
                rates->acetylation+=rates->acetylation_rate[i];
                break;
                
            case NO_NUC_NO_PIC:
                rates->deacetylation_rate[i]=(1-state->Pact[i])*DEACETYLATE;
                rates->deacetylation+=rates->deacetylation_rate[i];
                rates->pic_assembly_rate[i]=state->Pact[i]*PICASSEMBLY;
                rates->pic_assembly+=rates->pic_assembly_rate[i];                
                break;
                
            case PIC_NO_NUC: // Note: pic_disassembly_rate is a gene-specific constant, so it's defined in genotype
                rates->pic_disassembly_rate[i]=genotype->pic_disassembly[i][0];
                rates->pic_disassembly+=rates->pic_disassembly_rate[i]; 
                rates->transcript_init_rate[i]= 1;
                rates->transcript_init+=1;
                break;
        }
    }
    rates->subtotal+=rates->deacetylation;
    rates->subtotal+=rates->pic_assembly;
    rates->subtotal+=rates->acetylation;
    rates->subtotal+=rates->transport;
    rates->subtotal+=rates->mRNAdecay;
    rates->subtotal+=rates->pic_disassembly;
    rates->subtotal+=(float)rates->transcript_init*TRANSCRIPTINIT;
}

/* 
 * check to see if a fixed event ends within the time t
 *
 * returns:
 *  0 if there is no fixed event occuring before time t
 *  1 if a transcription event happens before time t
 *  2 if a translation event happens before time t
 *  3 if a gene replication event happens before time t
 */
 
int does_fixed_event_end(FixedEvent *mRNA_transl_time_end,
                         FixedEvent *mRNA_transcr_time_end,
                         FixedEvent *env0_time_end,
                         FixedEvent *env1_time_end,
			 float t) 
{
    int retval;
    float t1;
    float t2;
    float t3;
    float t4;

    if(mRNA_transcr_time_end == NULL && mRNA_transl_time_end==NULL && env0_time_end == NULL && env1_time_end == NULL)
    {
            retval =0;
    }
    else
    {
        t1 = mRNA_transcr_time_end ? mRNA_transcr_time_end->time : TIME_INFINITY;
        t2 = mRNA_transl_time_end ? mRNA_transl_time_end->time : TIME_INFINITY;
        t3 = env0_time_end ? env0_time_end->time : TIME_INFINITY;
        t4 = env1_time_end ? env1_time_end->time : TIME_INFINITY;

        if((t1 <= t2) && (t1 <= t) && (t1 <= t3) && (t1 <= t4)){
//			if (mRNA_transcr_time_end == NULL) retval = 0;
//      		else 
             retval = 1;	
        }
        else
        {
            if ((t2 <= t1) && (t2 <= t) && (t2 <= t3) && (t2 <= t4)) 
            {
//		        if (mRNA_transl_time_end == NULL)  retval = 0;
//		        else 
                retval = 2;
            }
            else
            {
                if ((t3 <= t1) && (t3 <= t) && (t3 <= t2) && (t3 <= t4)) 
                {
//			        if (env0_time_end == NULL)  retval = 0;
//			        else 
                    retval = 3;
                }
                else
                {
                    if ((t4 <= t1) && (t4 <= t) && (t4 <= t2) && (t4 <= t3)) 
                    {
//			        if (env1_time_end == NULL)  retval = 0;
//			        else 
                        retval = 4;
                    }
                    else
                    {
                        retval = 0;
                    }
                }
            }
        }
    }
    return retval;
}

/*
 * end transcription: update the mRNAs in the nucleus, cytoplasm
 * etc. accordingly and delete the event from the queue
 */
void end_transcription(float *dt,
                       float t,
                       CellState *state,
                       GillespieRates *rates,
                       int ngenes)
{
    int i, j, total;
    
//    LOG_ERROR("END TRANSCRIPTION\n");  
    
    /* recompute the delta-t based on difference between now and the
       time of transcription end */
    *dt = state->mRNA_transcr_time_end->time - t;

    if (verbose) 
    {
      total = 0;
      for (i=0; i < ngenes; i++) 
        for (j=0; j < MAX_COPIES; j++) 
          total += state->mRNA_transcr_num[i][j];

//      LOG_VERBOSE("\ntranscription event finishes out of %d possible t=%g dt=%g\n", total, t, *dt);
    }

    /* get the gene which is ending transcription */
    i = state->mRNA_transcr_time_end->gene_id;
    j = state->mRNA_transcr_time_end->copy;

    /* increase number of mRNAs in nucleus */
    (state->mRNA_nuclear_num[i])++;

    /* decrease the number of mRNAs undergoing transcription */
    (state->mRNA_transcr_num[i][j])--;

    /* delete the fixed even which has just occurred */
    delete_fixed_event_start(&(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));

    /* add rate KRNA to transport and Gillespie rates */ // update rates in calc_all_rates
  //  transport[i] += KRNA;
  //  rates->transport += KRNA;
  //  rates->transport_operations++;

  //  LOG_VERBOSE("add one new mRNA in nucleus, updating transport[%d]=%g, rates->transport=%g\n", i, transport[i], rates->transport);

}

void end_translation(Genotype *genotype, CellState *state, GillespieRates *rates, float *dt,float t)
{
    int i;
    int N_translation_event=0;
    
    *dt = state->mRNA_transl_time_end->time - t;         /* make dt window smaller */
    
    /* count current number of mRNAs that have recently arrived in cytoplasm */
    for (i=0; i<genotype->ngenes; i++) N_translation_event += state->mRNA_transl_cyto_num[i];

//    LOG_VERBOSE("translation event finishes out of %d possible t=%g dt=%g\n", N_translation_event, t, *dt); /* bug: dt can be negative */

    /* get identity of gene that has just finished translating */
    i=state->mRNA_transl_time_end->gene_id;   

    /* there is one less mRNA that has just finished translation */
    (state->mRNA_transl_cyto_num[i])--;   
    
    change_mRNA_cytoplasm(state->mRNA_transl_time_end->gene_id, genotype, state, rates); 

    /* delete the event that just happened */
//    LOG_VERBOSE("delete translation event that just happened at time=%g", t);    
    
    delete_fixed_event_start(&(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));

    /* there is one more mRNA that is now in cytoplasm */
    (state->mRNA_cyto_num[i])++;
}

/*
 * do the actual disassembling of the pre-initiation complex
 */
//void disassemble_PIC(CellState *state,
//                     Genotype *genotype,
//                     int gene_id,
//                     int gene_copy,
//                     GillespieRates *rates)
//{
//  float disassembly = rates->pic_disassembly_rate[gene_id];
//  LOG_ERROR("PIC DIS2\n");
// 
//  rates->pic_disassembly -= disassembly;
//  rates->pic_disassembly_operations++;
//  
//  state->active[gene_id][0]= NO_NUC_NO_PIC; // GO WITH UPDATE ALL RATES AT ONCE IN CALC_ALL_RATES
//}

/*
 * time course of [TF]s represented as array of TimeCourse lists.
 */
//void add_time_points(float time,
//                     float protein_conc[NPROTEINS],
//                     TimeCourse **timecoursestart,
//                     TimeCourse **timecourselast)
//{
//  int i;
//  
//  for (i=0; i < NPROTEINS; i++)
//    add_time_point(time, protein_conc[i], &(timecoursestart[i]), &(timecourselast[i]));
//}

//void add_integer_time_points(float time,
//                             int protein_conc[NPROTEINS],
//                             TimeCourse **timecoursestart,
//                             TimeCourse **timecourselast)
//{
//  int i;
//  
//  for (i=0; i < NPROTEINS; i++)
//    add_time_point(time, (float) protein_conc[i], &(timecoursestart[i]), &(timecourselast[i]));
//}

/*
 * compute tprime factor used in the integration of growth rate
 */
float compute_tprime(float c, 
                     float P, 
                     float alpha, 
                     float s_mRNA) 
{
  return (1/c) * log((c*P - alpha*s_mRNA)/(c*P - alpha*s_mRNA));
}

/*
 * get integral for growth rate
 */
float compute_integral(float alpha, 
                       float c, 
                       float gmax, 
                       float deltat, 
                       float s_mRNA, 
                       float P, 
                       float Pp,                       
                       float ect1) 
{
  return 1.0/(pow(c,2)*Pp) * gmax * (-alpha*ect1*s_mRNA + c*(P*ect1 + alpha*deltat*s_mRNA));
}

/*
 * return the instantaneous growth rate given the current cell state,
 * also return the integrated growth rate as a pointer
 */
float compute_growth_rate_dimer(float *integrated_growth_rate,
                                Genotype *genotype,
                                CellState *state,
                                float P_a,
                                float P_b,
                                float t, 
                                float dt,                 
                                float ect1_a,                 
                                float ect1_b,
				int env)
{
    int i;
    float instantaneous_growth_rate=0.0;  /* this is returned from the function */
    float total_alpha_s = 0.0;
    float dt_prime, dt_rest;
    
    int selection_gene_A=genotype->nproteins-2;
    int selection_gene_B=genotype->nproteins-1;

    float P_next_a=state->protein_conc[selection_gene_A];
    float P_next_b=state->protein_conc[selection_gene_B];
  

//  LOG_VERBOSE_NOCELLID("P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//                       P, P_next, c, t, t+deltat, s_mRNA);
                       
    switch (env)
    {   
        case 1: /* protein A is necessary!*/
            /* choose the appropriate piecewise linear integral */           
            if ( ((P_a >= Pp_a) && (P_next_a >= P_a)) || ((P_next_a >= Pp_a) && (P_a >= P_next_a)) )  /* P > Pp throughout */
            {          
//                if (verbose)
//                  LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);

                *integrated_growth_rate = gmax_a * dt;	    
                     
            }
            else if (((P_next_a <= Pp_a) && (P_next_a >= P_a)) || ((P_a <= Pp_a) && (P_a >= P_next_a))) /* P < Pp throughout */
            {   
  //		    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
                *integrated_growth_rate = compute_integral(genotype->translation[selection_gene_A], 
                                                            state->konvalues[selection_gene_A][KON_PROTEIN_DECAY_INDEX], 
                                                            gmax_a, 
                                                            dt,
                                                            state->mRNA_cyto_num[selection_gene_A],
                                                            P_a, Pp_a, ect1_a);
            }
            else if ((Pp_a > P_a) && (P_next_a > Pp_a)) /* P < Pp up until t' then P > Pp */
            {    
                dt_prime = compute_tprime(state->konvalues[selection_gene_A][KON_PROTEIN_DECAY_INDEX], 
                                          P_a, 
                                          genotype->translation[selection_gene_A], 
                                          state->mRNA_cyto_num[selection_gene_A]);
                
                dt_rest = dt - dt_prime;
//		LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
                
                *integrated_growth_rate = compute_integral(genotype->translation[selection_gene_A], 
                                                            state->konvalues[selection_gene_A][KON_PROTEIN_DECAY_INDEX], 
                                                            gmax_a, 
                                                            dt_prime, 
                                                            state->mRNA_cyto_num[selection_gene_A], 
                                                            P_a, 
                                                            Pp_a, 
                                                            ect1_a);
                
                *integrated_growth_rate += gmax_a * dt_rest;
            }
            else if ((P_a > Pp_a) && (Pp_a > P_next_a)) /* P > Pp up until t' then P < Pp */
            {   
                dt_prime = compute_tprime(state->konvalues[selection_gene_A][KON_PROTEIN_DECAY_INDEX], 
                                          P_a, 
                                          genotype->translation[selection_gene_A], 
                                          state->mRNA_cyto_num[selection_gene_A]);
                
                dt_rest = dt - dt_prime;
//		    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
                
                *integrated_growth_rate = gmax_a * dt_prime;
                
                *integrated_growth_rate += compute_integral(genotype->translation[selection_gene_A], 
                                                            state->konvalues[selection_gene_A][KON_PROTEIN_DECAY_INDEX], 
                                                            gmax_a, 
                                                            dt_rest, 
                                                            state->mRNA_cyto_num[selection_gene_A], 
                                                            P_a, 
                                                            Pp_a, 
                                                            ect1_a);
            }
            else 
            {
//              LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//		                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
//              LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
//		                       cell_id);
                exit(1);
            }

            *integrated_growth_rate-= penalty*compute_integral(genotype->translation[selection_gene_B], 
                                                                state->konvalues[selection_gene_B][KON_PROTEIN_DECAY_INDEX], 
                                                                1.0, 
                                                                dt, 
                                                                state->mRNA_cyto_num[selection_gene_B], 
                                                                P_b, 
                                                                1.0, 
                                                                ect1_b);
              /* compute instantaneous growth rate at t */
            if (P_next_a < Pp_a)
                instantaneous_growth_rate = gmax_a*P_next_a/Pp_a - penalty*P_next_b;
            else
                instantaneous_growth_rate = gmax_a - penalty*P_next_b;

            break;
    
        case 0: /* protein b is necessary! */

            if (((P_b >= Pp_b) && (P_next_b >= P_b)) || ((P_next_b >= Pp_b) && (P_b >= P_next_b))) /* P > Pp throughout */
            {          
    //          if (verbose)
    //		      LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);
                *integrated_growth_rate = gmax_b * dt;

            }
            else if (((P_next_b <= Pp_b) && (P_next_b >= P_b)) || ((P_b <= Pp_b) && (P_b >= P_next_b))) /* P < Pp throughout */
            {   
    //		    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
                *integrated_growth_rate = compute_integral(genotype->translation[selection_gene_B], 
                                                            state->konvalues[selection_gene_B][KON_PROTEIN_DECAY_INDEX],
                                                            gmax_b, 
                                                            dt, state->mRNA_cyto_num[selection_gene_B], 
                                                            P_b, Pp_b, ect1_b);
            }
            else if ((Pp_b > P_b) && (P_next_b > Pp_b)) /* P < Pp up until t' then P > Pp */
            {    
                dt_prime = compute_tprime(state->konvalues[selection_gene_B][KON_PROTEIN_DECAY_INDEX], 
                                          P_b, 
                                          genotype->translation[selection_gene_B], 
                                          state->mRNA_cyto_num[selection_gene_B]);

                dt_rest = dt - dt_prime;
    //		    LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
    //		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);

                *integrated_growth_rate = compute_integral(genotype->translation[selection_gene_B], 
                                                            state->konvalues[selection_gene_B][KON_PROTEIN_DECAY_INDEX], 
                                                            gmax_b, 
                                                            dt_prime, 
                                                            state->mRNA_cyto_num[selection_gene_B], 
                                                            P_b, 
                                                            Pp_b, 
                                                            ect1_b);

                *integrated_growth_rate += gmax_b * dt_rest;
            }
            else if ((P_b > Pp_b) && (Pp_b > P_next_b)) /* P > Pp up until t' then P < Pp */
            {   
                dt_prime = compute_tprime(state->konvalues[selection_gene_B][KON_PROTEIN_DECAY_INDEX], 
                                          P_b, 
                                          genotype->translation[selection_gene_B], 
                                          state->mRNA_cyto_num[selection_gene_B]);

                dt_rest = dt - dt_prime;
    //		    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
    //		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);

                *integrated_growth_rate = gmax_b * dt_prime;

                *integrated_growth_rate += compute_integral(genotype->translation[selection_gene_B], 
                                                            state->konvalues[selection_gene_B][KON_PROTEIN_DECAY_INDEX], 
                                                            gmax_b, 
                                                            dt_rest, 
                                                            state->mRNA_cyto_num[selection_gene_B], 
                                                            P_b, 
                                                            Pp_b, 
                                                            ect1_b);
            }
            else 
            {
    //          LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
    //		                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
    //          LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
    //		                       cell_id);
                exit(1);
            }

            *integrated_growth_rate -= penalty*compute_integral(genotype->translation[selection_gene_A], 
                                                                state->konvalues[selection_gene_A][KON_PROTEIN_DECAY_INDEX], 
                                                                1.0, 
                                                                dt, 
                                                                state->mRNA_cyto_num[selection_gene_A], 
                                                                P_a, 
                                                                1.0, 
                                                                ect1_a);	
            /* compute instantaneous growth rate at t */
            if (P_next_b < Pp_b)
                instantaneous_growth_rate = gmax_b*P_next_b/Pp_b - penalty*P_next_a;
            else
                instantaneous_growth_rate = gmax_b - penalty*P_next_a;

            break;
    }

//  LOG_VERBOSE_NOCELLID("growth rate (variable %g)-", *integrated_growth_rate);

    /* compute the total cost of translation across all genes  */
    for (i=0; i < genotype->ngenes; i++) 
    {
        total_alpha_s += genotype->translation[i] * state->mRNA_cyto_num[i];
    }

    /* add constant term for integrated rate */
    *integrated_growth_rate += -h * dt * (total_alpha_s);

    /* and instantaneous integrated rate */
    instantaneous_growth_rate += -h * (total_alpha_s);

  //  LOG_VERBOSE_NOFUNC("(constant %g) = (total %g)\n", (h*deltat*total_alpha_s), *integrated_growth_rate);

    /* make sure growth rates can't be negative */
    if (*integrated_growth_rate < 0.0)
        *integrated_growth_rate = 0.0;

    if (instantaneous_growth_rate < 0.0)
        instantaneous_growth_rate = 0.0;

  #if 0 /* currently disable printing out growth rate information */
    fprintf(fp_growthrate[0], "%g %g %g %g %g %g\n", t, instantaneous_growth_rate, *integrated_growth_rate, P, s_mRNA, c);
  #endif

    /* return the instantaneous growth rate */
    return (instantaneous_growth_rate);
}

/* 
 * update both the protein concentration and current cell size
 *
 * need some sort of control in case it declines to essentially zero.
 * Add in discrete, stochastic and/or zero values, but this may create
 * false attractor if time steps are small and rising tide keeps
 * getting rounded down
 */
void update_protein_conc_cell_size( Genotype *genotype,
                                    CellState *state,                                   
                                    GillespieRates *rates,
                                    float dt,                                                                      
                                    float t,                                  
                                    int env)
{
    int i,j;
    float ct, ect, ect1,ect1_a,ect1_b;
    float L_a, L_b;
    float instantaneous_growth_rate = 0.0;
    float integrated_growth_rate = 0.0;
    float adjusted_decay;

    /* store the concentration of the selection genes before updating*/
    L_a=0.0;    
    for(j=0;j<genotype->protein_pool[genotype->nproteins-2][0][0];j++)
        L_a+=state->gene_specific_protein_conc[genotype->protein_pool[genotype->nproteins-2][1][j]];
    
    L_b=0.0;
    for(j=0;j<genotype->protein_pool[genotype->nproteins-1][0][0];j++)
        L_b+=state->gene_specific_protein_conc[genotype->protein_pool[genotype->nproteins-1][1][j]];
    
    for (i=0; i < genotype->ngenes; i++) 
    {
        /* update protein decay rates due to dilution caused by growth */
        adjusted_decay = genotype->proteindecay[i] + state->growth_rate;

        /* if this results in a very small or zero decay rate, use protein aging term */
        /* NOTE: we need to update this for every gene that encodes this protein*/
        if (adjusted_decay > protein_aging)
            state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = adjusted_decay;
        else 
            state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = protein_aging;

      /* print out warning if decay rates get too low */
  //    if (state->konvalues[i][KON_PROTEIN_DECAY_INDEX] < 1e-10) {
  //      LOG_WARNING("protein=%02d, protein_decay=%g, genotype->proteindecay=%g, protein_aging=%g\n", i, adjusted_decay, 
  //                  genotype->proteindecay[i], protein_aging);
  //    }

  //    LOG_VERBOSE("prot decay[%d]=%g\n", i, state->konvalues[i][KON_PROTEIN_DECAY_INDEX]);

        ct = state->konvalues[i][KON_PROTEIN_DECAY_INDEX]*dt;
        ect = exp(-ct);
        if (fabs(ct)<EPSILON) ect1=ct;
        else ect1 = 1-ect; 

        /* record the ect for the selection genes*/
        if (i == genotype->ngenes-2) ect1_a=ect1;		
        if (i == genotype->ngenes-1) ect1_b=ect1;	

        /* get the new protein concentration for this gene */
        state->gene_specific_protein_conc[i]=ect*state->gene_specific_protein_conc[i]+state->konvalues[i][KON_SALPHC_INDEX]*ect1 ;        
    }
    
    /* now, use protein_pool to pool gene specific concentration*/
    for(i=0;i<genotype->nproteins;i++)
    {
        state->protein_conc[i]=0.0;
        
        for(j=0;j<genotype->protein_pool[i][0][0];j++)
            state->protein_conc[i]+=state->gene_specific_protein_conc[genotype->protein_pool[i][1][j]];
    }
    
    /* now find out the protein concentration at end of dt interval and compute growth rate */   
    instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genotype, 
                                                            state, 
                                                            L_a, 
                                                            L_b, 
                                                            t, 
                                                            dt, 
                                                            ect1_a,
                                                            ect1_b,
                                                            env);
    
    /* use the integrated growth rate to compute the cell size in the next timestep */
    state->cell_size = (state->cell_size)*exp(integrated_growth_rate);
    
    /* update the instantaneous growth rate for the beginning of the next timestep */
    state->growth_rate = instantaneous_growth_rate; 
//    fprintf(fp_cellsize[0], "%g %g\n", t, state->cell_size);
//
//    if ((output) && (*timecourselast)->time < t+dt-0.1) 
//        add_time_points(t+dt, otherdata, timecoursestart, timecourselast);
}

/*
 * START
 * Functions that handle each possible Gillespie event 
 *
 */
void transport_event(Genotype *genotype,
                     CellState *state,
                     GillespieRates *rates,
                     float dt,
                     float t,
                     long int *seed)
{
    int gene_id=-1;
    float x;
    float endtime = t + dt + TTRANSLATION;

    x=ran1(seed)*rates->transport; 

    /* choose gene product (mRNA) that gets transported to cytoplasm
       based on weighting in transport[] array */
    while (gene_id < genotype->ngenes-1 && x > 0.0) 
    {
        gene_id++;
        x -= rates->transport_rate[gene_id];
    }

    if (gene_id >= genotype->ngenes) 
    {
//        LOG_ERROR("[cell %03d] attempted to choose mRNA for gene=%d which doesn't exist\n", state->cell_id, gene_id);
        exit(0);
    } 

//    LOG_VERBOSE("do transport event mRNA from gene=%d from %d copies (x=%g)\n", gene_id, state->mRNA_nuclear_num[gene_id], x);

    (state->mRNA_nuclear_num[gene_id])--;   /* one less mRNAs in nucleus */
    (state->mRNA_transl_cyto_num[gene_id])++;   /* it has just arrived in cytoplasm, ready to be translated */

    /* add the endtime for translation */
//    LOG_VERBOSE("add translation event endtime=%f for mRNA encoded by gene=%d \n", endtime, gene_id);
    
    add_fixed_event_end(gene_id, -1, endtime, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));

//    rates->transport_rate[i] -= KRNA;   /* decrease transport frequency */ // UPDATE all rates in calc_all_rates
//
//    /* if below a threshold, make zero */
//    if (rates->transport_rate[i] < 0.1*KRNA) 
//      rates->transport_rate[i]=0.0;
//
//    /* adjust rates */
//    rates->transport -= KRNA;
//    rates->transport_operations++;
//
//    /* do similar threshold check */
//    if (rates->transport < 0.1*KRNA) 
//      rates->transport=0.0;
}

void mRNA_decay_event(GillespieRates *rates, CellState *state, Genotype *genotype, long int *seed)
{
    int gene_id = -1,j;
    float x=ran1(seed)*rates->mRNAdecay;

    /* loop through mRNA products, to choose the mRNA with the
       proportionally higher decay rate */
    while (gene_id < genotype->ngenes-1 && x > 0.0) 
    {
        gene_id++;
        x-= rates->mRNAdecay_rate[gene_id];
    }
    
//    if (x > temp_rate) 
//    { /* JM: had some rounding errors with rates->mRNAdecay. Calculate in calc_dt, hopefully fixed now */
//          LOG_WARNING("x=%g > temp_rate=%g out of rates->mRNAdecay=%g\n",
//                      x, temp_rate, rates->mRNAdecay);
//    }

    /* assume mRNA cytoplasm transport events equally likely */
    x = ran1(seed)*((float) (state->mRNA_cyto_num[gene_id] + state->mRNA_transl_cyto_num[gene_id]));
    
    /* decay mRNA in cytoplasm */
    if (x < (float)state->mRNA_cyto_num[gene_id]) 
    {
//        LOG_VERBOSE("mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
//                    gene_id, state->mRNA_cyto_num[gene_id], state->mRNA_transl_cyto_num[gene_id]);

        /* remove the mRNA from the cytoplasm count */
        (state->mRNA_cyto_num[gene_id])--;  
        
        change_mRNA_cytoplasm(gene_id, genotype, state, rates); 

    } 
    else 
    {
        /* decay mRNA in process of translating */
        x = ran1(seed)*((float) state->mRNA_transl_cyto_num[gene_id]);
        
//        LOG_VERBOSE("mRNA decay event gene %d not from %d copies in cytoplasm but %f from %d copies translating\n",
//                    gene_id, state->mRNA_cyto_num[gene_id], trunc(x), state->mRNA_transl_cyto_num[gene_id]);

        /* delete this fixed event: this mRNA will never be translated */
//        LOG_VERBOSE("delete fixed TRANSLATION EVENT at time =%d for gene=%d\n", (int) trunc(x), gene_id);
        delete_fixed_event(gene_id, -1, (int) trunc(x), &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));

        /* remove the mRNA from the count */
        (state->mRNA_transl_cyto_num[gene_id])--;
        
//        if (verbose) 
//            for (j=0; j < NGENES; j++) 
//            {
//                LOG_VERBOSE("%d copies of gene %d translating\n", state->mRNA_transl_cyto_num[j], j);
//            }
    }
}

void histone_acteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype, long int *seed)
{
    int gene_id=-1;
    float x = ran1(seed)*rates->acetylation;

    while(gene_id<genotype->ngenes-1 && x>0.0)
    {
        gene_id++;
        x-=rates->acetylation_rate[gene_id];
    }

  //  LOG_VERBOSE("acetylation event gene %d (copy %d)\nstate change from %d to 4\n",
  //              gene_id, 0, state->active[gene_id][0]);

  //  if (state->active[gene_id][0] != NUC_NO_PIC) {
  //    LOG_ERROR("acetylation event on gene %d (copy %d) attempted from state %d\n", gene_id, 0, state->active[gene_id][0]);
  //  }

    /* set state: eject nucleosome, but there is no PIC yet */
    state->active[gene_id][0] = NO_NUC_NO_PIC;
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype, long int *seed)
{
    int gene_id=-1; 
    float x = ran1(seed)*rates->deacetylation;

    /* choose a particular gene and copy to change state */
    while(gene_id<genotype->ngenes-1 && x>0.0)
    {
        gene_id++;
        x-=rates->deacetylation_rate[gene_id];
    }

  //  LOG_VERBOSE("deacetylation event gene %d (copy %d)\nstate change from %d to 1\n",
  //              gene_id, 0, state->active[gene_id][0]);
  //  
  //  if (state->active[gene_id][0] != NO_NUC_NO_PIC) {
  //    LOG_ERROR("deacetylation event attempted from state %d\n", state->active[gene_id][0]);
  //  }

    /* set state: nucleosome returns */
    state->active[gene_id][0] = NUC_NO_PIC;
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genotype, long int *seed)
{
    float x = ran1(seed)*rates->pic_assembly;

    int gene_id=-1; 

    /* choose a particular gene and copy to change state */
    while(gene_id<genotype->ngenes-1 && x>0.0)
    {
        gene_id++;
        x-=rates->pic_assembly_rate[gene_id];
    }

  //  LOG_VERBOSE("PIC assembly event gene %d copy %d\nstate change from %d to 6\n",
  //              gene_id,0, state->active[gene_id][0]);
  //
  //  if (state->active[gene_id][0] != NO_NUC_NO_PIC) {
  //    LOG_ERROR("PIC assembly event attempted from state %d\n", state->active[gene_id][0]);
  //  }

    /* turn gene fully on: ready for transcription and adjust rates */
    state->active[gene_id][0] = PIC_NO_NUC;

  //  (rates->transcript_init_num[0])++; //update rates in calc_all_rates

  //  rates->pic_disassembly_operations++;
}

void disassemble_PIC_event(Genotype *genotype, CellState *state,GillespieRates *rates, long int *seed)
{
    int gene_id=-1;
    float x=ran1(seed)*rates->pic_disassembly;

    /* choose an appropriate gene copy to disassemble the PIC from */
    while (gene_id < genotype->ngenes-1 && x>0) 
    {
        gene_id++;       
        x -= rates->pic_disassembly_rate[gene_id];
    }
    
//    if (gene_id>=NGENES*current_ploidy) { LOG_ERROR("error in PIC disassembly\n"); }
    
    /* now get the gene_id */     
//    LOG_VERBOSE("PIC disassembly event in copy %d of gene %d\n", gene_copy, gene_id);
    /* do the disassembling */
//    disassemble_PIC(state, genotype, gene_id, gene_copy, rates);
    
    state->active[gene_id][0]=NO_NUC_NO_PIC; // UPDATE ALL RATES TOGETHER IN CALC_ALL_RATES    
}

void transcription_init_event(GillespieRates *rates, CellState *state, Genotype *genotype, float dt, float t, long int *seed)
{
    int gene_id=-1;  
    int x=ran1(seed)*rates->transcript_init;

    while(gene_id<genotype->ngenes-1 && x>0)
    {
        gene_id++;
        x-=rates->transcript_init_rate[gene_id];
    }

  //  LOG_VERBOSE("transcription event gene %d, copy %d\n", gene_id, 0);
  //  LOG_ERROR("Gene id = %d\n", gene_id);
  //
  //  if (state->active[gene_id][0] != PIC_NO_NUC) {
  //    LOG_ERROR("transcription event attempted from state %d\n", state->active[gene_id][0]);
  //  }                           

    /* now that transcription of gene has been initiated, 
     * we add the time it will end transcription, 
     * which is dt+time of transcription from now */
    add_fixed_event_end(gene_id, 0, t+dt+TTRANSCRIPTION, 
                        &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));

    /* increase the number mRNAs being transcribed */
    (state->mRNA_transcr_num[gene_id][0])++;                      
}
/*
 * END
 * Functions that handle each possible Gillespie event 
 */

/*
 * given a rate array and position within the array, get the the gene
 * copy and location
 */

void clone_cell(Genotype *genotype_orig,                
                Genotype *genotype_clone,
                int clone_type)
{
    int i, j;
        
    if(clone_type!=4) /* not a mutation in rate constant*/
    {
        for(i=0; i< genotype_orig->ngenes;i++)
        {
            if(genotype_clone->re_calc[i][3]  /* only copy to places that mutated*/
               || genotype_orig->re_calc[i][3]) /* this argument is used when copy new genotype to current genotype, in try_fixation*/                
            {
                /* we copy binding sites info here, because recalc these info is more expensive*/
                for(j=0;j<genotype_orig->binding_sites_num[i];j++)
                {
                    genotype_clone->all_binding_sites[i][j].tf_id=genotype_orig->all_binding_sites[i][j].tf_id;
                    genotype_clone->all_binding_sites[i][j].Koff=genotype_orig->all_binding_sites[i][j].Koff;
                    genotype_clone->all_binding_sites[i][j].BS_pos=genotype_orig->all_binding_sites[i][j].BS_pos;
                    genotype_clone->all_binding_sites[i][j].N_hindered=genotype_orig->all_binding_sites[i][j].N_hindered;
                }
                genotype_clone->binding_sites_num[i]=genotype_orig->binding_sites_num[i];
                genotype_clone->max_hindered_sites[i]=genotype_orig->max_hindered_sites[i];
                genotype_clone->N_act_BS[i]=genotype_orig->N_act_BS[i];
                genotype_clone->N_rep_BS[i]=genotype_orig->N_rep_BS[i];
        
//                for(j=0;j<genotype_orig->max_N_rep_bound[i];j++)
//                {
//                    genotype_clone->N_configurations[i][j]=genotype_orig->N_configurations[i][j];
//                }        
//                genotype_clone->max_N_rep_bound[i]=genotype_orig->max_N_rep_bound[i];
//                genotype_clone->max_N_act_bound[i]=genotype_orig->max_N_act_bound[i];
                
                if(clone_type!=3) /* if the mutation was not in binding sequence */
                {
                    memcpy(&genotype_clone->cisreg_seq[i][0][0],&genotype_orig->cisreg_seq[i][0][0],CISREG_LEN*sizeof(char));
                }
            }
        }
    }
    
    if(clone_type!=1) /* not a substitution or indel*/
    {
        for (i=0; i < genotype_orig->ngenes; i++) 
        {              
            genotype_clone->mRNAdecay[i]=genotype_orig->mRNAdecay[i];
            genotype_clone->proteindecay[i]=genotype_orig->proteindecay[i];
            genotype_clone->translation[i]=genotype_orig->translation[i];
            genotype_clone->re_calc[i][0]=genotype_orig->re_calc[i][0];
            genotype_clone->re_calc[i][1]=genotype_orig->re_calc[i][1];
            genotype_clone->pic_disassembly[i][0]=genotype_orig->pic_disassembly[i][0];
            genotype_clone->which_protein[i]=genotype_orig->which_protein[i];
        }
        
        for(i=0;i<genotype_orig->nproteins;i++)
        {
            genotype_clone->activating[i][0]= genotype_orig->activating[i][0];  
            genotype_clone->protein_pool[i][0][0]=genotype_orig->protein_pool[i][0][0];

            for(j=0;j<MAXALLOC;j++)
            {
                genotype_clone->protein_pool[i][1][j]=genotype_orig->protein_pool[i][1][j];
            }         
        }
    }
    
    /* since there is no tag to mark which binding seq is mutated, we copy all*/
    for (i=0; i < genotype_orig->ntfgenes; i++) 
    {        
        for(j=0;j<TF_ELEMENT_LEN;j++)
        {    
            genotype_clone->tf_seq[i][0][j]=genotype_orig->tf_seq[i][0][j];
            genotype_clone->tf_seq_rc[i][0][j]=genotype_orig->tf_seq_rc[i][0][j];
        }
    }

    /* these are easy, so do it everytime*/
    genotype_clone->fitness=genotype_orig->fitness;
    genotype_clone->ngenes=genotype_orig->ngenes;
    genotype_clone->ntfgenes=genotype_orig->ntfgenes;
    genotype_clone->nproteins=genotype_orig->nproteins;
    genotype_clone->N_act=genotype_orig->N_act;
    genotype_clone->N_rep=genotype_orig->N_rep;
         
    for(i=0;i<genotype_orig->ngenes;i++)
    {
        genotype_clone->re_calc[i][2]=0;   /* do not recalc binding sites unless mutation changes it*/      
    }
    
    if(clone_type==6)/* if it is a fixation event*/
    {
        for(i=0;i<genotype_orig->ngenes;i++)
            genotype_clone->re_calc[i][3]=1;   
    }
    else
    {
        for(i=0;i<genotype_orig->ngenes;i++)
            genotype_clone->re_calc[i][3]=0;   /* do not copy info for this gene, unless mutation changes it*/
    }
}

/*
 * diagnostic function to dump a copy of the GillespieRates and some
 * of the CellState to the current error file
 */
void log_snapshot(Genotype *genotype, 
                  CellState *state,
                  GillespieRates *rates,                 
                  float x,
                  float t)
{
//  int i, p, nkon = 0;
//
//  LOG("snapshot at time=%g:\n x=%g, koff=%g = %d (tf_bound_num) * %g (koff/tf_bound_num)\n transport=%g\n decay=%g\n",
//      t, x, rates->koff, state->tf_bound_num, rates->koff/(float)state->tf_bound_num, 
//      rates->transport, rates->mRNAdecay);
//  LOG_NOFUNC(" rates->salphc=%g\n rates->max_salphc=%g rates->min_salphc=%g\n", rates->salphc, rates->max_salphc, rates->min_salphc);
//  LOG_NOFUNC(" konrate=%g\n", konrate);
//  LOG_NOFUNC(" pic_disassembly=%g\n kon=%g = %d * %g\n",
//             rates->pic_disassembly, rates->salphc+(konrate), kon_states->nkon, (rates->salphc+(konrate))/(float)kon_states->nkon);
//  
//  for (p=0; p < MAX_COPIES; p++) {
//    LOG_NOFUNC(" acetylation=%g (copy %d)\n deacetylation=%g (copy %d)\n PIC assembly=%g (copy %d)\n transcriptinit=%g (copy %d)\n",
//               (float)rates->acetylation_num[p]*ACETYLATE, p, (float)rates->deacetylation_num[p]*DEACETYLATE, p, 
//               (float)rates->pic_assembly_num[p]*PICASSEMBLY, p, (float)rates->transcript_init_num[p]*TRANSCRIPTINIT, p);
//  }
//  LOG_NOFUNC(" total rates=%g=%g+%g\n", rates->subtotal + (konrate), rates->subtotal, konrate);
//  LOG_NOFUNC(" total free=%d + total bound=%d + total hindered=%d = total sites=%d\n", 
//             kon_states->nkon, state->tf_bound_num, state->tf_hindered_num, genotype->binding_sites_num);
//  
//  LOG_NOFUNC("\n");
}

/*
 * run the model for a specified cell for a single timestep:
 *  - returns 0  if cell is not "dead" (i.e. rates haven't deteroriated to zero)
 *  - returns -1 if cell is "dead"
 */
void do_single_timestep(Genotype *genotype, 
                       CellState *state,                         
                       GillespieRates *rates, 
                       float *t,          
                       int maxbound2,
                       int maxbound3,                                      
		       int *env,
                       long int *seed) 
{    
    int event;     /* boolean to keep track of whether FixedEvent has ended */   
    float fixed_time; 
    float dt;
    float x;

  //  if (verbose) //VERBOSE JU CHANGE
  //    for (j=0; j < MAX_COPIES; j++) {
  //      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
  ////      LOG_ERROR("FDSTS rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
  //    }
  //    
  //    for (j=0; j < MAX_COPIES; j++) {
  //      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
  //      //LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
  ////       LOG_ERROR("rates->deacetylation_num[%d]=%d\n", j, rates->deacetylation_num[j]);
  //    }

  #if 0 /* currently disable printing out the TF occupancy and amount of rounding */
    print_tf_occupancy(state, genotype->all_binding_sites, *t);
    print_rounding(state, rates, *t);
  #endif

    x = expdev(seed);        /* draw random number */

    dt = x/rates->subtotal;

    if (dt < 0.0) {
//      LOG_ERROR("dt=%g is negative after calc_dt, t=%g\n", *dt, *t);
      exit(-1);
    }
    
//    LOG_VERBOSE("next stochastic event due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

    /* first check to see if a fixed event occurs in current t->dt window, or in tdevelopment if running for a fixed development time */
    fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;

    event = does_fixed_event_end(state->mRNA_transl_time_end,
                                 state->mRNA_transcr_time_end,
                                 state->env0_time_end,
                                 state->env1_time_end,
                                 fixed_time);
    while(event>0)
    {                
        do_fixed_event(genotype, state, rates, &dt, *t, event, env);
        
        update_protein_conc_cell_size(genotype, state, rates, dt, *t, *env);
    
        *t += dt;                  /* advance time by the dt */
        
        x -= dt*rates->subtotal;  /* we've been running with rates->subtotal for dt, so substract it from x*/

        /* update rates->subtotal and re-compute a new dt */
        calc_all_rates(genotype, state, rates, NO_KON_UPDATE);      
        
        dt = x/rates->subtotal;       
        
//        LOG_VERBOSE("next stochastic event (2) due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

        fixed_time = (*t+dt<tdevelopment)?(*t+dt):tdevelopment;    

        /* check to see there aren't more fixed events to do */
        event = does_fixed_event_end(state->mRNA_transl_time_end, 
                                     state->mRNA_transcr_time_end,
                                     state->env0_time_end,
                                     state->env1_time_end,
                                     fixed_time);
    } 

  /* no remaining fixed events to do in dt, now do stochastic events */
  
  /* if we haven't already reached end of development with last
     delta-t, if there is no fixed development time, we always execute
     this  */
          
    if (*t+dt < tdevelopment)
    { 
      /* if the total rates falls below zero, we do an emergency recalibration of cell */
//        if (!(rates->subtotal > 0.0)) 
//        {
////            log_snapshot(genotype, state, rates, x, *t);
////            LOG_ERROR("x should always be >0 t=%g (x=%g) rates->subtotal=%g, recalibrate cell!\n", *t, *x, rates->subtotal); 
//            calc_all_rates(genotype, state, rates, UPDATE_ALL);
////            log_snapshot(genotype, state, rates, x, *t);
//
//            /* if this still results in either zero or negative total rates,
//               this most likely due the cell being "dead" no TFs bound, no
//               activity etc.  We mark cell as "dead" in this case, and
//               remove from queue. */
//            if (!(rates->subtotal > 0.0)) 
//            {  
////                  LOG_ERROR("cell is effectively dead\n"); 
//                  return -1;        /* return cell status as "dead" */
//            }
//        }  

        do_Gillespie_event(genotype, state, rates, dt, *t, seed);

        update_protein_conc_cell_size(genotype, state, rates, dt, *t, *env);
        
        calc_all_rates(genotype,state,rates,NO_KON_UPDATE);
        
        /* Gillespie step: advance time to next event at dt */
        *t += dt;
//        LOG_VERBOSE("dt=%g t=%g\n", *dt, *t);
    } 
    else 
    { 
        /* we will reach the end of development in dt */
//        LOG_VERBOSE("finish at t=%g dt=%g\n", *t, *dt);

        /* do remaining dt */
        dt = tdevelopment - *t;

        /* final update of protein concentration */
        update_protein_conc_cell_size(genotype, state, rates, dt, *t, *env);          
                                       
        /* advance to end of development (this exits the outer while loop) */
        *t = tdevelopment;
    }
    
//    return 0;
}

/* while there are either transcription or translation events
     occuring in current t->dt window */
void do_fixed_event(Genotype *genotype, 
                   CellState *state, 
                   GillespieRates *rates, 
                   float *dt,
                   float t,
                   int event, 
                   int *env)
{  
    switch (event) 
    {
        case 1:     /* if a transcription event ends */
            end_transcription(dt, t, state, rates, genotype->ngenes);                 
            break;

        case 2:     /* if a translation event ends */ 
            end_translation(genotype, state, rates, dt, t);                            
            break;

        case 3:     /*change env from 0 to 1*/  
            state->cumulative_growth_rate*=log(state->cell_size/state->cell_size_copy)/duration_env0;
            state->cell_size_copy=state->cell_size;
            state->env_change++;
            *dt = state->env0_time_end->time - t;
            *env = 1;
            delete_fixed_event_start(&(state->env0_time_end),&(state->env0_time_end_last));	 
            add_fixed_event(0,0,t+*dt+duration_env1,&(state->env1_time_end),&(state->env1_time_end_last));
            break;

        case 4:	/*change env from 1 to 0*/
            state->cumulative_growth_rate*=log(state->cell_size/state->cell_size_copy)/duration_env1;
            state->cell_size_copy=state->cell_size;
            state->env_change++;
            *dt = state->env1_time_end->time - t;
            *env = 0;
            delete_fixed_event_start(&(state->env1_time_end),&(state->env1_time_end_last));	 
            add_fixed_event(0,0,t+*dt+duration_env0,&(state->env0_time_end),&(state->env0_time_end_last));
            break;

        default:
//            LOG_ERROR("event=%d should never get here\n", event);
            exit(-1);
            break;
    }    
}

int do_Gillespie_event(Genotype *genotype,
                        CellState *state,
                        GillespieRates *rates,
                        float dt,
                        float t,
                        long int *seed)
{
    float x;    
    
    x=ran1(seed)*(rates->subtotal);    
    
    if (verbose) log_snapshot(genotype, state, rates, x, t);    
    
    if (x < rates->transport)   /* transportation event */ 
    {   
//        LOG_ERROR("transport event\n");
        transport_event(genotype, state, rates, dt, t, seed);
    } 
    else 
    {
//        LOG_ERROR("x = %f\n", x); 
        x -= rates->transport;
       // LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
        
        if (x < rates->mRNAdecay)  /*STOCHASTIC EVENT: an mRNA decay event */
        { 
//            LOG_ERROR("decay event\n");
            mRNA_decay_event(rates, state, genotype, seed);
                                                    
        } 
        else 
        {
//            LOG_ERROR("x = %f\n", x);            
            x -= rates->mRNAdecay;
           // LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
          
            if (x < rates->pic_disassembly) /* STOCHASTIC EVENT: PIC disassembly*/
            {
//                LOG_ERROR("pic disassembly event\n");
                disassemble_PIC_event(genotype, state, rates, seed);
            } 
            else 
            {
//                LOG_ERROR("x = %f\n", x);
                x -= rates->pic_disassembly;
//                LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );

                if (x < rates->acetylation)  /* acetylation*/
                {
//                    LOG_ERROR("hist act event\n");
                    histone_acteylation_event(rates, state, genotype, seed);
                } 
                else 
                { 
                    x-= rates->acetylation;
                  
//                    LOG_ERROR("CHECK! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
                    
                    if (x < rates->deacetylation)/* STOCHASTIC EVENT: histone deacetylation */ 
                    {
//                        LOG_ERROR("deact event\n"); 
                        histone_deacteylation_event(rates, state, genotype, seed);
                    } 
                    else 
                    {
                        x -= rates->deacetylation; 
                        
                        if (x < rates->pic_assembly)/* STOCHASTIC EVENT: PIC assembly*/
                        {
//                            LOG_ERROR("pic assembly event\n");
                            assemble_PIC_event(rates, state, genotype, seed);                            
                        } 
                        else 
                        {
                            x -= rates->deacetylation;
                            
                            if (x < (float)rates->transcript_init * TRANSCRIPTINIT) /* STOCHASTIC EVENT: transcription initiation */
                            {
//                                LOG_ERROR("transcript init event time = %f\n", t);
                                transcription_init_event(rates, state, genotype, dt, t, seed);
                            } 
                            else 
                            {
                                /*
                                 * FALLBACK: shouldn't get here, previous
                                 * events should be exhaustive
                                 */

//                                LOG_ERROR("[cell %03d] t=%g no event assigned: x=%g, recalibrate cell\n", 
//                                          state->cell_id, t, x);

//                                log_snapshot(rates, state, genotype, x, t);
                                calc_all_rates(genotype, state, rates, UPDATE_ALL); 
//                                log_snapshot(rates, state, genotype, x, t);
                                
                                if (t > 1000.0) 
                                {
//                                    LOG_ERROR("cell has been running for %g which is > 1000 mins\n", t);
                                }
                                /* ignore this timestep, we throw out the dt by returning without updating t */
                                return 0;
                            }
                        }
                    }
                }
            }       
        }
    }
}


void free_fixedevent(CellState *state)
{
    FixedEvent *temp1, *temp2;

    temp1=state->env1_time_end;

    while(temp1){		
            temp1=temp1->next;
            temp2=temp1;
            free(temp2);	
    }

    temp1=state->env0_time_end;
    while(temp1){		
            temp1=temp1->next;
            temp2=temp1;
            free(temp2);	
    }

    temp1=state->mRNA_transcr_time_end;
    while(temp1){		
            temp1=temp1->next;
            temp2=temp1;
            free(temp2);	
    }

    temp1=state->mRNA_transl_time_end;
    while(temp1){		
            temp1=temp1->next;
            temp2=temp1;
            free(temp2);	
    }
}

void calc_avg_growth_rate(Genotype *genotype, 
                           CellState *state, 
                           float init_mRNA[NGENES],
                           float init_protein_conc[NGENES],
                           GillespieRates *rates,
                           float maxbound2,
                           float maxbound3,
                           long int *seed)
{   
    int i,j;
    int env;
    float random;
    float avg_GR1,avg_GR2;  
    float t;
    
    i=0;
    avg_GR1=0.0;
    avg_GR2=0.0;
    
    /* now calc growth rate under two environments*/

    while(i<N_replicates) /* constant env 1*/
    {	 
        env=1;
        for(j=0; j < genotype->ngenes-2; j++) /* loop through tf genes*/
        {   
            random=gasdev(seed)*0.1; // reduce sd to 0.1
            
            init_mRNA[j] = exp(0.91966*random-0.465902);
        
            init_protein_conc[j] = exp(1.25759*random+7.25669);
        }

        init_mRNA[genotype->ngenes-2]=1.0;
        init_mRNA[genotype->ngenes-1]=0.0;
        init_protein_conc[genotype->ngenes-2]=20000.0;
        init_protein_conc[genotype->ngenes-1]=1200.0;
        
        initialize_cell(state, genotype->ngenes, genotype->nproteins,
                        genotype->protein_pool,genotype->mRNAdecay, init_mRNA, init_protein_conc, seed);	
   	    
        calc_all_rates(genotype, state, rates, UPDATE_ALL);	
    	
        t = 0.0;
//printf("%d\n",i);
        while(t<tdevelopment)
        {
            do_single_timestep(genotype, 
                                state,                                            
                                rates, 
                                &t,                                        
                                maxbound2,
                                maxbound3,                                                                                         
                                &env,
                                seed);
        }
        
//        if(state->cell_size>1.2)
//        {
            i++;
            avg_GR1+=log(state->cell_size)/tdevelopment;
//        }  
        
//        iteration++;
        
        free_fixedevent(state);	   			   
    } 
    
    i=0;
    while(i<N_replicates) /* constant env 0*/
    {	 
        env=0;
        for(j=0; j < genotype->ngenes-2; j++) /* loop through tf genes*/
        {   
            random=gasdev(seed)*0.1; // reduce sd to 0.1
            
            init_mRNA[j] = exp(0.91966*random-0.465902);
        
            init_protein_conc[j] = exp(1.25759*random+7.25669);
        }

        init_mRNA[genotype->ngenes-2]=0.0;
        init_mRNA[genotype->ngenes-1]=1.0;
        init_protein_conc[genotype->ngenes-2]=1200.0;
        init_protein_conc[genotype->ngenes-1]=30000.0;
        
        initialize_cell(state, genotype->ngenes, genotype->nproteins,
                        genotype->protein_pool,genotype->mRNAdecay, init_mRNA, init_protein_conc, seed);	
   	    
        calc_all_rates(genotype, state, rates, UPDATE_ALL);	
    	
        t = 0.0;
//printf("%d\n",i);
        while(t<tdevelopment)
        {
            do_single_timestep(genotype, 
                                state,                                            
                                rates, 
                                &t,                                        
                                maxbound2,
                                maxbound3,                                                                                         
                                &env,
                                seed);
        }
        
//        if(state->cell_size>1.2)
//        {
            i++;
            avg_GR2+=log(state->cell_size)/tdevelopment;
//        }  
        
//        iteration++;
        free_fixedevent(state);	   			   
    } 
    genotype->avg_G1=avg_GR1;
    genotype->avg_G2=avg_GR2;
    genotype->fitness=(float)sqrt(avg_GR1*avg_GR2)/N_replicates;   
}

//void try_fixation(Genotype *current_genotype, Genotype *new_genotype, float ori_fitness, int *pfixation, int *OuterWhile, long int *seed)
//{	
//    if(*pfixation) return; // if other threads have reported fixation
//    else
//    {
//        float s, P_fix, ref;
//        
//        s = (new_genotype->fitness-ori_fitness)/ori_fitness;
//
//        if (fabs(s)<EPSILON){P_fix = 1/(float)POP_SIZE;}	
//        else{ P_fix = (1-exp(-s))/(1-exp(-s*(float)POP_SIZE)); }
//
//        ref=ran1(seed);
//
//        if(ref > P_fix) return;
//        
//        else 
//        {
//            int clone_type=5;            
//            clone_cell(new_genotype, current_genotype, clone_type);
//            *pfixation=1;
//            *OuterWhile+=1;
//            return;
//        }
//    }
//}

/* begin of mutation functions*/

void substitution(Genotype *genotype,long int *seed)
{
    int l_genome = genotype->ngenes*CISREG_LEN;
    int pos_n, pos_g;
    float random;
    char n;
    char *Genome;
    Genome= &genotype->cisreg_seq[0][0][0];
    
    random=ran1(seed)*l_genome;					
    pos_n=floor(random);		
    n=set_base_pair(ran1(seed));
    while (n == Genome[pos_n])
    {	
        n=set_base_pair(ran1(seed));
    }	
    Genome[pos_n]=n;   
    
    pos_g=pos_n/CISREG_LEN;
    genotype->re_calc[pos_g][0]=-1; /* we cannot copy tf distribution from elsewhere anymore*/
    genotype->re_calc[pos_g][1]=0;  /* or copy the distribution from this promoter */
    genotype->re_calc[pos_g][2]=1;  /* we need to recalc the binding sites on this promoter*/
    genotype->re_calc[pos_g][3]=1;
}

/* MAXCOPIES has to be 1*/
void insertion(Genotype *genotype,long int *seed)
{
    int inset_size=0;
    float random;
    int pos_g,pos_n,i;
    
    while(inset_size<=0 )
    {
        random=ran1(seed)*max_inset;			
        inset_size = round(random);
    }							
    
    random=ran1(seed)*genotype->ngenes;			
    pos_g=floor(random);   			
    random=ran1(seed)*CISREG_LEN;			
    pos_n=floor(random);                /* at which new seq will be inserted*/
    
    if (pos_n+inset_size>CISREG_LEN) inset_size=CISREG_LEN-pos_n;

    for(i=1;i<=inset_size;i++)
    {
        genotype->cisreg_seq[pos_g][0][CISREG_LEN-i]=genotype->cisreg_seq[pos_g][0][CISREG_LEN-inset_size-i];
    }
    
    for(i=pos_n;i<pos_n+inset_size;i++)
    {        					
        genotype->cisreg_seq[pos_g][0][i]=set_base_pair(ran1(seed));							
    }
    
    genotype->re_calc[pos_g][0]=-1;
    genotype->re_calc[pos_g][1]=0;
    genotype->re_calc[pos_g][2]=1;
    genotype->re_calc[pos_g][3]=1;
}

void partial_deletion(Genotype *genotype, long int *seed)
{
    int delet_size;
    float random;
    int pos_g, pos_n, i;
    
    while(delet_size<=0)
    {
        random=ran1(seed)*max_delet;			
        delet_size = round(random);
    }						
    random=ran1(seed)*genotype->ngenes;			
    pos_g=floor(random);
    random=ran1(seed)*CISREG_LEN;			
    pos_n=floor(random);                /* from which a seq will be deleted */

    if (pos_n+delet_size>CISREG_LEN)	/* if only the tail is deleted*/
    {
        delet_size=CISREG_LEN-pos_n;
        
        for(i=pos_n;i<pos_n+delet_size;i++)
        {					
            genotype->cisreg_seq[pos_g][0][i]=set_base_pair(ran1(seed));
        }				
    }
    else /* else, join the two fragments aside the deletion */
    {
        for(i=pos_n;i<CISREG_LEN-delet_size;i++)
        {
            genotype->cisreg_seq[pos_g][0][i]=genotype->cisreg_seq[pos_g][0][i+delet_size];
        }                				
        for(i++;i<CISREG_LEN;i++) /* and fill the gab by generating new seq */
        {				
            genotype->cisreg_seq[pos_g][0][i]= set_base_pair(ran1(seed));
        }						
    }	
    
    genotype->re_calc[pos_g][0]=-1;
    genotype->re_calc[pos_g][1]=0;
    genotype->re_calc[pos_g][2]=1;
    genotype->re_calc[pos_g][3]=1;
}

void whole_gene_deletion(Genotype *genotype,long int *seed) // any gene can be deleted
{
    float random;
    int pos_g, pos_g_copy, offset, i,j;
    char *temp1;
    float *mRNAdecay, *proteindecay, *translation, *pic_dis;
    int protein_id,protein_id_copy,gene_id,n_gene_of_proteinA, n_gene_of_proteinB, gene_of_proteinA, gene_of_proteinB;
    
    /* check first whether the fitness proteins have extra copies of genes*/
    n_gene_of_proteinA=genotype->protein_pool[genotype->nproteins-2][0][0];
    gene_of_proteinA=genotype->protein_pool[genotype->nproteins-2][1][0];
    n_gene_of_proteinB=genotype->protein_pool[genotype->nproteins-1][0][0];
    gene_of_proteinB=genotype->protein_pool[genotype->nproteins-1][1][0];
    
    if(n_gene_of_proteinA==1 && n_gene_of_proteinB==1)
    {
        random=ran1(seed)*genotype->ntfgenes;			
        pos_g=floor(random);
    }
    else 
    {
        if(n_gene_of_proteinA==1 && n_gene_of_proteinB!=1)
        {
            pos_g=gene_of_proteinA;
            while(pos_g==gene_of_proteinA)
            {
                random=ran1(seed)*genotype->ngenes;			
                pos_g=floor(random);
            }
        }
        else
        {
            if(n_gene_of_proteinA!=1 && n_gene_of_proteinB==1)
            {
                pos_g=gene_of_proteinB;
                while(pos_g==gene_of_proteinB)
                {
                    random=ran1(seed)*genotype->ngenes;			
                    pos_g=floor(random);
                }
            }
            else
            {
                random=ran1(seed)*genotype->ngenes;
                pos_g=floor(random);
            }
        }
    } 
    				
    temp1 = &genotype->cisreg_seq[pos_g][0][0];			
    offset=CISREG_LEN;

    /* remove it from cis_seq */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-pos_g-1);i++) 
    {				
        *temp1=*(temp1+offset);         /* move sequence to fill the gap */
        temp1++;				
    }   
    
    /* if a tf gene is deleted*/
    protein_id=genotype->which_protein[pos_g];
    if(protein_id<genotype->nproteins-2)
    {        
        /* remove from binding seq */
        temp1=&genotype->tf_seq[pos_g][0][0];
        offset=TF_ELEMENT_LEN;

        for(i=0;i<TF_ELEMENT_LEN*(genotype->ntfgenes-pos_g-1);i++)
        {
            *temp1=*(temp1+offset);
            temp1++;
        }

        /* remove from rc binding seq */    
        temp1=&genotype->tf_seq_rc[pos_g][0][0];

        for(i=0;i<TF_ELEMENT_LEN*(genotype->ntfgenes-pos_g-1);i++)
        {
            *temp1=*(temp1+offset);
            temp1++;
        }
        
        /* reduce ntfgenes*/
        genotype->ntfgenes--;
    }
    
    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
    for(i=pos_g;i<genotype->ngenes;i++)
    {
        genotype->pic_disassembly[i][0]=genotype->pic_disassembly[i+1][0];             /* shift elements in the array*/
        genotype->mRNAdecay[i]=genotype->mRNAdecay[i+1];
        genotype->proteindecay[i]=genotype->proteindecay[i+1];
        genotype->translation[i]=genotype->translation[i+1];
        gene_id=genotype->re_calc[i+1][0];
        genotype->re_calc[i][0]=(gene_id>pos_g)?gene_id-1:gene_id; /*note that deletion changes the ids of genes!!!*/
        genotype->re_calc[i][1]=genotype->re_calc[i+1][1];
        genotype->re_calc[i][2]=1;              /* need to recalc BS, because of the mismatch between gene id and all_binding_sites*/
        genotype->re_calc[i][3]=1;              /* copy info back to this site in clone_cell */
        
        /* now move the info about binding sites*/
        for(j=0;j<genotype->binding_sites_num[i+1];j++)
        {
            genotype->all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i+1][j].tf_id;
            genotype->all_binding_sites[i][j].Koff=genotype->all_binding_sites[i+1][j].Koff;
            genotype->all_binding_sites[i][j].BS_pos=genotype->all_binding_sites[i+1][j].BS_pos;
            genotype->all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i+1][j].N_hindered;
        }
        genotype->binding_sites_num[i]=genotype->binding_sites_num[i+1];
        genotype->max_hindered_sites[i]=genotype->max_hindered_sites[i+1];
        genotype->N_act_BS[i]=genotype->N_act_BS[i+1];
        genotype->N_rep_BS[i]=genotype->N_rep_BS[i+1];
        
//        for(j=0;j<genotype->max_N_rep_bound[i+1];j++)
//        {
//            genotype->N_configurations[i][j]=genotype->N_configurations[i+1][j];
//        }        
//        genotype->max_N_rep_bound[i]=genotype->max_N_rep_bound[i+1];
//        genotype->max_N_act_bound[i]=genotype->max_N_act_bound[i+1];
    }
    
    /* now change protein_pool*/
    if(genotype->protein_pool[protein_id][0][0]==1) /* if this is the only gene copy */
    {    
        protein_id_copy=protein_id;
        for(i=0;i<genotype->nproteins-protein_id;i++)  /* then we need to remove this protein from protein_pool*/
        {            
            genotype->protein_pool[protein_id_copy][0][0]=genotype->protein_pool[protein_id_copy+1][0][0];
            for(j=0;j<genotype->protein_pool[protein_id_copy][0][0];j++)
            {
                gene_id=genotype->protein_pool[protein_id_copy+1][1][j];
                genotype->protein_pool[protein_id_copy][1][j]=(gene_id>pos_g)?gene_id-1:gene_id;/*note that deletion changes the ids of genes!!!*/
            }            
            protein_id_copy++;
        }
        
        genotype->protein_pool[protein_id_copy-1][1][0]=0; /*added here for stability*/
         
        if(protein_id<genotype->nproteins-2) /*if a tf is deleted*/
        {    
            if(genotype->activating[protein_id][0] && (genotype->N_act!=0)) /* reduce the number of activator if necessary */
                genotype->N_act--;
            if(!genotype->activating[protein_id][0] && (genotype->N_rep!=0))
                genotype->N_rep--;

            protein_id_copy=protein_id;
            for(i=0;i<genotype->nproteins-protein_id-2;i++) /* also remove it from activating */
            {
                genotype->activating[protein_id_copy][0]=genotype->activating[protein_id_copy+1][0];
                protein_id_copy++;
            }
            
            /* in the case, all genes need to recalc binding sites*/
            for(i=0;i<pos_g;i++)
            {
                genotype->re_calc[i][2]=1; /* recalc BS */
                genotype->re_calc[i][3]=1; /* copy back the original BS in clone_cell*/
            }
        }
        
        genotype->nproteins--;                /* reduce the number of protein*/
        
        
    }  
    else /*if the protein has more than one genes, remove the one indicated by pos_g*/
    {
        i=0;
        while(genotype->protein_pool[protein_id][1][i]!=pos_g) i++; /* find where is this pos_g */
        for(;i<genotype->protein_pool[protein_id][0][0];i++)
        {
            gene_id=genotype->protein_pool[protein_id_copy][1][i+1];
            genotype->protein_pool[protein_id][1][i]=(gene_id>pos_g)?gene_id-1:gene_id; /*note that deletion changes the ids of genes!!!*/
        }
        genotype->protein_pool[protein_id][0][0]--;
    }
    
    /* change which_protein*/
    pos_g_copy=pos_g;
    for(i=0;i<genotype->ngenes-pos_g;i++)
    {
        protein_id_copy=genotype->which_protein[pos_g_copy+1];
        genotype->which_protein[pos_g_copy]=(protein_id_copy>protein_id)?protein_id_copy-1:protein_id_copy; /*the deletion also changes the ids of proteins*/
        pos_g_copy++;
    } 
        
    /* change ngenes*/
    genotype->ngenes--;
}

//void whole_gene_deletion(Genotype *genotype,long int *seed) // only tf genes get deleted. Haven't been debugged
//{
//    float random;
//    int pos_g, offset, i,j;
//    char *temp1;
//    float *mRNAdecay, *proteindecay, *translation, *pic_dis;
//    int protein_id;
//    
//    random=ran1(seed)*genotype->ntfgenes;			
//    pos_g=floor(random);				
//    temp1 = &genotype->cisreg_seq[pos_g][0][0];			
//    offset=CISREG_LEN;
//
//    /* remove it from cis_seq */
//    for(i=0;i<CISREG_LEN;i++) 
//    {				
//        *temp1=*(temp1+offset);         /* move sequence to fill the gap */
//        temp1++;				
//    }   
//    
//    /* remove from binding seq */
//    temp1=&genotype->tf_seq[pos_g][0][0];
//    offset=TF_ELEMENT_LEN;
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp1=*(temp1+offset);
//        temp1++;
//    }
//    
//    /* remove from rc binding seq */    
//    temp1=&genotype->tf_seq_rc[pos_g][0][0];
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp1=*(temp1+offset);
//        temp1++;
//    }
//    
//    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
//    for(i=pos_g;i<genotype->ngenes;i++)
//    {
//        genotype->pic_disassembly[i][0]=genotype->pic_disassembly[i+1][0];             /* shift elements in the array*/
//        genotype->mRNAdecay[i]=genotype->mRNAdecay[i+1];
//        genotype->proteindecay[i]=genotype->proteindecay[i+1];
//        genotype->translation[i]=genotype->translation[i+1];
//        genotype->re_calc[i][0]=genotype->re_calc[i+1][0];
//        genotype->re_calc[i][1]=genotype->re_calc[i+1][1];
//        genotype->re_calc[i][2]=genotype->re_calc[i+1][2];
//        genotype->re_calc[i][3]=1;              /* copy info back to this site in clone_cell */
//        
//        /* now move the info about binding sites*/
//        for(j=0;j<genotype->binding_sites_num[i+1];j++)
//        {
//            genotype->all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i+1][j].tf_id;
//            genotype->all_binding_sites[i][j].Koff=genotype->all_binding_sites[i+1][j].Koff;
//            genotype->all_binding_sites[i][j].BS_pos=genotype->all_binding_sites[i+1][j].BS_pos;
//            genotype->all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i+1][j].N_hindered;
//        }
//        genotype->binding_sites_num[i]=genotype->binding_sites_num[i+1];
//        genotype->max_hindered_sites[i]=genotype->max_hindered_sites[i+1];
//        genotype->N_act_BS[i]=genotype->N_act_BS[i+1];
//        genotype->N_rep_BS[i]=genotype->N_rep_BS[i+1];
//        
////        for(j=0;j<genotype->max_N_rep_bound[i+1];j++)
////        {
////            genotype->N_configurations[i][j]=genotype->N_configurations[i+1][j];
////        }        
////        genotype->max_N_rep_bound[i]=genotype->max_N_rep_bound[i+1];
////        genotype->max_N_act_bound[i]=genotype->max_N_act_bound[i+1];
//    }   
//    
//    /* change ngenes and ntfgenes*/
//    genotype->ngenes--;
//    genotype->ntfgenes--;
//    
//    /* now change protein_pool*/
//    protein_id=genotype->which_protein[pos_g];
//    
//    if(genotype->protein_pool[protein_id][0][0]==1) /* if this is the only gene copy */
//    {    
//        for(i=0;i<genotype->nproteins-protein_id;i++)  /* then we need to remove this protein from protein_pool*/
//        {
//            genotype->protein_pool[protein_id][0][0]=genotype->protein_pool[protein_id+1][0][0];
//            for(j=0;j<genotype->protein_pool[protein_id][0][0];j++)
//            {
//                genotype->protein_pool[protein_id][1][j]=genotype->protein_pool[protein_id+1][1][j];
//            }
//            protein_id++;
//        }
//        
//        if(genotype->activating[protein_id][0] && (genotype->N_act!=0)) /* reduce the number of activator if necessary */
//            genotype->N_act--;
//        if(!genotype->activating[protein_id][0] && (genotype->N_rep!=0))
//            genotype->N_rep--;
//                   
//        for(i=0;i<genotype->nproteins-protein_id-2;i++) /* also remove it from activating */
//        {
//            genotype->activating[protein_id][0]=genotype->activating[protein_id+1][0];
//        }
//        
//        genotype->nproteins--;                /* reduce the number of protein*/
//    }  
//    else
//    {
//        i=0;
//        while(genotype->protein_pool[protein_id][1][i]!=pos_g)i++; /* find where is this pos_g */
//        for(;i<genotype->protein_pool[protein_id][0][0];i++)
//        {
//            genotype->protein_pool[protein_id][1][i]=genotype->protein_pool[protein_id][1][i+1];
//        }
//        genotype->protein_pool[protein_id][0][0]--;
//    }
//}

//void gene_duplication(Genotype *genotype,long int *seed) //only tf genes can be duplicated
//{
//    float random;
//    int pos_g, i,j, protein_id;
//    char *temp1, *temp2;
//    
//    random=ran1(seed)*genotype->ntfgenes;
//    pos_g=floor(random);
//    
//    /* copy the promoter*/
//    temp1=&genotype->cisreg_seq[pos_g][0][0]; /* points to the gene to be duplicated*/
//    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0][CISREG_LEN-1]; /* points to the end of the 2nd selection gene */
//    
//    for(i=0;i<2*CISREG_LEN;i++) 
//    {
//        *(temp2+CISREG_LEN)=*temp2; /* shift the sequences of the selection genes CISREG_LEN bp */
//        temp2--;
//    }
//    
//    temp2=&genotype->cisreg_seq[genotype->ngenes-2][0][0]; /* point temp2 to the start of the original first selection gene */
//    
//    for(i=0;i<CISREG_LEN;i++) 
//    {
//        *temp2++=*temp1++;  /* put the duplicated gene at the original place of the first selection gene */
//    }
//    
//    /* copy the binding sequence*/
//    temp1=&genotype->tf_seq[pos_g][0][0];
//    temp2=&genotype->tf_seq[genotype->ntfgenes][0][0];   
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp2++=*temp1++;
//    }
//    
//    /* copy the rc binding sequence*/
//    temp1=&genotype->tf_seq_rc[pos_g][0][0];
//    temp2=&genotype->tf_seq_rc[genotype->ntfgenes][0][0];    
//    
//    for(i=0;i<TF_ELEMENT_LEN;i++)
//    {
//        *temp2++=*temp1++;
//    }  
//    
//    /* add it to protein_pool, but do not change nproteins*/    
//    protein_id=genotype->which_protein[pos_g];
//    
//    genotype->protein_pool[protein_id][1][genotype->protein_pool[protein_id][0][0]]=genotype->ntfgenes;
//    
//    genotype->protein_pool[protein_id][0][0]++; 
//    
//    /* add it to PIC_assembly, mRNAdecay, proteinDecay, and translation*/
//    for(i=genotype->ngenes;i>genotype->ntfgenes;i--)
//    {
//        genotype->pic_disassembly[i][0]=genotype->pic_disassembly[i-1][0]; /* shift the selection genes to make a slot*/
//        genotype->mRNAdecay[i]=genotype->mRNAdecay[i-1];
//        genotype->proteindecay[i]=genotype->proteindecay[i-1];
//        genotype->translation[i]=genotype->translation[i-1];
//        genotype->re_calc[i][0]=genotype->re_calc[i-1][0];
//        genotype->re_calc[i][1]=genotype->re_calc[i-1][1];
//        genotype->re_calc[i][2]=genotype->re_calc[i-1][2];
//        genotype->re_calc[i][3]=1;        /* copy info back to this site in clone_cell */                   
//        
//         /* now move the info about binding sites*/
//        for(j=0;j<genotype->binding_sites_num[i-1];j++)
//        {
//            genotype->all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i-1][j].tf_id;
//            genotype->all_binding_sites[i][j].Koff=genotype->all_binding_sites[i-1][j].Koff;
//            genotype->all_binding_sites[i][j].BS_pos=genotype->all_binding_sites[i-1][j].BS_pos;
//            genotype->all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i-1][j].N_hindered;
//        }
//        genotype->binding_sites_num[i]=genotype->binding_sites_num[i-1];
//        genotype->max_hindered_sites[i]=genotype->max_hindered_sites[i-1];
//        genotype->N_act_BS[i]=genotype->N_act_BS[i-1];
//        genotype->N_rep_BS[i]=genotype->N_rep_BS[i-1];
//        
////        for(j=0;j<genotype->max_N_rep_bound[i-1];j++)
////        {
////            genotype->N_configurations[i][j]=genotype->N_configurations[i-1][j];
////        }        
////        genotype->max_N_rep_bound[i]=genotype->max_N_rep_bound[i-1];
////        genotype->max_N_act_bound[i]=genotype->max_N_act_bound[i-1];
//    }
//    
//    genotype->pic_disassembly[genotype->ntfgenes][0]=genotype->pic_disassembly[pos_g][0];
//    genotype->mRNAdecay[genotype->ntfgenes]=genotype->mRNAdecay[pos_g];
//    genotype->proteindecay[genotype->ntfgenes]=genotype->proteindecay[pos_g];
//    genotype->translation[genotype->ntfgenes]=genotype->translation[pos_g];
//    genotype->re_calc[genotype->ntfgenes][0]=pos_g;        /* we can copy the tf distribution from pos_g */
//    genotype->re_calc[genotype->ntfgenes][1]=1;            /* we also copy the tf distribution from here */
//    genotype->re_calc[genotype->ntfgenes][2]=0;            /* we do not need to calculate the binding sites */
//    genotype->re_calc[pos_g][1]=1;                         /* we just duplicated pos_g, so we can copy the tf distribution from pos_g */
//    genotype->re_calc[genotype->ntfgenes][3]=1;            /* copy info back to this site in clone_cell */
//    
//    /* now copy the info about binding sites*/
//    for(j=0;j<genotype->binding_sites_num[pos_g];j++)
//        {
//            genotype->all_binding_sites[genotype->ntfgenes][j].tf_id=genotype->all_binding_sites[pos_g][j].tf_id;
//            genotype->all_binding_sites[genotype->ntfgenes][j].Koff=genotype->all_binding_sites[pos_g][j].Koff;
//            genotype->all_binding_sites[genotype->ntfgenes][j].BS_pos=genotype->all_binding_sites[pos_g][j].BS_pos;
//            genotype->all_binding_sites[genotype->ntfgenes][j].N_hindered=genotype->all_binding_sites[pos_g][j].N_hindered;
//        }
//        genotype->binding_sites_num[genotype->ntfgenes]=genotype->binding_sites_num[pos_g];
//        genotype->max_hindered_sites[genotype->ntfgenes]=genotype->max_hindered_sites[pos_g];
//        genotype->N_act_BS[genotype->ntfgenes]=genotype->N_act_BS[pos_g];
//        genotype->N_rep_BS[genotype->ntfgenes]=genotype->N_rep_BS[pos_g];
//        
////        for(j=0;j<genotype->max_N_rep_bound[pos_g];j++)
////        {
////            genotype->N_configurations[genotype->ntfgenes][j]=genotype->N_configurations[pos_g][j];
////        }        
////        genotype->max_N_rep_bound[genotype->ntfgenes]=genotype->max_N_rep_bound[pos_g];
////        genotype->max_N_act_bound[genotype->ntfgenes]=genotype->max_N_act_bound[pos_g];
//    
//    /* update gene numbers*/    
//    genotype->ntfgenes++;
//    genotype->ngenes++;
//}

void gene_duplication(Genotype *genotype,long int *seed) //any gene can be duplicated
{
    float random;
    int pos_g,pos_g_copy, i,j, protein_id;
    char *temp1, *temp2;
    
    
    random=ran1(seed)*genotype->ngenes;
    pos_g=floor(random);
    
    if(pos_g>=genotype->ngenes-2)
        pos_g_copy=pos_g+1; /* note that if the selection genes are to be duplicated, shifting sequences and info will cause problem*/
    else
        pos_g_copy=pos_g;
    
    /* copy the promoter*/
    temp1=&genotype->cisreg_seq[pos_g_copy][0][0]; /* points to the gene to be duplicated*/
    temp2=&genotype->cisreg_seq[genotype->ngenes-1][0][CISREG_LEN-1]; /* points to the end of the 2nd selection gene */
    
    for(i=0;i<2*CISREG_LEN;i++) 
    {
        *(temp2+CISREG_LEN)=*temp2; /* shift the sequences of the selection genes CISREG_LEN bp */
        temp2--;
    }
    
    temp2=&genotype->cisreg_seq[genotype->ngenes-2][0][0]; /* point temp2 to the start of the original first selection gene */
    
    for(i=0;i<CISREG_LEN;i++) 
    {
        *temp2++=*temp1++;  /* put the duplicated gene at the original place of the first selection gene */
    }
    
    /*if a tf gene is duplicated*/
    protein_id=genotype->which_protein[pos_g];
    
    if(protein_id<genotype->nproteins-2) /*note that the fitness proteins are always the last two protein*/
    {
        /* copy the binding sequence*/
        temp1=&genotype->tf_seq[pos_g][0][0];
        temp2=&genotype->tf_seq[genotype->ntfgenes][0][0];   

        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp2++=*temp1++;
        }

        /* copy the rc binding sequence*/
        temp1=&genotype->tf_seq_rc[pos_g][0][0];
        temp2=&genotype->tf_seq_rc[genotype->ntfgenes][0][0];    

        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp2++=*temp1++;
        }  
    }   
    
    /* add it to PIC_assembly, mRNAdecay, proteinDecay, translation, and which_protein*/
    for(i=genotype->ngenes;i>genotype->ngenes-2;i--)
    {
        genotype->pic_disassembly[i][0]=genotype->pic_disassembly[i-1][0]; /* shift the selection genes to make a slot*/
        genotype->mRNAdecay[i]=genotype->mRNAdecay[i-1];
        genotype->proteindecay[i]=genotype->proteindecay[i-1];
        genotype->translation[i]=genotype->translation[i-1];
        genotype->which_protein[i]=genotype->which_protein[i-1];
        genotype->re_calc[i][0]=genotype->re_calc[i-1][0];
        genotype->re_calc[i][1]=genotype->re_calc[i-1][1];
        genotype->re_calc[i][2]=genotype->re_calc[i-1][2];
        genotype->re_calc[i][3]=1;        /* copy info back to this site in clone_cell */                   
        
         /* now move the info about binding sites*/
        for(j=0;j<genotype->binding_sites_num[i-1];j++)
        {
            genotype->all_binding_sites[i][j].tf_id=genotype->all_binding_sites[i-1][j].tf_id;
            genotype->all_binding_sites[i][j].Koff=genotype->all_binding_sites[i-1][j].Koff;
            genotype->all_binding_sites[i][j].BS_pos=genotype->all_binding_sites[i-1][j].BS_pos;
            genotype->all_binding_sites[i][j].N_hindered=genotype->all_binding_sites[i-1][j].N_hindered;
        }
        genotype->binding_sites_num[i]=genotype->binding_sites_num[i-1];
        genotype->max_hindered_sites[i]=genotype->max_hindered_sites[i-1];
        genotype->N_act_BS[i]=genotype->N_act_BS[i-1];
        genotype->N_rep_BS[i]=genotype->N_rep_BS[i-1];
        
//        for(j=0;j<genotype->max_N_rep_bound[i-1];j++)
//        {
//            genotype->N_configurations[i][j]=genotype->N_configurations[i-1][j];
//        }        
//        genotype->max_N_rep_bound[i]=genotype->max_N_rep_bound[i-1];
//        genotype->max_N_act_bound[i]=genotype->max_N_act_bound[i-1];
    }
    
    genotype->pic_disassembly[genotype->ngenes-2][0]=genotype->pic_disassembly[pos_g_copy][0];
    genotype->mRNAdecay[genotype->ngenes-2]=genotype->mRNAdecay[pos_g_copy];
    genotype->proteindecay[genotype->ngenes-2]=genotype->proteindecay[pos_g_copy];
    genotype->translation[genotype->ngenes-2]=genotype->translation[pos_g_copy];
    genotype->which_protein[genotype->ngenes-2]=protein_id;
    
    /* things are a little different if the original selection genes are duplicated, because of their locations*/
    if(pos_g<genotype->ngenes-2)
    {       
        genotype->re_calc[genotype->ngenes-2][0]=pos_g;        /* we can copy the tf distribution from pos_g */
        genotype->re_calc[genotype->ngenes-2][1]=1;            /* we can also copy the tf distribution from here */
        genotype->re_calc[genotype->ngenes-2][2]=0;            /* we do not need to calculate the binding sites */
        genotype->re_calc[pos_g][1]=1;                         /* we just duplicated pos_g, so we can copy the tf distribution from pos_g */
        genotype->re_calc[genotype->ngenes-2][3]=1;            /* copy info back to this site in clone_cell */
    }
    else
    {         
        genotype->re_calc[genotype->ngenes-2][0]=genotype->re_calc[pos_g+1][0]; /* this extra copy can only copy from other copies of the original selection gene */
        genotype->re_calc[genotype->ngenes-2][1]=1;
        genotype->re_calc[genotype->ngenes-2][2]=0; 
        genotype->re_calc[genotype->ngenes-2][3]=1; 
        genotype->re_calc[pos_g+1][0]=genotype->ngenes-2; /* the original selection gene that is duplicated can copy info from here*/
        genotype->re_calc[pos_g+1][1]=1;
        genotype->re_calc[pos_g+1][2]=0;
        genotype->re_calc[pos_g+1][3]=1;
    }
    
    /* now copy the info about binding sites*/
    for(j=0;j<genotype->binding_sites_num[pos_g_copy];j++)
    {
        genotype->all_binding_sites[genotype->ngenes-2][j].tf_id=genotype->all_binding_sites[pos_g_copy][j].tf_id;
        genotype->all_binding_sites[genotype->ngenes-2][j].Koff=genotype->all_binding_sites[pos_g_copy][j].Koff;
        genotype->all_binding_sites[genotype->ngenes-2][j].BS_pos=genotype->all_binding_sites[pos_g_copy][j].BS_pos;
        genotype->all_binding_sites[genotype->ngenes-2][j].N_hindered=genotype->all_binding_sites[pos_g_copy][j].N_hindered;
    }
    genotype->binding_sites_num[genotype->ngenes-2]=genotype->binding_sites_num[pos_g_copy];
    genotype->max_hindered_sites[genotype->ngenes-2]=genotype->max_hindered_sites[pos_g_copy];
    genotype->N_act_BS[genotype->ngenes-2]=genotype->N_act_BS[pos_g_copy];
    genotype->N_rep_BS[genotype->ngenes-2]=genotype->N_rep_BS[pos_g_copy];
        
//        for(j=0;j<genotype->max_N_rep_bound[pos_g];j++)
//        {
//            genotype->N_configurations[genotype->ntfgenes][j]=genotype->N_configurations[pos_g][j];
//        }        
//        genotype->max_N_rep_bound[genotype->ntfgenes]=genotype->max_N_rep_bound[pos_g];
//        genotype->max_N_act_bound[genotype->ntfgenes]=genotype->max_N_act_bound[pos_g];
    
    /* update protein_pool*/
    /* first add it to protein_pool, but do not change nproteins*/    
    genotype->protein_pool[protein_id][1][genotype->protein_pool[protein_id][0][0]]=genotype->ngenes-2; /*the newly duplicated gene takes the original place of the first selection gene*/
    genotype->protein_pool[protein_id][0][0]++; 
    
    for(i=2;i>0;i--) /*update the id of the original selection genes stored in protein_pool*/
    {
        j=0;
        while(genotype->protein_pool[genotype->nproteins-i][1][j]!=genotype->ngenes-i)j++;
        genotype->protein_pool[genotype->nproteins-i][1][j]++;
    }
    
    /* update gene numbers*/ 
    genotype->ngenes++;    
    if(protein_id<genotype->nproteins-2)
        genotype->ntfgenes++;    
}

void mut_binding_sequence(Genotype *genotype,long int *seed)
{
    float random;
    int pos_g, pos_n, protein_id, i;
    char n;  
    
    char *tf_seq, *tf_seq_rc;    
    tf_seq=&genotype->tf_seq[0][0][0];
    tf_seq_rc=&genotype->tf_seq_rc[0][0][0];
    
    random=ran1(seed)*TF_ELEMENT_LEN;			
    pos_n=floor(random);
    n=set_base_pair(ran1(seed));
    
    while (n == tf_seq[pos_n])
    {	
        n=set_base_pair(ran1(seed));
    }	
    
    tf_seq[pos_n]=n;
    
    /* update the complement sequence*/
    switch (n)
    {
        case 'g':
            tf_seq_rc[pos_n]='c'; break;
        case 'c':
            tf_seq_rc[pos_n]='g'; break;
        case 'a':
            tf_seq_rc[pos_n]='t'; break;
        case 't':
            tf_seq_rc[pos_n]='a'; break;
    }
    
    /* if this tf gene has more than one copies, the mutation inreases nproteins*/
    pos_g=pos_n/TF_ELEMENT_LEN;
    
    protein_id=genotype->which_protein[pos_g];   
   
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        /* remove this copy of gene for the original protein*/
        i=0;

        while(genotype->protein_pool[protein_id][1][i]!=pos_g) i++;

        for(;i<genotype->protein_pool[protein_id][0][0];i++) 
        {
            genotype->protein_pool[protein_id][1][i]= genotype->protein_pool[protein_id][1][i+1]; /* rearrange data array */
        }

        genotype->protein_pool[protein_id][0][0]--; 
        
        /* create a new protein and link it to this gene*/
        genotype->which_protein[pos_g]=genotype->nproteins-2; /*put the new protein to the pos of the first selection gene*/
        
        genotype->protein_pool[genotype->nproteins-2][0][0]=1;
        
        genotype->protein_pool[genotype->nproteins-2][1][0]=pos_g;
        
        /* update acitivating*/
        if(genotype->activating[protein_id][0]) /* increase the number of activator */
            genotype->N_act++;
        else
            genotype->N_rep++;
        
        genotype->activating[genotype->nproteins-2][0]=genotype->activating[protein_id][0];
        
        /* finally, update protein numbers*/
        genotype->nproteins++;
        
        /* NOTE: this mutation does not change the number of genes*/
    } 
    
    for(i=0;i<genotype->ngenes;i++) 
    {
        genotype->re_calc[i][2]=1;   /* recalculate binding sites on every promoter */
        genotype->re_calc[i][3]=1;   /* copy info back to every gene in clone_cell */
    }
}

/* For the moment, only mRNA_decay, translation, protein_decay, and pic_disassembly 
 * will be mutated. We assume a mutation attacks the four constants with equal 
 * probability. 
 */
void mut_kinetic_constant(Genotype *genotype, float kdis[NUM_K_DISASSEMBLY],long int *seed)
{
    float random1, random2;
    int pos_kdis, pos_g, protein_id, i;
    
    random1=ran1(seed);
    
    random2=ran1(seed)*genotype->ngenes;
        
    pos_g=floor(random2); /* which gene */
    
    if(random1<=0.25) /* mut kdis */
    {        
        random2=ran1(seed)*NUM_K_DISASSEMBLY;
        
        pos_kdis=floor(random2);
        
        while(genotype->pic_disassembly[pos_g][0]==kdis[pos_kdis]) /* be sure to choose a different value*/
        {
            random2=ran1(seed)*NUM_K_DISASSEMBLY;
        
            pos_kdis=floor(random2);
        }
        
        genotype->pic_disassembly[pos_g][0]=kdis[pos_kdis];        
    }
    else if(random1<=0.5) /* mut mRNAdecay */
    {
        random2 = exp(0.4909*gasdev(seed)-3.20304);
            
        while(genotype->mRNAdecay[pos_g]==random2) /* be sure to choose a different value*/
        {
            random2 = exp(0.4909*gasdev(seed)-3.20304);
        }

        genotype->mRNAdecay[pos_g]=random2;
    }
    else if(random1<=0.75) /* mut translation */
    {
        random2= exp(0.7406*gasdev(seed)+4.56);
        
         while(genotype->translation[pos_g]==random2) /* be sure to choose a different value*/
        {
            random2= exp(0.7406*gasdev(seed)+4.56);
        }
        
        protein_id=genotype->which_protein[pos_g];
        
        if(protein_id>=genotype->nproteins-2) /*if this is a selection gene, we mutate its copies as well,*/
        {                                     /* because adding a new selection protein is complicate*/
            for(i=0;i<genotype->protein_pool[protein_id][0][0];i++) 
            {
                genotype->translation[genotype->protein_pool[protein_id][1][i]]=random2;
            }
        }
        else
        {
            genotype->translation[pos_g]=random2;
        }
    }
    else /* mut protein decay */
    {
        random2=genotype->proteindecay[pos_g];
        
        while(random2==genotype->proteindecay[pos_g]) /*make sure getting a different value*/
        {
            random2=-1.0;
            while (random2 < 0.0) 
            {
                if (ran1(seed) < 0.08421)
                    random2 = (float)EPSILON;
                else random2 = exp(0.7874*gasdev(seed)-3.7665);
            }
        }
        
        protein_id=genotype->which_protein[pos_g];
        
        if(protein_id>=genotype->nproteins-2) /*if this is a selection gene, we mutate its copies as well,*/
        {                                     /* because adding a new selection protein is complicate*/
            for(i=0;i<genotype->protein_pool[protein_id][0][0];i++) 
            {
                genotype->proteindecay[genotype->protein_pool[protein_id][1][i]]=random2;
            }
        }
        else
        {
            genotype->proteindecay[pos_g]=random2;
        }
        /* Theoretically,if this protein has more than one copy of genes, we need to make a new protein.
         * However, in order to reduce to burden of computing tf binding distribution, we treat the mutant as the 
         * original tf. We combine the concentration of the mutant and the original to compute binding distribution*/         
    }    
}

int mutate(Genotype *genotype, float kdis[NUM_K_DISASSEMBLY],long int *seed, char *mut_type2)
{
    char mut_type;
 
    draw_mutation(genotype->ngenes,&mut_type,seed);
    
//    mut_type='d';
    *mut_type2=mut_type;
    
//    mut_type='k';
    
    switch (mut_type)
    {
        case 's': //substitution        		
            substitution(genotype,seed);			
            return 1;
            break;
        		
        case 'i': // insertion        
            insertion(genotype,seed);	
            return 1;
            break;	
        		
        case 'p': // partial deletion        
            partial_deletion(genotype,seed);
            return 1;
            break;			
        		
        case 'w': // whole gene deletion. This mutation has two versions: only tf genes get deleted or any
                  // gene can be deleted            
            whole_gene_deletion(genotype,seed);
            return 2;
            break;
        		
        case 'd': // Whole gene duplication also has two versions                  
            gene_duplication(genotype,seed);
            return 2;
            break;
        
        case 'c': //binding sequence        
            mut_binding_sequence(genotype,seed);
            return 3;
            break;  
            
        case 'k': //mutations in kinetic constants        
            mut_kinetic_constant(genotype, kdis,seed);
            return 4;
            break;        
    }
}

/* this function calculates the probability of different mutations based on
 * the current genotype. It then modify the value of mut_type
 */
void draw_mutation(int ngenes, char *mut_type, long int *seed)
{
    float random,random2;
    float tot_mut_rate=0.0;
    float tot_subs_rate, tot_indel_rate, tot_dup_rate, tot_sil_rate, tot_kin_rate;   
    
    /* calc total susbtitution rate*/
    tot_subs_rate=(float)ngenes*CISREG_LEN*SUBSTITUTION;
    tot_mut_rate+=tot_subs_rate;
    
    /* indel rate*/
    tot_indel_rate=(float)ngenes*CISREG_LEN*INDEL;
    tot_mut_rate+=tot_indel_rate;
    
    /* duplication rate*/
    tot_dup_rate=(float)(ngenes-2)*DUPLICATION;
    tot_mut_rate+=tot_dup_rate;
    
    /* silencing rate*/
    tot_sil_rate=(float)(ngenes-2)*SILENCING;
    tot_mut_rate+=tot_sil_rate;
    
    /* mut in kinetic constants and binding seq */
    tot_kin_rate=(tot_subs_rate+tot_indel_rate)*MUTKINETIC;
    tot_mut_rate+=tot_kin_rate;
    
    random=ran1(seed);
    
    if(random<=tot_kin_rate/tot_mut_rate)
    {
        random2=ran1(seed);
        
        if(random2<0.5)
            *mut_type='c';                              /* mut binding seq*/
        else
            *mut_type='k';                              /* mut kinetic const*/
    }
    else
    {
        random-=tot_kin_rate/tot_mut_rate;
        
        if(random<=tot_sil_rate/tot_mut_rate)
            *mut_type='w';                              /* while gene deletion*/        
        else
        {
            random -= tot_sil_rate/tot_mut_rate;
            
            if(random<= tot_dup_rate/tot_mut_rate)
                *mut_type='d';                          /* gene duplication */
            else
            {
                random-=tot_dup_rate/tot_mut_rate;
                
                if(random<=tot_indel_rate/tot_mut_rate)
                    *mut_type='i';                      /* indel*/
                else
                    *mut_type='s';                      /* substituion*/
            }
                
        }
    }
    
}
/* end of mutation functions*/

void initialize_cache(Genotype *genotype)
{
    int j,k;    
    
    for(j=0;j<NGENES;j++)
    {
        genotype->re_calc[j][0]=-1; /* we cannot copy distribution for this gene from elsewhere*/
        genotype->re_calc[j][1]=0;  /* we cannot copy distribution from here*/
        genotype->re_calc[j][2]=1;  /* we need to calc binding sites for this gene*/
        genotype->re_calc[j][3]=1;  /* we need to copy info for this gene in clone_cell */
    }

    /* alloc space or protein_pool */
    for(j=0;j<NPROTEINS;j++)
    {
        genotype->protein_pool[j][0]=malloc(sizeof(int));
        genotype->protein_pool[j][1]=malloc(MAXALLOC*sizeof(int));

        genotype->protein_pool[j][0][0]=0;
        for(k=0;k<MAXALLOC;k++)
        {
            genotype->protein_pool[j][1][k]=-1;
        }
    }
    
    for(j=0;j<NGENES;j++)
    {
        genotype->all_binding_sites[j] = malloc(MAXELEMENTS*sizeof(AllTFBindingSites)); 
//        genotype->N_configurations[j] = malloc(MAXELEMENTS*sizeof(int));
        
        if (!(genotype->all_binding_sites[j])) 
        {
    //        LOG_ERROR_NOCELLID("initial setting of all_binding_sites failed.\n");
            exit(1);
        }
    } 
}


void init_run_pop(//Genotype genotype[N_para_threads+1],
                  //CellState state[N_para_threads+1],
//                  float temperature,   /* in Kelvin */
                  float kdis[NUM_K_DISASSEMBLY],
                  FILE *OUTPUT)
//                  int no_fixed_dev_time)
{  
    int i;
    int fixation = 0; 
    int maxbound2, maxbound3; 
    Genotype genotype_ori;
    CellState state_ori;
    float init_mRNA[NGENES]; 
    float init_protein_conc[NGENES];    
    GillespieRates rates_ori;	
    maxbound2 = MAXBOUND;
    maxbound3 = 10*MAXBOUND;
    float avg_G1;
    float avg_G2;
    double t1,t2;
    
    

    for(i=0;i<TF_ELEMENT_LEN-NMIN+1;i++)
    {      
        Koff[i]=NUMSITESINGENOME*kon*0.25*KR/exp(-((float)i/3.0-1.0));      
    }     

    initialize_cache(&genotype_ori);
    
    initialize_genotype(&genotype_ori, kdis); 
    
    calc_avg_growth_rate(&genotype_ori,
                            &state_ori,
                            init_mRNA,
                            init_protein_conc,                                                                
                            &rates_ori,
                            maxbound2,
                            maxbound3,
                            &master_seed);        
    
    avg_G1=genotype_ori.avg_G1;
    avg_G2=genotype_ori.avg_G2;
    
    calc_avg_growth_rate(&genotype_ori,
                            &state_ori,
                            init_mRNA,
                            init_protein_conc,                                                                
                            &rates_ori,
                            maxbound2,
                            maxbound3,
                            &master_seed); 
    
    genotype_ori.fitness=(float)sqrt((genotype_ori.avg_G1+avg_G1)*(genotype_ori.avg_G2+avg_G2))/(2*N_replicates);

//    printf("rep %d, fitness=%f\n",N_replicates,genotype[current_genotype].fitness);
    
    i=0;
    
//    omp_set_num_threads(1);
    t1=omp_get_wtime();
    printf("thread_num=%d \n", omp_get_max_threads());
    #pragma omp parallel
    {
        int ID=omp_get_thread_num();
//        int ID=0;
        int clone_type=5; 
        long seed=master_seed+ID;        
        Genotype genotype_ori_copy;
        Genotype genotype_offspring;
        CellState state_offspring;
        GillespieRates rate_offspring;
        float init_mRNA_offspring[NGENES]; 
        float init_protein_conc_offspring[NGENES];
        float s,P_fix;
        char mut_type;
        
        initialize_cache(&genotype_offspring);
        initialize_cache(&genotype_ori_copy);
        
        while(i<MAX_MUT_STEP)
        {
            clone_type=5; /*copy all info*/           
            clone_cell(&genotype_ori, &genotype_ori_copy, clone_type);
            
            #pragma omp single
            {   
                fixation=0; 
                fprintf(OUTPUT,"Step %d, fitness=%f\n",ID,i,genotype_ori.fitness);
                printf("thread:%d, Step %d, fitness=%f\n",ID,i,genotype_ori.fitness);                   
            }
            
            while(!fixation)
            {
// //               #pragma omp critical
                clone_cell(&genotype_ori_copy,&genotype_offspring,clone_type); 
             
                clone_type=mutate(&genotype_offspring,kdis,&seed,&mut_type); 
                
                calc_all_binding_sites(&genotype_offspring);

                calc_avg_growth_rate(&genotype_offspring,
                                        &state_offspring,
                                        init_mRNA_offspring,
                                        init_protein_conc_offspring,                                               
                                        &rate_offspring,
                                        maxbound2,
                                        maxbound3,
                                        &seed); 
                
//                printf("fitness=%f\n",genotype[ID].fitness);

                #pragma omp critical
                {
                    if(!fixation)
                    { 
                        s=(genotype_offspring.fitness-genotype_ori_copy.fitness)/genotype_ori_copy.fitness;
                        
                        if (fabs(s)<EPSILON)
                            P_fix =(float)1/POP_SIZE;	
                        else 
                            P_fix =(float)(1-exp(-s))/(1-exp(-s*POP_SIZE)); 
                        
                        if(P_fix>ran1(&seed))
                        {
                            printf("mutation=%c\n",mut_type);
                            fprintf(OUTPUT,"mutation=%c\n",mut_type);
                            fixation=1;
                            i++;
                            clone_type=6;
                            clone_cell(&genotype_offspring, &genotype_ori, clone_type);
                            
                            /*double replicates to increase accuracy*/
                            calc_avg_growth_rate(&genotype_ori,
                                                    &state_ori,
                                                    init_mRNA,
                                                    init_protein_conc,                                                                
                                                    &rates_ori,
                                                    maxbound2,
                                                    maxbound3,
                                                    &master_seed);
                            genotype_ori.fitness=(float)sqrt((genotype_ori.avg_G1+genotype_offspring.avg_G1)*(genotype_ori.avg_G2+genotype_offspring.avg_G2))/(2*N_replicates);
                        }
                    }
                }
//                try_fixation(&genotype_ori, &genotype_offspring, genotype_ori_copy.fitness, &fixation, &i, &seed);   
            }
            #pragma omp barrier
        }
    } 
    t2=omp_get_wtime();
    fprintf(OUTPUT,"runtime=%f\n",t2-t1);
    printf("%f",t2-t1);  
}

//void print_time_course(TimeCourse *start,
//                       int i,
//                       int j)
//{
//  FILE *fpout;
//  char filename[80];
//  
//  /* do the normal thing on the first cell */
//  if (POP_SIZE == 1)
//    sprintf(filename, "%s/p%d.dat", output_directory, i);
//  else
//    sprintf(filename, "%s/p%d-%d.dat", output_directory, j, i);
//  if ((fpout = fopen(filename,"w"))==NULL) {
//    LOG_ERROR_NOCELLID("error: Can't open %s file\n", filename);
//  }
//  while (start) {
//    fprintf(fpout,"%g %g\n", start->time, start->concentration);
//    start = start->next;
//  }
//  fclose(fpout);  
//}
////
//void print_all_protein_time_courses(TimeCourse *timecoursestart[2][NPROTEINS],
//                                    TimeCourse *timecourselast[2][NPROTEINS])
//{
////  int i, j;
////  for (j = 0; j < 1; j++) {
////    for (i=0; i < NPROTEINS; i++) {
////      if ((output)) print_time_course(timecoursestart[j][i], i, j);
////    }
////  }
//print_time_course(timecoursestart[0][SELECTION_GENE_A],1,0);
//print_time_course(timecoursestart[0][SELECTION_GENE_B],2,0);
//}
