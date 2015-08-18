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
const float COOPERATIVITY=1.0;        /* dGibbs, relative to 1 additional specific nt */
const float COOPERATIVE_DISTANCE=11;  /* distance co-operativity operates (changed from 20) */ 
const float NUMSITESINGENOME = 1.3e+6; /* updated from 1.8e+6 */

const float mN = 3.3e-10 * 5e6;    /* Lynch et al. (2008), Tsai et al. (2008)  */

const int avg_protein_conc = 12000;
const float penalty = 0.0;


const float avg_inset=5.0;
const float avg_delet=5.0;
/* below are default options, can be changed on command line*/

float kon=1e-4;              /* lower value is so things run faster */
                             /* actual value should be kon=0.2225 is
                                based on 1 molecule taking
                                240seconds=4 minutes and 89% of the
                                proteins being in the nucleus*/
float kon_after_burnin=1e-4; /* lower value value after burn is so things run faster */
float Koff[TF_ELEMENT_LEN-NMIN+1];

int burn_in = 1;             /* disable burn-in by default */

float tdevelopment = 30.0;/* default  development time: can be changed at runtime */
float timemax = -1.0;      /* set an upper limit to development time (default to -1.0=no limit) */
int current_ploidy = 1;    /* ploidy can be changed at run-time: 1 = haploid, 2 = diploid */
int output = 0;
long seed = -1808;        
int dummyrun = 10;          /* used to change seed */

float growth_rate_scaling = 2.0; /* set default growth rate scaling factor */
float duration_env0 = 10.0; // in minutes
float duration_env1 = 1.0;
int N_replicates=10;
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
  float hc, gpeak, Ltf;
  gpeak = 0.005776*growth_rate_scaling;  /* in min^-1 based on doubling time of 120 min: ln(2)/(120 min)=0.005776 */
  Pp_a = 12000;
  Pp_b = 12000;              /* mean gene expression of all proteins is 12064.28 */
  Ltf= 1418;               /* mean gene expression of only TFs is 1418 */
  hc = (gpeak/avg_protein_conc)*(1-(log(2-2*0.2)/log(2)));      /* in min^-1 cost of doubling gene expression, based on Wagner (2005) 
                                                   * using {s=0.2, N=500} matches {s=10^-5, N=10^7} combination (both Ns=100) */
  h = hc/0.023;            /* using c=0.023/min from mean of distribution from Belle et al (2006)*/
  gmax_a = gpeak + hc*(Pp_a+(TFGENES*Ltf));    /* compute the gmax coefficient based on gpeak and other parameters */
  gmax_b = gpeak + hc*(Pp_b+(TFGENES*Ltf));
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

  LOG_VERBOSE_NOCELLID("len=%d, NGENES=%d, ploidy=%d, current_element=%d\n", len, NGENES, ploidy, current_element); 
  for (i=0; i<len/ploidy; i++) {
    first = (i / current_element)*ploidy*current_element + i % current_element;
//    second = first + current_element;
//    third = second + current_element;
//    fourth = third + current_element;
    LOG_VERBOSE_NOCELLID("first=%d, second=%d, third=%d, fourth=%d\n", first, second, third, fourth); 
    x = ran1(&seed);
    
    Seq[first] = set_base_pair(x);
    /* clone the randomly chosen sequence for all other sites */
//    Seq[second] = Seq[first];
//    Seq[third] = Seq[first];
//    Seq[fourth] = Seq[first];
  }
  LOG_VERBOSE_NOCELLID("length: %d, sequence is %s\n", strlen(Seq), Seq);
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
void print_tf_occupancy(CellState *state,
                        AllTFBindingSites *all_binding_sites,
                        float t)
{
  int bound_count[NGENES][MAX_COPIES];
  int i, j, gene_id, gene_copy;

  for (i = 0; i < NGENES; i++) 
    for (j = 0; j < MAX_COPIES; j++) 
      bound_count[i][j] = 0;
  
  /* for all the currently bound sites */
  for (j = 0; j < state->tf_bound_num; j++) {
    gene_id = all_binding_sites[state->tf_bound_indexes[j]].cisreg_id;
    gene_copy = all_binding_sites[state->tf_bound_indexes[j]].gene_copy;
    bound_count[gene_id][gene_copy]++;
  }
  fprintf(fp_tfsbound[state->cell_id], "%g %d ", t, state->tf_bound_num);
  for (i = 0; i < NGENES; i++) 
    for (j = 0; j < MAX_COPIES; j++) 
      fprintf(fp_tfsbound[state->cell_id], "%d ", bound_count[i][j]);

  fprintf(fp_tfsbound[state->cell_id], "\n");

}

void print_rounding(CellState *state, GillespieRates *rates, float t)
{
  fprintf(fp_rounding[state->cell_id], "%g %d %d %d %d %d %d %d\n", 
          t, rates->koff_operations, rates->transport_operations, rates->mRNAdecay_operations, 
          rates->pic_disassembly_operations, rates->salphc_operations, rates->max_salphc_operations, rates->min_salphc_operations);
}
#endif

void initialize_genotype_fixed(Genotype *genotype, 
                               float kdis[],
                               int genotype_id)
{
  int i, j, p;

  LOG_NOCELLID("[genotype %03d] activators vs repressors ", genotype_id);
  
  for (i=0; i < NGENES; i++) {
    genotype->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    while (genotype->mRNAdecay[i]<0.0)
      genotype->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    genotype->proteindecay[i]=-1.0;
    while (genotype->proteindecay[i] < 0.0) {
      if (ran1(&seed) < 0.08421)
        genotype->proteindecay[i] = 0.0;
      else genotype->proteindecay[i] = exp(0.7874*gasdev(&seed)-3.7665);
    }
    /* dilution no longer done here, because it is now variable (function of instantaneous growth rate) */
    genotype->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    while (genotype->translation[i] < 0.0)
      genotype->translation[i] = exp(0.7406*gasdev(&seed)+4.56);

    /* make the activations the same in each copy */
    genotype->N_act=0;
    genotype->N_rep=0;
    
    if (ran1(&seed)<PROB_ACTIVATING) {
      genotype->N_act++;  
      for (p=0; p < MAX_COPIES; p++) 
        genotype->activating[i][p] = 1;
    } else {
      genotype->N_rep++;
      for (p=0; p < MAX_COPIES; p++) 
        genotype->activating[i][p] = 0;
    }

    for (p=0; p < MAX_COPIES; p++) 
      LOG_NOFUNC("%d ", genotype->activating[i][p]);

    j = trunc(NUM_K_DISASSEMBLY * ran1(&seed));
    
    for (p=0; p < MAX_COPIES; p++) 
      genotype->pic_disassembly[i][p] = kdis[j];
  }
  LOG_NOFUNC("\n");
 
}

/*
 * initialize the genotype, this initializes random cis-regulatory
 * sequences for each individual, but the same random TF sequence,
 * hindrance positions, replication times, etc.  (full list below)
 */
void initialize_genotype(Genotype *genotype, 
                         Genotype *clone,
                         float kdis[],
                         int genotype_id)
{ 
  int i,j,k,p;
  
  initialize_sequence((char *)genotype->cisreg_seq, CISREG_LEN*MAX_COPIES*NGENES, MAX_COPIES, NGENES);
  
  initialize_sequence((char *)genotype->tf_seq, TF_ELEMENT_LEN*MAX_COPIES*TFGENES, MAX_COPIES, TFGENES);
 
  /* We now generate the complementary sequence of BS that are on the non-template strand.
   * The complementary sequence is used to search for BS that on the non-template strand.  
   * We also assume that all the TFs have strong orientation preference.*/  
    for(i=0;i<TFGENES;i++)
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
       
    initialize_genotype_fixed(genotype, kdis, genotype_id);

  /* start number of copies of gene at current_ploidy */
  for (p=0; p < NGENES; p++) {
    genotype->copies[p] = current_ploidy;
  }

  calc_all_binding_sites(genotype);
}

/*
 * compute the list binding sites for specified gene and gene copy
 */
void calc_all_binding_sites_copy(Genotype *genotype, int gene_id, float Koff[TF_ELEMENT_LEN-NMIN+1])
{
  int i, j, k, match,match_rc;
  int N_hindered_BS=0;
  int N_binding_sites=0;
  
  genotype->N_act_BS[gene_id]=0;
  genotype->N_rep_BS[gene_id]=0;
  genotype->max_hindered_sites[gene_id]=0;
  
  //make the code looks cleaner
  #define p_all_BS genotype->all_binding_sites[gene_id]

  //some helper pointer 
  char *tf_seq;
  char *cis_seq;
  char *tf_seq_rc; 
  cis_seq=&(genotype->cisreg_seq[gene_id][0][0]);  
  
  for (i=0; i < CISREG_LEN-TF_ELEMENT_LEN; i++) /* scan forwards */
  {  
        if(i>=TF_ELEMENT_LEN+2*HIND_LENGTH)
        {
            N_hindered_BS=0;
            
            for(j=0;j<CISREG_LEN*(genotype->N_act+genotype->N_rep);j++)
            {
               if(p_all_BS[j].BS_pos<i+1 && p_all_BS[j].BS_pos> i-TF_ELEMENT_LEN-2*HIND_LENGTH+1)
                   N_hindered_BS++; 
            }
        }
      
        for (k=0; k < TFGENES; k++) /* only loop through TF genes */
        {   
            tf_seq=&(genotype->tf_seq[k][0][0]);
            tf_seq_rc=&(genotype->tf_seq_rc[k][0][0]); 
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
                      LOG_ERROR_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d failed.\n", max_binding_site_alloc);
                      exit(1);
                    }
                    else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to binding_sites_num = %d succeeded\n", max_binding_site_alloc);
                }
                
                p_all_BS[N_binding_sites].tf_id = k;
                p_all_BS[N_binding_sites].Koff = Koff[TF_ELEMENT_LEN-match];
                p_all_BS[N_binding_sites].N_hindered = N_hindered_BS;
                p_all_BS[N_binding_sites].BS_pos = i ;
                N_binding_sites++;
                N_hindered_BS++;

                if(genotype->activating[k][0]==1) genotype->N_act_BS[gene_id]++;
            }

            /*find BS on the non-template strand*/
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
                    LOG_ERROR_NOCELLID("realloc of all_binding_sites to N_binding_sites = %d failed.\n", max_binding_site_alloc);
                    exit(1);
                  }
                  else LOG_VERBOSE_NOCELLID("realloc of all_binding_sites to N_binding_sites = %d succeeded\n", max_binding_site_alloc);
                }
                /************************************************************************************************************/
                p_all_BS[N_binding_sites].tf_id = k;
                p_all_BS[N_binding_sites].Koff = Koff[TF_ELEMENT_LEN-match_rc];
                p_all_BS[N_binding_sites].N_hindered = N_hindered_BS;
                p_all_BS[N_binding_sites].BS_pos = i;
                N_binding_sites++;
                N_hindered_BS++;

                if(genotype->activating[k][0]==1) genotype->N_act_BS[gene_id]++;
               
            }
        }
        
        genotype->max_hindered_sites[gene_id]=(genotype->max_hindered_sites[gene_id]>N_hindered_BS)? genotype->max_hindered_sites[gene_id]:N_hindered_BS;
  } 
  genotype->binding_sites_num[gene_id]=N_binding_sites; 
  genotype->N_rep_BS[gene_id]=N_binding_sites-(genotype->N_act_BS[gene_id]);
}

/*
 * compute the list of binding sites for the specified number of gene
 * copies
 */
void calc_all_binding_sites(Genotype *genotype)
{
  int p; 
  int max_binding_site_alloc= MAXELEMENTS;
  int gene_id;  
   
  for(gene_id=0;gene_id<NGENES;gene_id++)
  {
    genotype->all_binding_sites[gene_id] = malloc(max_binding_site_alloc*sizeof(AllTFBindingSites));
    
    if (!(genotype->all_binding_sites[gene_id])) 
    {
        LOG_ERROR_NOCELLID("initial setting of all_binding_sites failed.\n");
        exit(1);
    }
    
    calc_all_binding_sites_copy(genotype,gene_id);
  } 
}

int mod(int a, int b) // b is the base
{
    if(a>=0)
        return (a%b);
    else
        return (abs(a+b)%b);
}

float calc_ratio_act_to_rep(AllTFBindingSites *BS_info,
                           int max_N_hindered_BS,
                           int N_BS,
                           int N_act_BS,
                           int N_rep_BS, 
                           int activating[NGENES][MAX_COPIES],
                           float protein_conc[NGENES])
{
    double ratio_matrices[max_N_hindered_BS+1][N_rep_BS+1][N_act_BS+1];   
    double transition_matrix[N_rep_BS+1][N_act_BS+1];
    double sum,prob_act_over_rep=0.0;
    float product_of_freq;    
    float Kon[TFGENES];      
    
    int pos_of_last_record;    
    int pos_next_record;
    int i,j,k,m,n; 
    
    /*calc Kon based on TF concentration*/
    for(i=0;i<TFGENES;i++)
    {
        Kon[i]=kon*protein_conc[i];
    }
    
    /* initializing matrices to all zeros */
    for(i=0;i<max_N_hindered_BS+1;i++)
    {
        for(j=0;j<N_rep_BS+1;j++)
        {
            for(k=0;k<N_act_BS+1;k++)
            {
                ratio_matrices[i][j][k]=0.0;
            }
        }
    }
    
    for(j=0;j<N_rep_BS+1;j++)
    {
        for(k=0;k<N_act_BS+1;k++)
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
    }
    else
    {
        ratio_matrices[pos_next_record][1][0]=Kon[BS_info[0].tf_id];       
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
                if(m-BS_info[m].N_hindered!=0)
                {
                  // suppose the first dimension of ratio_matrices is 10 (0-9), then the 11th ratio matrix should be put in 0 and  
                  // the 10th record is at 9. Note in gcc mod does not follow the conventional mathematical definition 
                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1); //find the closest BS that is not hindered                                              

                  for(i=0;i<N_rep_BS+1;i++)
                  {
                      transition_matrix[i][0]=ratio_matrices[pos_of_last_record][i][N_act_BS];

                      for(j=1;j<N_act_BS+1;j++)
                      {
                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i][j-1];
                      }
                  }
                }
                else
                {
                    for(i=0;i<N_rep_BS+1;i++)
                    {
                        for(j=0;j<N_act_BS+1;j++)
                        {
                            transition_matrix[i][j]=0.0;
                        }
                    }

                    transition_matrix[0][1]=1.0;
                }
                    
                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);  //find last record              

                for(i=0;i<N_rep_BS+1;i++)
                {
                    for(j=0;j<N_act_BS+1;j++)
                    {
                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
                                                      product_of_freq*transition_matrix[i][j];                            
                    }
                }
                break;

            case 0: // a BS of repressors
                if(m-BS_info[m].N_hindered!=0)
                {
                  pos_of_last_record=mod(pos_next_record-BS_info[m].N_hindered-1,max_N_hindered_BS+1);                                 

                  for(j=0;j<N_act_BS+1;j++)
                  {
                      transition_matrix[0][j]=ratio_matrices[pos_of_last_record][N_rep_BS][j];

                      for(i=1;i<N_rep_BS+1;i++)
                      {
                          transition_matrix[i][j]=ratio_matrices[pos_of_last_record][i-1][j];
                      }
                  }
                }
                else
                {
                    for(i=0;i<N_rep_BS+1;i++)
                    {
                        for(j=0;j<N_act_BS+1;j++)
                        {
                            transition_matrix[i][j]=0.0;
                        }
                    }

                    transition_matrix[1][0]=1.0;
                }
                    
                pos_of_last_record=mod(pos_next_record-1,max_N_hindered_BS+1);

                for(i=0;i<N_rep_BS+1;i++)
                {
                    for(j=0;j<N_act_BS+1;j++)
                    {
                        ratio_matrices[pos_next_record][i][j]=BS_info[m].Koff*ratio_matrices[pos_of_last_record][i][j]+
                                                      product_of_freq*transition_matrix[i][j];                            
                    }
                }
                break;
        }
    }

    sum=0.0;

    for(i=0;i<N_rep_BS+1;i++)
    {
        for(j=0;j<N_act_BS+1;j++)
        {
            sum+=ratio_matrices[pos_next_record][i][j];
        }
    }

    for(i=0;i<N_rep_BS+1;i++)
    {
        j=round(fabs((i-0.31)/0.33)); // need at least one act to transcribe                     
        for(;j<N_act_BS+1;j++)
        {
            prob_act_over_rep+=ratio_matrices[pos_next_record][i][j];
        }
    } 

    return prob_act_over_rep/sum;
     // end of the forward algorithm   
}

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
    printf("Out of memory\n");
    exit(1);
  }
  newtime->gene_id = i;
  newtime->copy = p;
  newtime->time = t;
  LOG_VERBOSE_NOCELLID("adding event at time=%f for gene=%d (copy=%d)\n", t, i, p);
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
    printf("Out of memory\n");
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
    printf("Out of memory\n");
    exit(1);
  }
  newtime->gene_id = gene_id;
  newtime->copy = gene_copy;
  newtime->time = t;
  LOG_VERBOSE_NOCELLID("adding event end at time=%f for gene=%d (copy=%d)\n", t, gene_id, gene_copy);
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
    LOG_ERROR_NOCELLID("In %d elements, couldn't find element %d to delete in gene %d (copy %d)\n",
                       j+1, i, gene_id, gene_copy);
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
                     int id,
                     int copies[NGENES],
                     float mRNAdecay[NGENES],
                     float meanmRNA[NGENES],
                     float init_protein_conc[NPROTEINS])
{
  int i, j, k, totalmRNA;
  float t;

  /* initialize the ID of the cell */
  state->cell_id = id;

  /* start cell size at 0.5 */
  state->cell_size = 1.0; //04/20 KUN

  /* initialize growth rate to zero (could also be based on 120 min doubling, i.e. 0.00578) */
  state->growth_rate = 0.0;

  state->mRNA_transcr_time_end = NULL;
  state->mRNA_transcr_time_end_last = NULL;
  state->mRNA_transl_time_end = NULL;
  state->mRNA_transl_time_end_last = NULL;
  state->env0_time_end = NULL;
  state->env0_time_end_last = NULL;
  state->env1_time_end = NULL;
  state->env1_time_end_last = NULL;

  for (i=0; i < NGENES; i++) {

    for (j=0; j < MAX_COPIES; j++) {
      state->active[i][j] = NUC_NO_PIC;// was ON_WITH_NUCLEOSOME
    }

    totalmRNA = (int) poidev(meanmRNA[i],&seed);
    state->mRNA_nuclear_num[i] = (int) bnldev(STARTNUCLEUS, totalmRNA, &seed);
    state->mRNA_cyto_num[i] = totalmRNA - state->mRNA_nuclear_num[i];
    state->mRNA_transl_cyto_num[i] = 0;

    for (k=0; k<state->mRNA_cyto_num[i]; k++) {
      t = expdev(&seed) / mRNAdecay[i];
      if (t < TTRANSLATION) {
        (state->mRNA_cyto_num[i])--;
        (state->mRNA_transl_cyto_num[i])++;
        LOG_VERBOSE("add translation event time=%g for gene=%d\n", (TTRANSLATION-t), i);
        add_fixed_event(i, -1, TTRANSLATION-t, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));
      }
    } 

    int total_mRNA_transcribing = (int) poidev(meanmRNA[i]*TTRANSCRIPTION*mRNAdecay[i], &seed);
    
    /* split it up evenly between the copies */
    int mRNA_copy1 = trunc(total_mRNA_transcribing/current_ploidy);
    int mRNA_copy2 = total_mRNA_transcribing - mRNA_copy1;

    for (j=0; j < MAX_COPIES; j++) {
      if (j < current_ploidy)  {
        state->mRNA_transcr_num[i][j] = (j==0) ? mRNA_copy1 : mRNA_copy2;
        LOG_VERBOSE_NOCELLID("initializing state->mRNA_transcr_num[%2d][%2d]=%d\n", i, j, state->mRNA_transcr_num[i][j]);
        for (k=0; k < state->mRNA_transcr_num[i][j]; k++)
          add_fixed_event(i, j, ran1(&seed)*TTRANSCRIPTION, &(state->mRNA_transcr_time_end), &(state->mRNA_transcr_time_end_last));
      } else {
        state->mRNA_transcr_num[i][j] = 0;
      }
    }
  }
  for (i=0; i < NPROTEINS; i++) {
    state->protein_conc[i] = init_protein_conc[i];
  }
    
  for (j=0; j < MAX_COPIES; j++) {  
    int pos = 0;
    for (i=0; i < NGENES; i++) {
      if (genotype->copies[i] > j) {  
        state->state_change_ids[ACETYLATION_STATE][j][pos] = i;
        LOG_VERBOSE("Initializing statechange gene=%d, ploidy=%d state_change_ids[%d][%d]=%d\n", i, j, j, 
                    pos, state->state_change_ids[ACETYLATION_STATE][j][pos]);
        pos++;
      }
    }
  }
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
  
  LOG_VERBOSE("change_mRNA_cytoplasm[%d]: mRNA=%d, transl rate=%g, protein decay=%g, salphc=%g\n", 
              i, state->mRNA_cyto_num[i], genotype->translation[i], state->konvalues[i][KON_PROTEIN_DECAY_INDEX], salphc); 

  state->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / state->konvalues[i][KON_PROTEIN_DECAY_INDEX];
  state->konvalues[i][KON_SALPHC_INDEX] = salphc;
}

/* 
 * for specified gene_id and gene_copy, tests whether criterion for
 * transcription is met
 */
//int ready_to_transcribe(int gene_id,
//                        int gene_copy,
//                        int *tf_bound_indexes,
//                        int tf_bound_num,
//                        AllTFBindingSites *all_binding_sites,
//                        int activating[NGENES][MAX_COPIES],
//                        int *on)
//{
//  int i, off;
//  
//  *on=off=0;
//  for (i=0; i < tf_bound_num; i++) {
//    if (gene_id==all_binding_sites[tf_bound_indexes[i]].cisreg_id &&
//        gene_copy==all_binding_sites[tf_bound_indexes[i]].gene_copy)
//      {
//        if (activating[all_binding_sites[tf_bound_indexes[i]].tf_id][gene_copy]) (*on)++;
//        else off++;
//    }
//  }
//  if ((float)off <= 0.33442*(float)(*on) + 0.31303) 
//    return (1);
//  else 
//    return (0);
//}

/* 
 * returns true if at least one activator is bound
 */
//int is_one_activator(int gene_id,
//                     int gene_copy,
//                     int *tf_bound_indexes,
//                     int tf_bound_num,
//                     AllTFBindingSites *all_binding_sites,
//                     int activating[NGENES][MAX_COPIES])
//{
//  int i;
//  
//  for (i=0; i < tf_bound_num; i++)
//    if (gene_id==all_binding_sites[tf_bound_indexes[i]].cisreg_id && 
//        gene_copy==all_binding_sites[tf_bound_indexes[i]].gene_copy &&
//        (activating[all_binding_sites[tf_bound_indexes[i]].tf_id][gene_copy])) 
//      return (1);
//  return (0);
//}

/* 
 * calculate the rates based upon the initialization of the genotype
 * and the cell state.  Note this is only appropriate if nothing is
 * bound
 */
void calc_from_state(Genotype *genotype,
                     CellState *state,
                     GillespieRates *rates,                     
                     float transport_rate[NGENES],
                     float mRNAdecay[NGENES]) 
{
//  int i, j, k;
//  float salphc; 
//  float protein_decay;
//
//  for (i=0; i < NPROTEINS; i++) {
//    /* if protein decay is otherwise going to fall below cut-off, use aging term */
//    
//    protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;
//    
//    salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);
//    
//    state->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / (protein_decay);
//    
//    state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
//    
//    state->konvalues[i][KON_SALPHC_INDEX] = salphc;
//    
//    LOG_VERBOSE("protein decay[%d]=%g\n", i, kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX]);
//  }  
//  
///* these rates are calced in calc_dt */
//  rates->transport=0.0;
//  rates->transport_operations = 0;
//  rates->mRNAdecay=0.0;   
//  rates->mRNAdecay_operations = 0;
//  rates->pic_disassembly=0.0;
//  rates->pic_disassembly_operations = 0;
//
//  for (i=0; i < NGENES; i++) {  	
////    LOG("after initializing kon_states for gene=%d, site_count=%d, sites_per_gene=%d, nkon=%d\n", 
////        i, kon_states->kon_list[i]->site_count, genotype->sites_per_gene[i], kon_states->nkon);
//    transport_rate[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
//    rates->transport += transport_rate[i];
//    rates->transport_operations++;
////    LOG_VERBOSE("mRNA_nuclear_num=%d, initializing transport[%d]=%g\n", state->mRNA_nuclear_num[i], i, transport[i]);
//  }
////  LOG_VERBOSE("initializing rates->transport=%g\n", rates->transport);
//
//  /* start all genes in acteylated state */
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
//  
//  /* first initialize everything at zero */
//  for (j=0; j < MAX_COPIES; j++) {
////    rates->acetylation_num[j]=0;
////    rates->deacetylation_num[j]=0;
////    rates->pic_assembly_num[j]=0;
//    rates->transcript_init_num[j]=0;
////    rates->pic_disassembly_num[j]=0;
//  }

//  /* now set the per copy acetylation rate  */
//  for (i=0; i < NGENES; i++) {
//    for (j=0; j < genotype->copies[i]; j++) {
////      rates->acetylation_num[j]++;
//    }
//  }

//  if (verbose) 
//    for (j=0; j < MAX_COPIES; j++) {
//      LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
//      LOG_ERROR("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
//    }
}

void calc_all_rates(Genotype *genotype,
                    CellState *state,
                    GillespieRates *rates,                     
                    int UPDATE_KONVALUES)
{
    int i, j, k;
    float salphc; 
    float protein_decay;

    if(UPDATE_KONVALUES) /*Since these values are also updated in update_protein_conc_cell_size, this is only used in initilization*/
    {
        for (i=0; i < NPROTEINS; i++) 
        {
          /* if protein decay is otherwise going to fall below cut-off, use aging term */

          protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;

          salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);

          state->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / (protein_decay);

          state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);

          state->konvalues[i][KON_SALPHC_INDEX] = salphc;

          LOG_VERBOSE("protein decay[%d]=%g\n", i, kon_states->konvalues[i][KON_PROTEIN_DECAY_INDEX]);
        }
    }
    
    /* reset rates and operations */
    rates->transport=0.0;
//    rates->transport_operations = 0;
    rates->mRNAdecay=0.0;   
//    rates->mRNAdecay_operations = 0;
    rates->pic_disassembly=0.0;
//    rates->pic_disassembly_operations = 0;
    
    for(i=0;i<NGENES;i++)
    {
        rates->acetylation_rate[i]=0.0;
        rates->deacetylation_rate[i]=0.0;
        rates->pic_assembly_rate[i]=0.0;
        rates->pic_disassembly_rate[i]0.0;
        rates->transport_rate[i]=0.0;
        rates->mRNAdecay_rate[i]=0.0;
    }

    /* update RNA transport rate*/
    for (i=0; i < NGENES; i++) 
    {  	
  //    LOG("after initializing kon_states for gene=%d, site_count=%d, sites_per_gene=%d, nkon=%d\n", 
  //        i, kon_states->kon_list[i]->site_count, genotype->sites_per_gene[i], kon_states->nkon);
        rates->transport_rate[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
        rates->transport += rates->transport_rate[i];
//        rates->transport_operations++;
    //    LOG_VERBOSE("mRNA_nuclear_num=%d, initializing transport[%d]=%g\n", state->mRNA_nuclear_num[i], i, transport[i]);
    }
    rates->subtotal+=rates->transport;
  //  LOG_VERBOSE("initializing rates->transport=%g\n", rates->transport);
    
    /* update mRNAdecay rate based on the total number of mRNAs in both
     cytoplasm (mRNA_cyto_num) and ones that have only just recently arrived
     (mRNA_transl_cyto_num) */
    for (i=0; i < NGENES; i++) 
    {
        rates->mRNAdecay_rate[i] = genotype->mRNAdecay[i] * ((float) state->mRNA_cyto_num[i] + (float) state->mRNA_transl_cyto_num[i]);
        rates->mRNAdecay += rates->mRNAdecay_rate[i];
//        rates->mRNAdecay_operations++;      
    }
    rates->subtotal+=rates->mRNAdecay;
    
    /* start all genes in acteylated state */
//    for (j=0; j < MAX_COPIES; j++) 
//    {  
//        int pos = 0;
//        
//        for (i=0; i < NGENES; i++) 
//        {
//            if (genotype->copies[i] > j) 
//            {  
//                state->state_change_ids[ACETYLATION_STATE][j][pos] = i;
//                
//                LOG_VERBOSE("Initializing statechange gene=%d, ploidy=%d state_change_ids[%d][%d]=%d\n", i, j, j, 
//                            pos, state->state_change_ids[ACETYLATION_STATE][j][pos]);
//                pos++;
//            }
//        }
//    }
    
    /* calc TF binding prob distribution*/
    for(i=0;i<NGENES;i++)
    {
        state->Pact[i]=calc_ratio_act_to_rep(genotype->all_binding_sites[i],
                                            genotype->max_hindered_sites[i],
                                            genotype->binding_sites_num[i],
                                            genotype->N_act_BS[i],
                                            genotype->N_rep_BS[i],
                                            genotype->activating,                            
                                            state->protein_conc);    
    }
    
    /* calc other rates*/
    for(i=0;i<NGENES;i++)
    {
        switch (state->active[i][0])
        {
            case NUC_NO_PIC:
                rates->acetylation_rate[i]=state->Pact[i]*ACETYLATE;
                rates->subtotal+=rates->acetylation_rate[i];
                break;
                
            case NO_NUC_NO_PIC:
                rates->deacetylation_rate[i]=(1-state->Pact[i])*DEACETYLATE;
                rates->pic_assembly_rate[i]=state->Pact[i]*PICASSEMBLY;
                rates->subtotal+=rates->deacetylation_rate[i];
                rates->subtotal+=rates->pic_assembly_rate[i];
                break;
                
            case PIC_NO_NUC: // Note: pic_disassembly_rate is a gene-specific constant, so it's defined in genotype
                rates->pic_disassembly_rate[i]=genotype->pic_disassembly[i][0];
                rates->pic_disassembly+=rates->pic_disassembly_rate[i]; 
                rates->transcript_init_rate[i]= 1;
                rates->transcript_init+=1;
                break;
        }        
    }
    
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
 * calculate the length of the next (Gillespie) timestep, dt
 */
//calc_dt(x, dt, rates, mRNAdecay, 
//          genotype->mRNAdecay,
//          state->mRNA_cyto_num,
//          state->mRNA_transl_cyto_num, 
//          state->Pact,
//          state->protein_conc,
//          state->cell_id);
void calc_dt(float *x,
             float *dt,
             GillespieRates *rates, 
             Genotype *genotype,
             float mRNAdecay[],
             float mRNAdecayrates[],
             int mRNA_cyto_num[],
             int mRNA_transl_cyto_num[],
             float Pact[NGENES],
             float protein_conc[NPROTEINS],
             int cell_id)
{
//  float tbound1, tbound2;
//  int i, j;
//
//  /* reset the subtotal rate (excludes konrate) for current step */
//  rates->subtotal=0.0;
//  /* reset mRNA decay rate */
//  rates->mRNAdecay=0.0;
//  /* hence reset the number of rounding operations for this particular rate */
//  rates->mRNAdecay_operations=0;
//
//  /* update mRNAdecay rate based on the total number of mRNAs in both
//     cytoplasm (mRNA_cyto_num) and ones that have only just recently arrived
//     (mRNA_transl_cyto_num) */
//  for (i=0; i < NGENES; i++) {
//    mRNAdecay[i] = mRNAdecayrates[i] * ((float) mRNA_cyto_num[i] + (float) mRNA_transl_cyto_num[i]);
//    rates->mRNAdecay += mRNAdecay[i];
//    rates->mRNAdecay_operations++;
//    Pact[i]=calc_ratio_act_to_rep(genotype->all_binding_sites[i],
//                                            genotype->max_hindered_sites[i],
//                                            genotype->binding_sites_num[i],
//                                            genotype->N_act_BS[i],
//                                            genotype->N_rep_BS[i],
//                                            genotype->activating,                            
//                                            protein_conc);    
//  }
//
//  /* recompute and cache the total rate in data structure */
//  rates->subtotal += rates->transport;
//  rates->subtotal += rates->mRNAdecay;
//  rates->subtotal += rates->pic_disassembly;
//  //LOG_ERROR_NOCELLID("rates subtotal = %f, salphc = %f\n", rates->subtotal, rates->salphc);
//  //printf("rates subtotal = %f\n rates koff = %f\n rates transport = %f\n rates mrna = %f\n rates pic = %f\n rates salphc = %f \n",rates->subtotal, rates->koff, rates->transport, rates->mRNAdecay, rates->pic_disassembly, rates->salphc);
//  
//  /* 
//   * convert the counts back into rates using the constants 
//   */
//  for (j=0; j < MAX_COPIES; j++) {
//    rates->subtotal += (float) rates->acetylation_num[j] * ACETYLATE;
//    rates->subtotal += (float) rates->deacetylation_num[j] * DEACETYLATE;
//    rates->subtotal += (float) rates->pic_assembly_num[j] * PICASSEMBLY;
//    rates->subtotal += (float) rates->transcript_init_num[j] * TRANSCRIPTINIT;    
//  } 
//  //printf("rates subtotal = %f\n",rates->subtotal);
//  //system("PAUSE");
//  
//  *dt = *x/rates->subtotal;
}

/*
 * end transcription: update the mRNAs in the nucleus, cytoplasm
 * etc. accordingly and delete the event from the queue
 */
void end_transcription(float *dt,
                       float t,
                       CellState *state,
//                       float transport[NGENES],
                       GillespieRates *rates)
{
  int i, j, total;
  
  /* recompute the delta-t based on difference between now and the
     time of transcription end */
  *dt = state->mRNA_transcr_time_end->time - t;

  if (verbose) {
    total = 0;
    for (i=0; i < NGENES; i++) 
      for (j=0; j < MAX_COPIES; j++) 
        total += state->mRNA_transcr_num[i][j];
    
    LOG_VERBOSE("\ntranscription event finishes out of %d possible t=%g dt=%g\n", total, t, *dt);
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

  /* add rate KRNA to transport and Gillespie rates */
//  transport[i] += KRNA;
  rates->transport += KRNA;
  rates->transport_operations++;

//  LOG_VERBOSE("add one new mRNA in nucleus, updating transport[%d]=%g, rates->transport=%g\n", i, transport[i], rates->transport);

}

void end_translation(CellState *state,float *dt,float t)
{
    int i;
    int N_translation_event;
    
    *dt = state->mRNA_transl_time_end->time - t;         /* make dt window smaller */
              
    N_translation_event=0; 

    /* count current number of mRNAs that have recently arrived in cytoplasm */
    for (i=0; i<NGENES; i++) N_translation_event += state->mRNA_transl_cyto_num[i];

    LOG_VERBOSE("translation event finishes out of %d possible t=%g dt=%g\n", N_translation_event, t, *dt); /* bug: dt can be negative */

    /* get identity of gene that has just finished translating */
    i=state->mRNA_transl_time_end->gene_id;   

    /* there is one less mRNA that has just finished translation */
    (state->mRNA_transl_cyto_num[i])--;   

    /* delete the event that just happened */
    LOG_VERBOSE("delete translation event that just happened at time=%g", t);
    
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
 * update the activity state of 'gene_id' based on the new balance of
 * bound transcription factors and the current state of the gene
 */
//void revise_activity_state(int gene_id,
//                           int gene_copy,
//                           Genotype *genotype,
//                           CellState *state,
//                           GillespieRates *rates)
//{
//  int transcriptrule=-1, oldstate=-1, numactive=-1;
//
//  transcriptrule = ready_to_transcribe(gene_id, gene_copy, 
//                                       state->tf_bound_indexes, 
//                                       state->tf_bound_num,
//                                       genotype->all_binding_sites,
//                                       genotype->activating,
//                                       &numactive);
//
//  /* get last state of transcription initiation */
//  oldstate = state->active[gene_id][gene_copy];
//  //printf("gene id = %d\n", gene_id);
//  //if(gene_id == 0){
//            // system("PAUSE");
//  //}
//  /*
//   * first set of rules:
//   * ACTIVATING TFs exceed REPRESSING TFs 
//   */
//   LOG_ERROR("transcript rule = %d, oldstate = %d, transcriptrule & oldstate = %d, OFF_FULL = %d\n", transcriptrule, oldstate, (transcriptrule) && oldstate, OFF_FULL);
//  /* OFF_FULL -> ON_WITH_NUCLEOSOME */
//  if ((transcriptrule) && oldstate==OFF_FULL){
//                       LOG_ERROR("UP ACE NUM!\n");
//                       //printf("A\n");
//    state->active[gene_id][gene_copy] = ON_WITH_NUCLEOSOME;
//    state->state_change_ids[ACETYLATION_STATE][gene_copy][rates->acetylation_num[gene_copy]] = gene_id;
//    (rates->acetylation_num[gene_copy])++; // rates will be updated in calc_dt
//  }
//  /* OFF_NO_PIC -> ON_NO_PIC */
//  if ((transcriptrule) && oldstate==OFF_NO_PIC) {
//                       //printf("B\n");
//    state->active[gene_id][gene_copy] = ON_NO_PIC;
//    remove_from_array(gene_id, DEACETYLATION_STATE,  state->state_change_ids[DEACETYLATION_STATE][gene_copy], 
//                      &(rates->deacetylation_num[gene_copy]), (int) 1);
//    if (numactive){
//      state->state_change_ids[PICASSEMBLY_STATE][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id;
//      (rates->pic_assembly_num[gene_copy])++;
//      LOG_ERROR("PIC NUM = %d\n", rates->pic_assembly_num[gene_copy]);
//    }
//  }
//  /* OFF_PIC -> ON_FULL */
//  if ((transcriptrule) && oldstate==OFF_PIC) {
//                      // printf("C\n");
//    state->active[gene_id][gene_copy] = ON_FULL;
//  }
//
//  /*
//   * second set of rules:
//   * REPRESSING TFs exceed ACTIVATING TFs 
//   */
//
//  /* ON_WITH_NUCLEOSOME -> OFF_FULL */
//  if (!(transcriptrule) && oldstate==ON_WITH_NUCLEOSOME) {
//                       // printf("D\n");
//    state->active[gene_id][gene_copy] = OFF_FULL;
//    LOG_VERBOSE("removing gene=%d, copy=%d from state_change_ids[ACETYLATION_STATE][%d]\n", gene_id, gene_copy, gene_copy);
//    remove_from_array(gene_id, ACETYLATION_STATE, state->state_change_ids[ACETYLATION_STATE][gene_copy], 
//                      &(rates->acetylation_num[gene_copy]), (int) 1);
//  }
//  
//  /* ON_NO_PIC -> OFF_NO_PIC */
//  if (!(transcriptrule) && oldstate==ON_NO_PIC){    
//                       // printf("E\n");      
//    state->active[gene_id][gene_copy] = OFF_NO_PIC;
//    remove_from_array(gene_id, PICASSEMBLY_STATE, state->state_change_ids[PICASSEMBLY_STATE][gene_copy], 
//                      &(rates->pic_assembly_num[gene_copy]), (int) 0);
//    state->state_change_ids[DEACETYLATION_STATE][gene_copy][rates->deacetylation_num[gene_copy]] = gene_id;
//    (rates->deacetylation_num[gene_copy])++;
//  }
//  /* ON_FULL -> OFF_PIC  */
//  if (!(transcriptrule) && oldstate==ON_FULL) {
//                       // printf("F\n");
//    state->active[gene_id][gene_copy] = OFF_PIC;
//  }
////system("PAUSE");
//  /* do remaining transitions:
//   * OFF_PIC -> OFF_NO_PIC
//   * ON_FULL -> ON_NO_PIC 
//   */
//  if ((state->active[gene_id][gene_copy]==OFF_PIC || state->active[gene_id][gene_copy]==ON_FULL) && numactive==0)
//    disassemble_PIC(state, genotype, gene_id, gene_copy, rates);
//
//  if (verbose && (oldstate!=state->active[gene_id][gene_copy])) {
//    LOG_VERBOSE("state change from %d to %d in gene %d, copy %d\n", oldstate, state->active[gene_id][gene_copy], gene_id, gene_copy);
//  }
//}

/*
 * remove the bound transcription factor from the specified
 * site_id_to_unbind
 */
//void remove_tf_binding(Genotype *genotype,
//                       CellState *state,
//                       GillespieRates *rates,
//                       KonStates *kon_states,
//                       int site_id_to_unbind,
//                       float koffvalues[],
//                       float t)
//{
//  int i, j, k, bound, site_id, gene_id, gene_copy;
//
//  i = 0;
//
//  /* given site 'site_id_to_unbind', look for the index in the list of bound sites */
//  while ((state->tf_bound_indexes[i] != site_id_to_unbind) && (i < state->tf_bound_num)) 
//    i++;
//  if (i == state->tf_bound_num) {  /* couldn't find the site */
//    LOG_ERROR("t=%g could not find site %d with %d possibilities\n Bound sites are\n",
//              t, site_id_to_unbind, state->tf_bound_num);
//    for (j = 0; j < state->tf_bound_num; j++)  {
//      LOG_NOFUNC("%d\n", state->tf_bound_indexes[j]);
//    }
//  }
//  else {
//    j = 0;
//    /* loop through the sterically hindered sites */
//    while (j < state->tf_hindered_num) {
//
//      /* check all sites hindered by binding to location 'site_id_to_unbind' */
//      if (state->tf_hindered_indexes[j][1] == site_id_to_unbind) {
//        k = bound = 0;
//
//        /* is anything else hindering the same site? */
//        while (bound == 0 && k < state->tf_hindered_num) {
//          if (state->tf_hindered_indexes[j][0] == state->tf_hindered_indexes[k][0] && j != k) 
//            bound=1;
//          k++;
//        }
//
//        /* if nothing else is hindering this site then allow site_id_to_unbind to be (re-)bound */
//        if (bound==0) {
//          site_id = state->tf_hindered_indexes[j][0];
//          LOG_VERBOSE("Site %d left_edge_pos %d on gene %d freed from steric hindrance\n",
//                      site_id, genotype->all_binding_sites[site_id].left_edge_pos, genotype->all_binding_sites[site_id].cisreg_id);
//
//          /* adjust rates by returning kon to pool */
//          add_kon(state->protein_conc[genotype->all_binding_sites[site_id].tf_id],
//                  kon_states->konvalues[genotype->all_binding_sites[site_id].tf_id][KON_SALPHC_INDEX],
//                  genotype->all_binding_sites[site_id].tf_id,
//                  site_id,
//                  rates,
//                  kon_states);
//        }
//
//        /* now we have one less sterically hindered site */
//        (state->tf_hindered_num)--;
//
//        if (j < state->tf_hindered_num) {
//          /* shorten array by moving the end of the array to the hole opened up by removed site */
//          state->tf_hindered_indexes[j][0] = state->tf_hindered_indexes[state->tf_hindered_num][0];
//          state->tf_hindered_indexes[j][1] = state->tf_hindered_indexes[state->tf_hindered_num][1];
//        }
//      } else { /* only increment if we haven't shortened array */
//        j++;
//      }
//    }    
//
//    /* reduce the koff rate by the amount */
//    rates->koff -= koffvalues[i];
//    rates->koff_operations++;
//
//    /* one less bound site */
//    (state->tf_bound_num)--;
//
//    /* shift end of array to hole opened up */
//    state->tf_bound_indexes[i] = state->tf_bound_indexes[state->tf_bound_num];
//
//    /* likewise with koffvalues */
//    koffvalues[i] = koffvalues[state->tf_bound_num];
//
//    /* find the gene and copy whose cisreg region has an unbinding event */
//    gene_id = genotype->all_binding_sites[site_id_to_unbind].cisreg_id;
//    gene_copy = genotype->all_binding_sites[site_id_to_unbind].gene_copy;
//    LOG_VERBOSE("Add site %d at left_edge_pos %d on gene %d copy %d freed by unbinding\n",
//                site_id_to_unbind, genotype->all_binding_sites[site_id_to_unbind].left_edge_pos, gene_id, gene_copy);
//
//    /* adjust kon */
//    add_kon(state->protein_conc[genotype->all_binding_sites[site_id_to_unbind].tf_id],
//            kon_states->konvalues[genotype->all_binding_sites[site_id_to_unbind].tf_id][KON_SALPHC_INDEX],
//            genotype->all_binding_sites[site_id_to_unbind].tf_id,
//            site_id_to_unbind,
//            rates,
//            kon_states);
//
//    /* adjust the state of the gene */
//    revise_activity_state(gene_id, gene_copy, genotype, state, rates);
//
//    /* when TF unbinds adjust the co-operativity at close sites */
//    scan_nearby_sites(site_id_to_unbind, genotype->all_binding_sites, state, rates, koffvalues, t);
//  }
//}

/*
 * do transcription factor binding at specified site_id
 */
//void attempt_tf_binding(Genotype *genotype,
//                        CellState *state,
//                        GillespieRates *rates,
//                        float **koffvalues,
//                        KonStates *kon_states,
//                        int *maxbound2,
//                        int *maxbound3,
//                        int site_id,
//                        float t)
//{
//  int gene_id, gene_copy, k, posdiff,i;
//
//  LOG_VERBOSE("kon1 event at site %d out of %d possible, %d TFs previously bound binding_sites_num=%d\n",
//              site_id, kon_states->nkon, state->tf_bound_num, genotype->binding_sites_num);
//
//  /* if we have run out of space, double memory  */
//  if (state->tf_bound_num >= *maxbound2){
//    (*maxbound2) *= 2;
//    state->tf_bound_indexes = realloc(state->tf_bound_indexes, (*maxbound2)*sizeof(int));
//    
//    /* do the copy */
//    *koffvalues = realloc(*koffvalues, (*maxbound2)*sizeof(float));
//    
//    for(i=*maxbound2/2;i<*maxbound2;i++){
//        state->tf_bound_indexes[i]=-1;
//        (*koffvalues)[i]=0.0;
//    }
//    /* check return value */
//    if (!state->tf_bound_indexes || !(*koffvalues)) {
//      LOG_ERROR("memory allocation error resetting maxbound2=%d\n", *maxbound2);
//      exit(1);
//    }
//  }
//
//  /* append the site to end of indexes */
//  state->tf_bound_indexes[state->tf_bound_num] = site_id;
//  LOG_VERBOSE("remove site %3d on gene %2d (copy %d)\n", 
//              site_id, genotype->all_binding_sites[site_id].cisreg_id, genotype->all_binding_sites[site_id].gene_copy);
//
//  /* remove the site_id from the kon pool */
//  remove_kon(site_id,
//             genotype->all_binding_sites[site_id].tf_id,
//             rates, 
//             kon_states->konvalues[genotype->all_binding_sites[site_id].tf_id][KON_SALPHC_INDEX],
//             kon_states,
//             state->protein_conc[genotype->all_binding_sites[site_id].tf_id]);
//
//  /* recompute the koffvalues */
//  calc_koff(site_id, genotype->all_binding_sites, state, &((*koffvalues)[state->tf_bound_num]), t);
//
//  LOG_VERBOSE("new koff = %g is number %d\n",
//              (*koffvalues)[state->tf_bound_num], (state->tf_bound_num+1));
//
//  /* adjust rates by adding the new koffvalue to rates->koff */
//  rates->koff += (*koffvalues)[state->tf_bound_num];
//  rates->koff_operations++;
//  /* append site_id to list of bound sites */
//  state->tf_bound_indexes[state->tf_bound_num] = site_id;
//  /* increment number of bound TFs */
//  (state->tf_bound_num)++;
//  /* adjust co-operative binding in context of new TF */
//  scan_nearby_sites(site_id, genotype->all_binding_sites, state, rates, *koffvalues, t);
//  /* get the gene that the TF is binding to */
//  gene_id = genotype->all_binding_sites[site_id].cisreg_id;
//  /* get the copy that the TF is binding to */
//  gene_copy = genotype->all_binding_sites[site_id].gene_copy;
//  
//  /* update steric hindrance data structures */
//  /* JM: this cycles over all sites, not just bound ones, in order to
//     record redundancy in steric hindrance*/
//  for (k = 0; k < genotype->binding_sites_num; k++) {
//    /* if we are on the same gene and not the same binding site */
//    if (gene_id == genotype->all_binding_sites[k].cisreg_id &&
//        gene_copy == genotype->all_binding_sites[k].gene_copy &&
//        !(k==site_id)) {
//
//      /* check distance from current binding site (k) to the original (site_id) */
//      LOG_VERBOSE("site_id=%d k=%d\n", genotype->all_binding_sites[site_id].left_edge_pos, genotype->all_binding_sites[k].left_edge_pos);
//      posdiff = genotype->all_binding_sites[site_id].left_edge_pos - genotype->all_binding_sites[k].left_edge_pos;
//
//      /* if within HIND_LENGTH, we prevent future binding by adding to steric hindrance */
//      if (abs(posdiff) < HIND_LENGTH) {
//        /* if not enough memory, reallocate */
//        if (state->tf_hindered_num >= *maxbound3) {
//          (*maxbound3) *= 2;
//          state->tf_hindered_indexes = realloc(state->tf_hindered_indexes,2*(*maxbound3)*sizeof(int));
//         
//          for(i=*maxbound3/2;i<*maxbound3;i++){
//            state->tf_hindered_indexes[1][i]=-1;
//            state->tf_hindered_indexes[0][i]=-1;
//          }
//        }
//        /* record hindrance: 'site_id' blocks 'k' */
//        state->tf_hindered_indexes[state->tf_hindered_num][0] = k;
//        state->tf_hindered_indexes[state->tf_hindered_num][1] = site_id;
//        /* update list of hindered count */
//        (state->tf_hindered_num)++;
//        LOG_VERBOSE("%d steric hindrance sites after %d blocks site %d\n", state->tf_hindered_num, site_id, k);
//
//        /* remove the kon from pool */
//        remove_kon(k,
//                   genotype->all_binding_sites[k].tf_id,
//                   rates,
//                   kon_states->konvalues[genotype->all_binding_sites[k].tf_id][KON_SALPHC_INDEX],
//                   kon_states,
//                   state->protein_conc[genotype->all_binding_sites[k].tf_id]);
//      }
//    }
//  }
//  LOG_VERBOSE("tf_bound_num=%d tf_hindered_num=%d maxbound2=%d maxbound3=%d\n",
//              state->tf_bound_num, state->tf_hindered_num, *maxbound2, *maxbound3);
//
//  /* gene activity may change as a result of binding */
//  revise_activity_state(gene_id, gene_copy, genotype, state, rates);
//}

/*
 * time course of [TF]s represented as array of TimeCourse lists.
 */
void add_time_points(float time,
                     float protein_conc[NPROTEINS],
                     TimeCourse **timecoursestart,
                     TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i < NPROTEINS; i++)
    add_time_point(time, protein_conc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

void add_integer_time_points(float time,
                             int protein_conc[NPROTEINS],
                             TimeCourse **timecoursestart,
                             TimeCourse **timecourselast)
{
  int i;
  
  for (i=0; i < NPROTEINS; i++)
    add_time_point(time, (float) protein_conc[i], &(timecoursestart[i]), &(timecourselast[i]));
}

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
                                float alpha_a, 
                                float s_mRNA_a,
                                float alpha_b, 
                                float s_mRNA_b,
                                float all_alpha[NGENES],
                                int all_s_mRNA[NGENES],
                                float P_a,
                                float P_b,
                                float P_next_a,
                                float P_next_b,
                                float t, 
                                float deltat,
                                float c_a,                                
                                float ect1_a,
                                float c_b,                                
                                float ect1_b,
				int env) 
{
  int i;
  float instantaneous_growth_rate=0.0;  /* this is returned from the function */
  float total_alpha_s = 0.0;
  float deltatprime, deltatrest;



//  LOG_VERBOSE_NOCELLID("P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//                       P, P_next, c, t, t+deltat, s_mRNA);
                       
  switch (env){ // indicator of environment
  
	  case 1:{ // protein A is necessary!
	  
			 /* choose the appropriate piecewise linear integral */
		  if (((P_a >= Pp_a) && (P_next_a >= P_a)) || ((P_next_a >= Pp_a) && (P_a >= P_next_a))) {          /* P > Pp throughout */
//		    if (verbose)
//		      LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);
		
		    *integrated_growth_rate = gmax_a * deltat;	    
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);	     
		  } 
		  
		  else if (((P_next_a <= Pp_a) && (P_next_a >= P_a)) || ((P_a <= Pp_a) && (P_a >= P_next_a))) {   /* P < Pp throughout */
//		    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_a, c_a, gmax_a, deltat, s_mRNA_a, P_a, Pp_a, ect1_a);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);

		  } 
		  
		  
		  else if ((Pp_a > P_a) && (P_next_a > Pp_a)) {    /* P < Pp up until t' then P > Pp */
		    deltatprime = compute_tprime(c_a, P_a, alpha_a, s_mRNA_a);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_a, c_a, gmax_a, deltatprime, s_mRNA_a, P_a, Pp_a, ect1_a);
		    *integrated_growth_rate += gmax_a * deltatrest;
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);

		  } 
		  
		  
		  else if ((P_a > Pp_a) && (Pp_a > P_next_a)) {   /* P > Pp up until t' then P < Pp */
		    deltatprime = compute_tprime(c_a, P_a, alpha_a, s_mRNA_a);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = gmax_a * deltatprime;
		    *integrated_growth_rate += compute_integral(alpha_a, c_a, gmax_a, deltatrest, s_mRNA_a, P_a, Pp_a, ect1_a);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_b, c_b, 1.0, deltat, s_mRNA_b, P_b, 1.0, ect1_b);

		  } 
		  
		  
		  else {
//		    LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//		                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
//		    LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
//		                       cell_id);		
			
		    exit(1);
		  }
		
		  /* compute instantaneous growth rate at t */
		  if (P_next_a < Pp_a)
		    instantaneous_growth_rate = gmax_a*P_next_a/Pp_a - penalty*P_next_b;
		  else
		    instantaneous_growth_rate = gmax_a - penalty*P_next_b;
			
			break;
	}
	case 0:{ // protein b is necessary!
		 /* choose the appropriate piecewise linear integral */
		  if (((P_b >= Pp_b) && (P_next_b >= P_b)) || ((P_next_b >= Pp_b) && (P_b >= P_next_b))) {          /* P > Pp throughout */
//		    if (verbose)
//		      LOG_VERBOSE_NOCELLID("case 1: P=%g, P_next=%g > Pp=%g\n", P, P_next, Pp);
		    *integrated_growth_rate = gmax_b * deltat;
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  else if (((P_next_b <= Pp_b) && (P_next_b >= P_b)) || ((P_b <= Pp_b) && (P_b >= P_next_b))) {   /* P < Pp throughout */
//		    LOG_VERBOSE_NOCELLID("case 2: P=%g, P_next=%g < Pp=%g\n", P, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_b, c_b, gmax_b, deltat, s_mRNA_b, P_b, Pp_b, ect1_b);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  
		  else if ((Pp_b > P_b) && (P_next_b > Pp_b)) {    /* P < Pp up until t' then P > Pp */
		    deltatprime = compute_tprime(c_b, P_b, alpha_b, s_mRNA_b);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 3: P=%g < Pp=%g until t'=%g (deltatprime=%g) then P_next=%g > Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = compute_integral(alpha_b, c_b, gmax_b, deltatprime, s_mRNA_b, P_b, Pp_b, ect1_b);
		    *integrated_growth_rate += gmax_b * deltatrest;
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  
		  else if ((P_b > Pp_b) && (Pp_b > P_next_b)) {   /* P > Pp up until t' then P < Pp */
		    deltatprime = compute_tprime(c_b, P_b, alpha_b, s_mRNA_b);
		    deltatrest = deltat - deltatprime;
//		    LOG_VERBOSE_NOCELLID("case 4: P=%g > Pp=%g until t'=%g (deltatprime=%g) then P_next=%g < Pp=%g\n", 
//		                         P, Pp, t+deltatprime, deltatprime, P_next, Pp);
		    *integrated_growth_rate = gmax_b * deltatprime;
		    *integrated_growth_rate += compute_integral(alpha_b, c_b, gmax_b, deltatrest, s_mRNA_b, P_b, Pp_b, ect1_b);
		    *integrated_growth_rate -= penalty*compute_integral(alpha_a, c_a, 1.0, deltat, s_mRNA_a, P_a, 1.0, ect1_a);
		  } 
		  
		  
		  else {
//		    LOG_ERROR_NOCELLID("[cell %03d] P=%g, P_next=%g, c=%g, t=%g (in min) t+deltat=%g (in min), s_mRNA=%g\n", 
//		                       cell_id, P, P_next, c, t, t+deltat, s_mRNA);
//		    LOG_ERROR_NOCELLID("[cell %03d] growth rate computation error: should not reach here.  Exiting\n", 
//		                       cell_id);
		    exit(1);
		  }
		
		  /* compute instantaneous growth rate at t */
		  if (P_next_b < Pp_b)
		    instantaneous_growth_rate = gmax_b*P_next_b/Pp_b - penalty*P_next_a;
		  else
		    instantaneous_growth_rate = gmax_b - penalty*P_next_a;
		
		  break;
	}
  }

  LOG_VERBOSE_NOCELLID("growth rate (variable %g)-", *integrated_growth_rate);

  /* compute the total cost of translation across all genes  */
  for (i=0; i < NGENES; i++) {
    total_alpha_s += all_alpha[i] * all_s_mRNA[i];
  }

  /* add constant term for integrated rate */
  *integrated_growth_rate += -h * deltat * (total_alpha_s);

  /* and instantaneous integrated rate */
  instantaneous_growth_rate += -h * (total_alpha_s);

  LOG_VERBOSE_NOFUNC("(constant %g) = (total %g)\n", (h*deltat*total_alpha_s), *integrated_growth_rate);

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
  int i;
  float ct, ect, ect1,ect1_a,ect1_b;
  float L_a, L_b;
  float instantaneous_growth_rate = 0.0;
  float integrated_growth_rate = 0.0;
  float adjusted_decay;
 
  rates->max_salphc = rates->min_salphc = 0.0;
  for (i=0; i < NPROTEINS; i++) {
//   printf("p%i=%f\n",i,protein_conc[i]);
    if (i == SELECTION_GENE_A)  /* if we are looking at the selection gene, record protein concentration before update */
        L_a = state->protein_conc[i];
    if (i == SELECTION_GENE_B)
	L_b = state->protein_conc[i];
	
    /* update protein decay rates due to dilution caused by growth */
    adjusted_decay = genotype->proteindecay[i] + state->growth_rate;
 
    /* if this results in a very small or zero decay rate, use protein aging term */
    if (adjusted_decay > protein_aging)
      state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = adjusted_decay;
    else 
      state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = protein_aging;

    /* print out warning if decay rates get too low */
    if (state->konvalues[i][KON_PROTEIN_DECAY_INDEX] < 1e-10) {
      LOG_WARNING("protein=%02d, protein_decay=%g, genotype->proteindecay=%g, protein_aging=%g\n", i, adjusted_decay, 
                  genotype->proteindecay[i], protein_aging);
    }

    LOG_VERBOSE("prot decay[%d]=%g\n", i, state->konvalues[i][KON_PROTEIN_DECAY_INDEX]);

    ct = state->konvalues[i][KON_PROTEIN_DECAY_INDEX]*dt;
    ect = exp(-ct);
    if (fabs(ct)<EPSILON) ect1=ct;
    else ect1 = 1-ect; 
	
    if (i == SELECTION_GENE_A) ect1_a=ect1;		
    if (i == SELECTION_GENE_B) ect1_b=ect1;		
	
	
    /* get the new protein concentration for this gene */
    protein_conc[i] = state->konvalues[i][KON_SALPHC_INDEX]*ect1 + ect*state->protein_conc[i];

    /* recompute the konvalues and max and min salphc */
    state->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - state->konvalues[i][KON_SALPHC_INDEX]) / state->konvalues[i][KON_PROTEIN_DECAY_INDEX];
    
  }

  if (*env==1) {  /* now find out the protein concentration at end of dt interval and compute growth rate */ // a is the required
      
      instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genotype->translation[SELECTION_GENE_A], state->mRNA_cyto_num[SELECTION_GENE_A], 
                                                            genotype->translation[SELECTION_GENE_B], state->mRNA_cyto_num[SELECTION_GENE_B],
                                                            genotype->translation, state->mRNA_cyto_num, 
                                                            L_a, L_b, state->protein_conc[SELECTION_GENE_A], state->protein_conc[SELECTION_GENE_B],t, dt, 
                                                            state->konvalues[SELECTION_GENE_A][KON_PROTEIN_DECAY_INDEX], ect1_a,
                                                            state->konvalues[SELECTION_GENE_B][KON_PROTEIN_DECAY_INDEX], ect1_b,1);
      /* us the integrated growth rate to compute the cell size in the next timestep */
       
      state->cell_size = (state->cell_size)*exp(integrated_growth_rate);
      
      fprintf(fp_cellsize[0], "%g %g\n", t, state->cell_size);
    }
    
    if (*env==0) {  /* now find out the protein concentration at end of dt interval and compute growth rate */ // b is the required
         
      instantaneous_growth_rate = compute_growth_rate_dimer(&integrated_growth_rate, 
                                                            genotype->translation[SELECTION_GENE_A], state->mRNA_cyto_num[SELECTION_GENE_A], 
                                                            genotype->translation[SELECTION_GENE_B], state->mRNA_cyto_num[SELECTION_GENE_B], 
                                                            genotype->translation, state->mRNA_cyto_num, 
                                                            L_a, L_b, state->protein_conc[SELECTION_GENE_A],state->protein_conc[SELECTION_GENE_B], t, dt, 
                                                            state->konvalues[SELECTION_GENE_A][KON_PROTEIN_DECAY_INDEX], ect1_a,
                                                            state->konvalues[SELECTION_GENE_B][KON_PROTEIN_DECAY_INDEX], ect1_b,0);
      /* use the integrated growth rate to compute the cell size in the next timestep */
       
      state->cell_size = (state->cell_size)*exp(integrated_growth_rate);
      
      fprintf(fp_cellsize[0], "%g %g\n", t, state->cell_size);
    }
  
//  if ((output) && (*timecourselast)->time < t+dt-0.1) 
//    add_time_points(t+dt, otherdata, timecoursestart, timecourselast);

  /* update the instantaneous growth rate for the beginning of the next timestep */
  state->growth_rate = instantaneous_growth_rate;  
}

/*
 * START
 * Functions that handle each possible Gillespie event 
 *
 */
void transport_event(GillespieRates *rates,
                     CellState *state,
                     Genotype *genotype,
                     float dt,
                     float t,
                     float x)
{
    int i;
    float temp_rate;
    float endtime = t + dt + TTRANSLATION; 

//    update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, rates, t, env);

    i = -1;
    temp_rate = 0.0;  

    /* choose gene product (mRNA) that gets transported to cytoplasm
       based on weighting in transport[] array */
    while (i < NGENES && x > temp_rate) 
    {
        i++;
        x -= transport_rate[i];
    }

    if (i >= NGENES) 
    {
        LOG_ERROR("[cell %03d] attempted to choose mRNA for gene=%d which doesn't exist\n", state->cell_id, i);
        exit(0);
    } 

    LOG_VERBOSE("do transport event mRNA from gene=%d from %d copies (x=%g)\n", i, state->mRNA_nuclear_num[i], x);

    (state->mRNA_nuclear_num[i])--;   /* one less mRNAs in nucleus */
    (state->mRNA_transl_cyto_num[i])++;   /* it has just arrived in cytoplasm, ready to be translated */

    /* add the endtime for translation */
    LOG_VERBOSE("add translation event endtime=%f for mRNA encoded by gene=%d \n", endtime, i);
    
    add_fixed_event_end(i, -1, endtime, &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));

//    rates->transport_rate[i] -= KRNA;   /* decrease transport frequency */ // UPDATE all rates in calc_all_rates
//
//    /* if below a threshold, make zero */
//    if (rates->transport_rate[i] < 0.1*KRNA) 
//      rates->transport_rate[i]=0.0;
//
//    /* adjust rates */
//    rates->transport -= KRNA;
    rates->transport_operations++;
//
//    /* do similar threshold check */
//    if (rates->transport < 0.1*KRNA) 
//      rates->transport=0.0;
}

void mRNA_decay_event(GillespieRates *rates, CellState *state, Genotype *genotype)
{
    int i = -1, j;
    float temp_rate = 0.0;

    /* loop through mRNA products, to choose the mRNA with the
       proportionally higher decay rate */
    while (i < NGENES-1 && x > temp_rate) 
    {
        i++;
        temp_rate += rates->mRNAdecay_rate[i];
    }
    
    if (x > temp_rate) 
    { /* JM: had some rounding errors with rates->mRNAdecay. Calculate in calc_dt, hopefully fixed now */
          LOG_WARNING("x=%g > temp_rate=%g out of rates->mRNAdecay=%g\n",
                      x, temp_rate, rates->mRNAdecay);
    }

    /* assume mRNA cytoplasm transport events equally likely */
    x = ran1(&seed)*((float) (state->mRNA_cyto_num[i] + state->mRNA_transl_cyto_num[i]));

//    update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, rates, t, env);
    
    /* decay mRNA in cytoplasm */
    if (x < (float)state->mRNA_cyto_num[i]) 
    {
        LOG_VERBOSE("mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
                    i, state->mRNA_cyto_num[i], state->mRNA_transl_cyto_num[i]);

        /* remove the mRNA from the cytoplasm count */
        (state->mRNA_cyto_num[i])--;  
        
        change_mRNA_cytoplasm(i, genotype, state, rates); 

    } 
    else 
    {
        /* decay mRNA in process of translating */
        x = ran1(&seed)*((float) state->mRNA_transl_cyto_num[i]);
        
        LOG_VERBOSE("mRNA decay event gene %d not from %d copies in cytoplasm but %f from %d copies translating\n",
                    i, state->mRNA_cyto_num[i], trunc(x), state->mRNA_transl_cyto_num[i]);

        /* delete this fixed event: this mRNA will never be translated */
        LOG_VERBOSE("delete fixed TRANSLATION EVENT at time =%d for gene=%d\n", (int) trunc(x), i);
        delete_fixed_event(i, -1, (int) trunc(x), &(state->mRNA_transl_time_end), &(state->mRNA_transl_time_end_last));

        /* remove the mRNA from the count */
        (state->mRNA_transl_cyto_num[i])--;
        
        if (verbose) 
            for (j=0; j < NGENES; j++) 
            {
                LOG_VERBOSE("%d copies of gene %d translating\n", state->mRNA_transl_cyto_num[j], j);
            }
    }
}

void histone_acteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype)
{
  int gene_id=-1;
  float x = ran1(&seed)*rates->acetylation;

  while(gene_id<NGENES && x>0.0)
  {
      gene_id++;
      x-=rates->acetylation_rate[gene_id];
  }

  LOG_VERBOSE("acetylation event gene %d (copy %d)\nstate change from %d to 4\n",
              gene_id, 0, state->active[gene_id][0]);
  
  if (state->active[gene_id][0] != NUC_NO_PIC) {
    LOG_ERROR("acetylation event on gene %d (copy %d) attempted from state %d\n", gene_id, 0, state->active[gene_id][0]);
  }

  /* update protein concentration and cell size */
//  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, rates, t,
////				timecoursestart, timecourselast,
//                                env);
  
  /* set state: eject nucleosome, but there is no PIC yet */
  state->active[gene_id][0] = NO_NUC_NO_PIC;
  
//  remove_from_array(gene_id, ACETYLATION_STATE, state->state_change_ids[ACETYLATION_STATE][gene_copy], &(rates->acetylation_num[gene_copy]), (int) 1);
  
//  if (is_one_activator(gene_id, gene_copy, state->tf_bound_indexes, state->tf_bound_num, 
//                       genotype->all_binding_sites, genotype->activating)) 
//  {
//    state->state_change_ids[PICASSEMBLY_STATE][gene_copy][rates->pic_assembly_num[gene_copy]] = gene_id; 
//    (rates->pic_assembly_num[gene_copy])++;
//  }
}

void histone_deacteylation_event(GillespieRates *rates, CellState *state, Genotype *genotype)
{
  int gene_id=-1; 
  float x = ran1(&seed)*rates->deacetylation;

  /* choose a particular gene and copy to change state */
  while(gene_id<NGENES && x>0.0)
  {
      gene_id++;
      x-=rates->deacetylation_rate[gene_id];
  }

  LOG_VERBOSE("deacetylation event gene %d (copy %d)\nstate change from %d to 1\n",
              gene_id, 0, state->active[gene_id][0]);
  
  if (state->active[gene_id][0] != NO_NUC_NO_PIC) {
    LOG_ERROR("deacetylation event attempted from state %d\n", state->active[gene_id][gene_copy]);
  }

//  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, rates, t,
////				 timecoursestart, timecourselast,
//                                env);
  /* set state: nucleosome returns */
  state->active[gene_id][0] = NUC_NO_PIC;
//  remove_from_array(gene_id, DEACETYLATION_STATE, state->state_change_ids[DEACETYLATION_STATE][gene_copy], &(rates->deacetylation_num[gene_copy]), (int) 1);
}

void assemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genotype)
{
  float x = ran1(&seed)*rates->pic_assembly;

  int gene_id=-1; 

  /* choose a particular gene and copy to change state */
  while(gene_id<NGENES && x>0.0)
  {
      gene_id++;
      x-=rates->pic_assembly_rate[gene_id];
  }
//  get_gene(rates->pic_assembly_num, (int)trunc(x), &gene_loc, &gene_copy);
//  int gene_id = state->state_change_ids[PICASSEMBLY_STATE][gene_copy][gene_loc];

  LOG_VERBOSE("PIC assembly event gene %d copy %d\nstate change from %d to 6\n",
              gene_id,0, state->active[gene_id][0]);

  if (state->active[gene_id][0] != NO_NUC_NO_PIC) {
    LOG_ERROR("PIC assembly event attempted from state %d\n", state->active[gene_id][0]);
  }

//  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, rates, kon_states,t, 
////				timecoursestart, timecourselast,
//                                env);
  
  /* turn gene fully on: ready for transcription and adjust rates */
  state->active[gene_id][0] = PIC_NO_NUC;
//  remove_from_array(gene_id, PICASSEMBLY_STATE, state->state_change_ids[PICASSEMBLY_STATE][gene_copy], 
//                    &(rates->pic_assembly_num[gene_copy]), (int) 1);
//  state->state_change_ids[TRANSCRIPTINIT_STATE][gene_copy][rates->transcript_init_num[gene_copy]] = gene_id;
  (rates->transcript_init_num[0])++;
//  state->state_change_ids[PICDISASSEMBLY_STATE][gene_copy][rates->pic_disassembly_num[gene_copy]] = gene_id;
//
//  (rates->pic_disassembly_num[gene_copy])++;
//
//  rates->pic_disassembly += genotype->pic_disassembly[gene_id][gene_copy];
//  rates->pic_disassembly_operations++;
}

void disassemble_PIC_event(GillespieRates *rates, CellState *state, Genotype *genotype,
//			TimeCourse **timecoursestart, TimeCourse **timecourselast,
                           float dt, float t, float x)
{
//    int gene_copy, gene_loc;
    int gene_id=-1;

    /* choose an appropriate gene copy to disassemble the PIC from */
    while (gene_id < NGENES*current_ploidy && x>0) 
    {
        gene_id++;       
        x -= rates->pic_disassembly_rate[gene_id];
    }
    
    if (gene_id==NGENES*current_ploidy) { LOG_ERROR("error in PIC disassembly\n"); }
    
    /* now get the gene_id */     
//    LOG_VERBOSE("PIC disassembly event in copy %d of gene %d\n", gene_copy, gene_id);
    /* do the disassembling */
//    disassemble_PIC(state, genotype, gene_id, gene_copy, rates);
    
    state->active[gene_id][0]=NO_NUC_NO_PIC; // UPDATE ALL RATES TOGETHER IN CALC_ALL_RATES    
}

void transcription_init_event(GillespieRates *rates, CellState *state, Genotype *genotype, float dt, float t)
{
  int gene_id=-1;  
  int x=ran1(&seed)*rates->transcript_init;
  
  while(gene_id<NGENES && x>0)
  {
      gene_id++;
      x-=rates->transcript_init_rate[gene_id];
  }
//  x /= TRANSCRIPTINIT;
//  LOG_ERROR("TRANSCRIPT INIT\n");
//  /* choose the gene and copy that gets transcribed */
//  get_gene(rates->transcript_init_num, (int)trunc(x), &gene_loc, &gene_copy);
//  gene_id = state->state_change_ids[TRANSCRIPTINIT_STATE][gene_copy][gene_loc];
  LOG_VERBOSE("transcription event gene %d, copy %d\n", gene_id, 0);
  LOG_ERROR("Gene id = %d\n", gene_id);

  if (state->active[gene_id][0] != PIC_NO_NUC) {
    LOG_ERROR("transcription event attempted from state %d\n", state->active[gene_id][0]);
  }

//  update_protein_conc_cell_size(state->protein_conc, state, genotype, dt, rates, kon_states, t, 
////				timecoursestart, timecourselast,
//                                env);

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
void get_gene(int rate_array[MAX_COPIES],
              int pos,
              int *gene_loc,
              int *gene_copy)
{
  int i = 0;
  int total_rate = 0;
  *gene_copy = -1;   /* haven't found the copy yet */

  while (i < MAX_COPIES && *gene_copy < 0) {
    if (pos < (total_rate + rate_array[i])) {
      *gene_copy = i;
      *gene_loc = pos - total_rate;
    } else {
      total_rate += rate_array[i];
      i++;
    }
  }
}

/*
 * recalibrate the rates and cached data structures (KonStates) from
 * the current cell state
 */
void recalibrate_cell(GillespieRates *rates,
                      CellState *state,
                      Genotype *genotype,
                      float mRNAdecay[NGENES],
                      float transport_rate[NGENES],
                      float dt) 
{
  int i, j; 
  float protein_decay;
  float salphc = 0.0;

  /* reset all rates and operations */
  rates->transport=0.0;
  rates->transport_operations=0;
  rates->mRNAdecay=0.0;
  rates->mRNAdecay_operations=0;
  rates->pic_disassembly=0.0;
  rates->pic_disassembly_operations=0;
  rates->subtotal=0.0;
  
  /* regenerate kon_states */
  for (i=0; i < NPROTEINS; i++) {
    /* if protein decay is otherwise going to fall below cut-off, use aging term */
    protein_decay = genotype->proteindecay[i] >= protein_aging ? genotype->proteindecay[i] : protein_aging;
    salphc = (float) (state->mRNA_cyto_num[i]) * genotype->translation[i] / (protein_decay);
    state->konvalues[i][KON_DIFF_INDEX] = (state->protein_conc[i] - salphc) / (protein_decay);
    state->konvalues[i][KON_PROTEIN_DECAY_INDEX] = (protein_decay);
    state->konvalues[i][KON_SALPHC_INDEX] = salphc;
  }

  for (i=0; i < NGENES; i++) {
    /* transport rates */
    transport_rate[i] = KRNA * (float) (state->mRNA_nuclear_num[i]);
    rates->transport += transport_rate[i];
    rates->transport_operations++;

    /* regenerate decay rates */
    mRNAdecay[i] = genotype->mRNAdecay[i] * ((float) state->mRNA_cyto_num[i] + (float) state->mRNA_transl_cyto_num[i]);
    rates->mRNAdecay += mRNAdecay[i];
    rates->mRNAdecay_operations++;
    
    state->Pact[i]=calc_ratio_act_to_rep(genotype->all_binding_sites[i],
                                            genotype->max_hindered_sites[i],
                                            genotype->binding_sites_num[i],
                                            genotype->N_act_BS[i],
                                            genotype->N_rep_BS[i],
                                            genotype->activating,                            
                                            state->protein_conc);
  }

  /* recompute and cache the total rate in data structure */
  rates->subtotal += rates->transport;
  rates->subtotal += rates->mRNAdecay;
//  rates->subtotal += rates->pic_disassembly;

  /* 
   * convert the counts back into rates using the constants 
   */
  for (j=0; j < MAX_COPIES; j++) {
    rates->subtotal += (float) rates->acetylation_num[j] * ACETYLATE;
    rates->subtotal += (float) rates->deacetylation_num[j] * DEACETYLATE;
    rates->subtotal += (float) rates->pic_assembly_num[j] * PICASSEMBLY;
    rates->subtotal += (float) rates->transcript_init_num[j] * TRANSCRIPTINIT;    
  } 
}

void clone_cell(Genotype *genotype_orig,                
                Genotype *genotype_clone,
				int clone_type)
{
  int i, p;

  /* clone the cis-regulatory sequence */
  if(clone_type==4 || clone_type == 1)
  memcpy(genotype_clone->cisreg_seq, genotype_orig->cisreg_seq, sizeof(char [NGENES][MAX_COPIES][CISREG_LEN]));

  /* clone all other aspects of genotype */
  if(clone_type==4 || clone_type == 3)
    {
        for (i=0; i < NGENES; i++) {

            for (p=0; p < MAX_COPIES; p++) {

              genotype_clone->pic_disassembly[i][p] = genotype_orig->pic_disassembly[i][p];

            }   

            genotype_clone->mRNAdecay[i] =  genotype_orig->mRNAdecay[i];
            genotype_clone->proteindecay[i] =  genotype_orig->proteindecay[i];
            genotype_clone->translation[i] =  genotype_orig->translation[i];
        }
    }	
  
  if(clone_type == 4)
    {	
        for (i=0; i < NGENES; i++) {
            genotype_clone->copies[i] = genotype_orig->copies[i];

            for (p=0; p < MAX_COPIES; p++){genotype_clone->activating[i][p]=  genotype_orig->activating[i][p];}		
        }

        for (i=0; i < TFGENES; i++) {genotype_clone->hindrance_positions[i] = genotype_orig->hindrance_positions[i];} /* clone the hindrance positions */
    }	
  /* clone the TF sequence */
  if(clone_type==4 || clone_type==2) memcpy(genotype_clone->tf_seq, genotype_orig->tf_seq,sizeof(char [TFGENES][MAX_COPIES][TF_ELEMENT_LEN]));
}

/*
 * diagnostic function to dump a copy of the GillespieRates and some
 * of the CellState to the current error file
 */
void log_snapshot(GillespieRates *rates,
                  CellState *state,
                  Genotype *genotype,
                  float mRNAdecay[NGENES],
                  float transport[NGENES],
                  float x,
                  float t)
{
  int i, p, nkon = 0;

  LOG("snapshot at time=%g:\n x=%g, koff=%g = %d (tf_bound_num) * %g (koff/tf_bound_num)\n transport=%g\n decay=%g\n",
      t, x, rates->koff, state->tf_bound_num, rates->koff/(float)state->tf_bound_num, 
      rates->transport, rates->mRNAdecay);
  LOG_NOFUNC(" rates->salphc=%g\n rates->max_salphc=%g rates->min_salphc=%g\n", rates->salphc, rates->max_salphc, rates->min_salphc);
  LOG_NOFUNC(" konrate=%g\n", konrate);
  LOG_NOFUNC(" pic_disassembly=%g\n kon=%g = %d * %g\n",
             rates->pic_disassembly, rates->salphc+(konrate), kon_states->nkon, (rates->salphc+(konrate))/(float)kon_states->nkon);
  
  for (p=0; p < MAX_COPIES; p++) {
    LOG_NOFUNC(" acetylation=%g (copy %d)\n deacetylation=%g (copy %d)\n PIC assembly=%g (copy %d)\n transcriptinit=%g (copy %d)\n",
               (float)rates->acetylation_num[p]*ACETYLATE, p, (float)rates->deacetylation_num[p]*DEACETYLATE, p, 
               (float)rates->pic_assembly_num[p]*PICASSEMBLY, p, (float)rates->transcript_init_num[p]*TRANSCRIPTINIT, p);
  }
  LOG_NOFUNC(" total rates=%g=%g+%g\n", rates->subtotal + (konrate), rates->subtotal, konrate);
  LOG_NOFUNC(" total free=%d + total bound=%d + total hindered=%d = total sites=%d\n", 
             kon_states->nkon, state->tf_bound_num, state->tf_hindered_num, genotype->binding_sites_num);
  
  LOG_NOFUNC("\n");
}

/*
 * run the model for a specified cell for a single timestep:
 *  - returns 0  if cell is not "dead" (i.e. rates haven't deteroriated to zero)
 *  - returns -1 if cell is "dead"
 */
int do_single_timestep(Genotype *genotype, 
                       CellState *state,                         
                       GillespieRates *rates, 
                       float *t,                       
                       float *transport_rate,
                       float mRNAdecay[NGENES],
                       float *x,
                       float *dt,                      
                       int maxbound2,
                       int maxbound3,
                       int no_fixed_dev_time, 
                       int burn_in,
		       int *env) 
{
    int i,j;
    int event;     /* boolean to keep track of whether FixedEvent has ended */
    int total;     /* total possible translation events */
    float fixed_time;   

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

    *x = expdev(&seed);        /* draw random number */

    *dt = *x/rates->subtotal;

    if (*dt < 0.0) {
      LOG_ERROR("dt=%g is negative after calc_dt, t=%g\n", *dt, *t);
      exit(-1);
    }
    
    LOG_VERBOSE("next stochastic event due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

    /* first check to see if a fixed event occurs in current t->dt window, or in tdevelopment if running for a fixed development time */
    fixed_time = no_fixed_dev_time ? (*t+*dt) : fmin(*t+*dt, tdevelopment);

    event = does_fixed_event_end(state->mRNA_transl_time_end,
                                 state->mRNA_transcr_time_end,
                                 state->env0_time_end,
                                 state->env1_time_end,
                                 fixed_time);
    while(event>0)
    {                
        do_fixed_event(genotype, state, rates, dt,*t, event, *env);
        
        update_protein_conc_cell_size(state->protein_conc, state, genotype, *dt, rates, t, *env);
    
        *t += *dt;                  /* advance time by the dt */
        
        *x -= *dt*rates->subtotal; //note that rates->subtotal hasn't been updated, so it's a historical rate

        /* update rates->subtotal and re-compute a new dt */
        calc_all_rates(genotype, state, rates, NO_KON_UPDATE);      
        
        *dt = *x/rates->subtotal;       
        
        LOG_VERBOSE("next stochastic event (2) due at t=%g dt=%g x=%g\n", *t+*dt, *dt, *x);

        fixed_time = no_fixed_dev_time ? (*t+*dt) : fmin(*t+*dt, tdevelopment);    

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
          
    if (*t+*dt < tdevelopment || no_fixed_dev_time)
    { 
      /* if the total rates falls below zero, we do an emergency recalibration of cell */
        if (!(rates->subtotal > 0.0)) 
        {
            log_snapshot(rates, state, genotype, *x, *t);
            LOG_ERROR("x should always be >0 t=%g (x=%g) rates->subtotal=%g, konrate=%g, recalibrate cell!\n", *t, *x, rates->subtotal, *konrate); 
            recalibrate_cell(rates, state, genotype, *dt); 
            log_snapshot(rates, state, genotype, mRNAdecay, transport_rate, *x, *t);

            /* if this still results in either zero or negative total rates,
               this most likely due the cell being "dead" no TFs bound, no
               activity etc.  We mark cell as "dead" in this case, and
               remove from queue. */
            if (!(rates->subtotal > 0.0)) 
            {  
                  LOG_ERROR("cell is effectively dead\n"); 
                  return -1;        /* return cell status as "dead" */
            }
        }  

        do_Gillespie_event(genotype, state, rates, *dt, *t);

        update_protein_conc_cell_size(genotype, state, rates, *dt, *t, *env);
        
        /* Gillespie step: advance time to next event at dt */
        *t += *dt;
        LOG_VERBOSE("dt=%g t=%g\n", *dt, *t);
    } 
    else 
    { 
        /* we will reach the end of development in dt */
        LOG_VERBOSE("finish at t=%g dt=%g\n", *t, *dt);

        /* do remaining dt */
        *dt = tdevelopment - *t;

        /* final update of protein concentration */
        update_protein_conc_cell_size(genotype, state, rates,*dt, *t, *env);    
                                       
        /* advance to end of development (this exits the outer while loop) */
        *t = tdevelopment;
    }
    
    return 0;
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
            end_transcription(dt, t, state, rates);
            LOG_ERROR("END TRANSCRIPTION\n");            
            break;

        case 2:     /* if a translation event ends */ // I don't think this event is useful!
            end_translation(state, dt, t);
            change_mRNA_cytoplasm(state->mRNA_transl_time_end->gene_id, genotype, state, rates); /* the number of mRNAs in cytoplasm affects binding */                
            break;

        case 3:     /*change env from 0 to 1*/  
            *dt = state->env0_time_end->time - t;
            *env = 1;
            delete_fixed_event_start(&(state->env0_time_end),&(state->env0_time_end_last));	 
            add_fixed_event(0,0,t+*dt+duration_env1,&(state->env1_time_end),&(state->env1_time_end_last));
            break;

        case 4:	/*change env from 1 to 0*/
            *dt = state->env1_time_end->time - t;
            *env = 0;
            delete_fixed_event_start(&(state->env1_time_end),&(state->env1_time_end_last));	 
            add_fixed_event(0,0,t+*dt+duration_env0,&(state->env0_time_end),&(state->env0_time_end_last));
            break;

        default:
            LOG_ERROR("event=%d should never get here\n", event);
            exit(-1);
            break;
    }    
}
 /* 
     * choose a new uniform random number weighted by the rate of all
     * Gillespie events, note that konrate is *not* included in
     * rates->subtotal, so it needs to be added here
     */
void do_Gillespie_event(Genotype *genotype,
                        CellState *state,
                        GillespieRates *rates,
                        float dt,
                        float t)
{
    float x;    
    
    x=ran1(&seed)*(rates->subtotal);    
    
    if (verbose) log_snapshot(rates, state,genotype, x, t);    
    
    if (x < rates->transport)   /* transportation event */ 
    {   
        LOG_ERROR("transport event\n");
        transport_event(rates, state, genotype, dt, t, x, env);
                
//        for (j=0; j < MAX_COPIES; j++) 
//        {
//            //LOG_VERBOSE("rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
//            LOG_ERROR("2rates->acetylation_num[%d]=%d\n", j, rates->acetylation_num[j]);
//        }

    } 
    else 
    {
        LOG_ERROR("x = %f\n", x); 
        x -= rates->transport;
       // LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
        
        if (x < rates->mRNAdecay)  /*STOCHASTIC EVENT: an mRNA decay event */
        { 
            LOG_ERROR("decay event\n");
            mRNA_decay_event(rates, state, genotype);
                                                    
        } 
        else 
        {
            LOG_ERROR("x = %f\n", x);            
            x -= rates->mRNAdecay;
           // LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
          
            if (x < rates->pic_disassembly) /* STOCHASTIC EVENT: PIC disassembly*/
            {
                LOG_ERROR("pic disassembly event\n");
                disassemble_PIC_event(rates, state, genotype);
            } 
            else 
            {
                LOG_ERROR("x = %f\n", *x);
                x -= rates->pic_disassembly;
//                LOG_ERROR("CHECKBOB! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );

                if (x < rates->acetylation)  /* acetylation*/
                {
                    LOG_ERROR("hist act event\n");
                    histone_acteylation_event(rates, state, genotype);
                } 
                else 
                { 
                    x-= rates->acetylation;
                  
//                    LOG_ERROR("CHECK! x = %f,  sum_rate_counts = %d, deacetylate = %f\n", *x, sum_rate_counts(rates->deacetylation_num), DEACETYLATE );
                    
                    if (x < rates->deacetylation)/* STOCHASTIC EVENT: histone deacetylation */ 
                    {
                        LOG_ERROR("deact event\n"); 
                        histone_deacteylation_event(rates, state, genotype);
                    } 
                    else 
                    {
                        x -= rates->deacetylation; 
                        
                        if (x < rates->pic_assembly)/* STOCHASTIC EVENT: PIC assembly*/
                        {
                            LOG_ERROR("pic assembly event\n");
                            assemble_PIC_event(rates, state, genotype);                            
                        } 
                        else 
                        {
                            x -= rates->deacetylation;
                            
                            if (x < (float)rates->transcript_init * TRANSCRIPTINIT) /* STOCHASTIC EVENT: transcription initiation */
                            {
                                LOG_ERROR("transcript init event time = %f\n", *t);
                                transcription_init_event(rates, state, genotype, *dt, *t);
                            } 
                            else 
                            {
                                /*
                                 * FALLBACK: shouldn't get here, previous
                                 * events should be exhaustive
                                 */

                                LOG_ERROR("[cell %03d] t=%g no event assigned: x=%g, rates->subtotal+konrate=%g, recalibrate cell\n", 
                                          state->cell_id, *t, *x, rates->subtotal + *konrate);

                                log_snapshot(rates, state, genotype, *x, *t);
                                recalibrate_cell(rates, state, genotype, *dt); 
                                log_snapshot(rates, state, genotype, *x, *t);
                                
                                if (*t > 1000.0) 
                                {
                                    LOG_ERROR("cell has been running for %g which is > 1000 mins\n", *t);
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

float calc_avg_growth_rate(int current_genotype,
                           Genotype *genotype, 
                           CellState *state, 
                           float init_mRNA[NGENES],
                           float init_protein_conc[NGENES],
                           float *t,                         
                           float transport_rate[NGENES],
                           float mRNAdecay[NGENES],
                           float *x,
                           float *dt,                           
                           GillespieRates *rates,
                           float maxbound2,
                           float maxbound3,
                           float temperature,
                           int no_fixed_dev_time)
{	
    int cell_status=0; 
    int i,j;
    int env=0;
    float avg_GR=0;   

    for(i=0;i<N_replicates;i++)
    {	 
        env=0;
        for (j=0; j < NGENES; j++) {
            init_protein_conc[j] = exp(1.25759*gasdev(&seed)+7.25669);
            init_mRNA[j] = exp(0.91966*gasdev(&seed)-0.465902);
        }

        initialize_cell(state, current_genotype, genotype->copies, genotype->mRNAdecay, init_mRNA, init_protein_conc);	
   	    
        calc_from_state(genotype, state, rates, kon_states, transport_rate, mRNAdecay);	
    	
        *t = 0.0;

        do_single_timestep(genotype, 
                           state,                          
                           rates, 
                           t,                           
                           transport_rate,
                           mRNAdecay,
                           x,
                           dt,                          
                           maxbound2,
                           maxbound3,
                           no_fixed_dev_time,                           
                           &env);

        while(*t<tdevelopment){
            cell_status = do_single_timestep(genotype, 
                                             state,                                            
                                             rates, 
                                             t,                                            
                                             transport_rate,
                                             mRNAdecay,
                                             x,
                                             dt,                                             
                                             maxbound2,
                                             maxbound3,
                                             no_fixed_dev_time,                                            
                                             &env);
            if(cell_status==-1) {
                printf("dead!\n");
                break;				
            }		 		    
       }
       avg_GR += log(state->cell_size)/tdevelopment;	   	   
       free_fixedevent(state);	   			   
    }  
    
    return avg_GR/N_replicates;
}

int try_fixation(float current_fitness, float new_fitness)
{	
    float s, P_fix, ref;
    s = (new_fitness-current_fitness)/current_fitness;

    if (fabs(s)<EPSILON){P_fix = 1/(float)POP_SIZE;}	
    else{ P_fix = (1-exp(-s))/(1-exp(-s*(float)POP_SIZE)); }

    ref=ran1(&seed);

    if(ref > P_fix) return 0;
    else return 1;	
}

int mutate(Genotype *genotype)
{
   int l_genome = NGENES*MAX_COPIES*CISREG_LEN;
 //  int l_BS = TFGENES*MAX_COPIES*TF_ELEMENT_LEN;
   char mut='c';
//   char str_buf[CISREG_LEN];
   char *temp1, *temp2;
   char n;   
   int pos_n,i,pos_g, pos_c;
   float random;
 //  int l_del;
 //  int extra_n;
   int offset;
   int inset_size=0, delet_size=0;   
  //two helper pointers 
   char *Genome;
   Genome= (char*)genotype->cisreg_seq;     
   char *BS;
   BS=(char*)genotype->tf_seq;      
 // TODO: write a code to draw mutations
    switch (mut)
    {
        case 's': //substitution        		
            random=ran1(&seed)*l_genome;					
            pos_n=floor(random);		
            n=set_base_pair(random);
            while (n == Genome[pos_n]){	
                random=ran1(&seed);
                n=set_base_pair(random);
            }	
            Genome[pos_n]=n;			
            return 1;
            break;
        		
        case 'i': // insertion        
            while(inset_size<=0)
            {
                random=gasdev(&seed)+avg_inset;			
                inset_size = round(random);
            }							
            random=ran1(&seed)*NGENES;			
            pos_g=floor(random);			
            random=ran1(&seed)*MAX_COPIES;			
            pos_c=floor(random);			
            random=ran1(&seed)*CISREG_LEN;			
            pos_n=floor(random);					
            if (pos_n+inset_size>CISREG_LEN) inset_size=CISREG_LEN-pos_n;

            for(i=pos_n;i<pos_n+inset_size;i++){
                random = ran1(&seed);					
                genotype->cisreg_seq[pos_g][pos_c][i]=set_base_pair(random);							
            }			
            return 1;
            break;	
        		
        case 'p': // partial deletion        
            while(delet_size<=0)
            {
                random=gasdev(&seed)+avg_delet;			
                delet_size = round(random);
            }						
            random=ran1(&seed)*NGENES;			
            pos_g=floor(random);			
            random=ran1(&seed)*MAX_COPIES;			
            pos_c=floor(random);			
            random=ran1(&seed)*CISREG_LEN;			
            pos_n=floor(random);

            if (pos_n+delet_size>=CISREG_LEN-1)	// if only the tail is deleted, treat it like insertion
            {
                delet_size=CISREG_LEN-1-pos_n;				
                for(i=pos_n;i<pos_n+delet_size;i++){								
                    random = ran1(&seed);					
                    genotype->cisreg_seq[pos_g][pos_c][i]=set_base_pair(random);
                }				
            }
            else // else, join the two fragments aside the deletion 
            {
                for(i=pos_n;i<CISREG_LEN-delet_size;i++){genotype->cisreg_seq[pos_g][pos_c][i]=genotype->cisreg_seq[pos_g][pos_c][i+delet_size];}                				
                for(i++;i<CISREG_LEN;i++){					
                    random=ran1(&seed);					
                    genotype->cisreg_seq[pos_g][pos_c][i]= set_base_pair(random);
                }						
            }			
            return 1;
            break;			
        		
        case 'w': // whole gene (only tf genes get deleted) deletion.       			
            random=ran1(&seed)*TFGENES;			
            pos_g=floor(random);				
            temp1 = &genotype->cisreg_seq[pos_g][0][0];			
            offset=MAX_COPIES*CISREG_LEN;

            for(i=0;i<CISREG_LEN*MAX_COPIES*(TFGENES-pos_g-1);i++){				
                *temp1=*(temp1+offset);
                temp1++;				
            }					
//            NGENES--;
//            TFGENES--;
//            NPROTEINS--;
//            SELECTION_GENE_A=NGENES-2;
//            SELECTION_GENE_B=NGENES-1;            
            return 1;
            break;
        		
        case 'd': // only tf genes get duplicated                   
//            random=ran1(&seed)*n_gene;
//            pos_g=floor(random);
//            char new_genome[n_gene+1][n_copy][l_element];
//            temp1=(char*)new_genome;
//            temp2=&Seq2[pos_g][0][0];//
//            for(i=0;i<len;i++){ //copy the original genome to the new genome
//
//                    temp1[i]=Seq1[i];
//            }
//            extra_n=n_copy*l_element;
//            for(i++;i<len+extra_n;i++){
//                    temp1[i]=*temp2;
//                    temp2++;				
//            }
//            random=ran1(&seed)*NGENES;
//            pos_g=floor(random);
//            temp1=&genotype->cisreg_seq[pos_g][0][0];
//            temp2=Genome+l_genome;
//            for(i=0;i<CISREG_LEN*MAX_COPIES;i++){*temp2++=*temp1++;}
//	    n_gene++; // only TF-gene get mutated in this way
            return 1;
            break;
        
        case 'c': //binding sequence        
            random=ran1(&seed)*TF_ELEMENT_LEN;			
            pos_n=floor(random);
            n=set_base_pair(random);				
            while (n == BS[pos_n]){	
                random=ran1(&seed);
                n=set_base_pair(random);
            }	
            BS[pos_n]=n;
            return 2;
            break;  
            
        case 'k': //mutations in kinetic constants        
            return 3;
            break;        
    }
}

void init_run_pop(Genotype genotype[N_para_threads],
                  CellState state[N_para_threads],
                  float temperature,   /* in Kelvin */
                  float kdis[NUM_K_DISASSEMBLY],
                  int output_binding_sites,
                  int no_fixed_dev_time)
{  
  int i;
  int current_genotype = 0;
  int new_genotype = 1;
  float current_fitness;
  float new_fitness;
  int fixation = 0;
  int maxbound2, maxbound3; 
  float init_mRNA[N_para_threads][NGENES], init_protein_conc[N_para_threads][NGENES];
  float t[N_para_threads];             /* time of last event */  
  float transport_rate[N_para_threads][NGENES];  /* transport rates of each mRNA */
  float mRNAdecay[N_para_threads][NGENES];  /* mRNA decay rates */
  float x[N_para_threads];                  /* random number */
  float dt[N_para_threads];                 /* delta-t */  
  int clone_type=4;
 
  GillespieRates rates[N_para_threads];
	
  maxbound2 = MAXBOUND;
  maxbound3 = 10*MAXBOUND;
  
  for(i=0;i<TF_ELEMENT_LEN-NMIN+1;i++)
  {      
      Koff[i]=NUMSITESINGENOME*kon*0.25*KR/exp(-((float)i/3.0-1.0));      
  }
  //output=1;
       
  initialize_genotype(&genotype[current_genotype], &genotype[current_genotype], kdis, current_genotype); // checked
  
  current_fitness = calc_avg_growth_rate(current_genotype,
                                         &genotype[current_genotype],
                                         &state[current_genotype],
                                         &init_mRNA[current_genotype][0],
                                         &init_protein_conc[current_genotype][0],
                                         &t[current_genotype],
                                         &transport_rate[current_genotype][0],
                                         &mRNAdecay[current_genotype][0],
                                         &x[current_genotype],
                                         &dt[current_genotype],
                                         &rates[current_genotype],
                                         maxbound2,
                                         maxbound3,
                                         temperature,
                                         no_fixed_dev_time);										 

for(i=0;i<MAX_MUT_STEP;i++){

    printf("step:%d, gr=%f \n",i,current_fitness);

    while(!fixation){
        clone_cell(&genotype[current_genotype],&genotype[new_genotype],clone_type); // use type to decide which element in genotype needs clone
        clone_type=mutate(&genotype[new_genotype]); // type 0: genome type 1: tf-element type 3: kinetic rate constant type 4: everything

        calc_all_binding_sites(&genotype[new_genotype]);

        new_fitness = calc_avg_growth_rate(new_genotype,
                                 &genotype[new_genotype],
                                 &state[new_genotype],
                                 &init_mRNA[new_genotype][0],
                                 &init_protein_conc[new_genotype][0],
                                 &t[new_genotype],                                 
                                 &transport_rate[new_genotype][0],
                                 &mRNAdecay[new_genotype][0],
                                 &x[new_genotype],
                                 &dt[new_genotype],                                
                                 &rates[new_genotype],
                                 maxbound2,
                                 maxbound3,
                                 temperature,
                                 no_fixed_dev_time);

        printf("n_gr=%f\n",new_fitness);
        fixation = try_fixation(current_fitness, new_fitness); // returns 0 if fails to fix				
    }		
    current_genotype = 1 - current_genotype; 
    new_genotype = 1 - new_genotype;
    current_fitness = new_fitness;
    fixation=0;		
 }
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

