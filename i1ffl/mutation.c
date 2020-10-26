/* 
 * This file contains functions to generate mutations and maintain data structure 
 * 
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
#include <stdlib.h>
#include <math.h>
#include "RngStream.h"
#include "netsim.h"
#include "numerical.h"
#include "mutation.h"

enum mutation_bound {BOUND_EXCLUDED=0,BOUND_INCLUDED=1};

/*Mutation rates*/
float DUPLICATION=1.5e-7;   
float SILENCING=1.5e-7; 
static const float SUBSTITUTION = 3.5e-10; // per bp 
static const float MUT_Kd=3.5e-9; //per gene
static const float MUT_identity=3.5e-9; //per gene
static const float MUT_binding_seq=3.5e-9; //per gene
static const float MUT_ACT_to_INT=3.5e-9; //per gene
static const float MUT_mRNA_decay=9.5e-12; //per codon
static const float MUT_protein_decay=9.5e-12; //per codon
static const float MUT_protein_syn_rate=9.5e-12; //per codon
static const float MUT_GENE_LENGTH_RATE=1.2e-11; //per codon
static const float MUT_cooperation=0.0; //per gene
static const float MUT_effector_to_TF=0.0; //per gene
static const float MUT_locus_specific_tf_behavior=0.0; //5.25e-10 per codon

/*Mutational effect*/
float miu_ACT_TO_INT_RATE=1.57;
float miu_Kd=-5.0;       
float miu_protein_syn_rate=0.021; 
static const float sigma_ACT_TO_INT_RATE=0.773; 
static const float sigma_mRNA_decay=0.396; 
static const float sigma_protein_decay=0.739; 
static const float sigma_protein_syn_rate=0.76; 
static const float sigma_Kd=0.776;
static const float miu_mRNA_decay=-1.19;
static const float miu_protein_decay=-1.88;
static const float mutational_regression_rate=0.5;

/******************************************************************************
 * 
 *                          Private function prototypes
 *
 *****************************************************************************/
static float mut_make_new_value(float, float, float, float, float, RngStream, Mutation *, int);

static void update_protein_pool(Genotype *, int, int, char);

static void rebuild_protein_pool(Genotype *);

static void update_node_family_pool(Genotype *, int, char);

static void rebuild_node_family_pool(Genotype *);

static void update_output_protein_pool(Genotype *, int, char);

static void regroup_cisreg_clusters(Genotype *);

static void rebuild_cisreg_cluster_pool(Genotype *);

/*******************************************************************************
 *
 *                              Global functions 
 * 
 ******************************************************************************/
/*single nucleic acid substitution in cis-reg*/
void mut_substitution(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_nucleotide, which_gene,i,N_BS_bf_mutation,flag_difference,cluster_id;
    char *Genome,nucleotide; 
    AllTFBindingSites *container; // used to store binding sites before substitution
                                  // we compare whether substituion changes any binding sites,
                                  // in order to determine whether mutation creates a unique cis-reg     
                                  // sequence whose binding configuration always needs computation     
    /*points to the current cis-reg*/
    Genome= &genotype->cisreg_seq[0][0];  
    /*this is the nucleic acid that's going to be mutated*/        
    which_nucleotide=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
    /*this is gene hit by the mutation*/
    which_gene=which_nucleotide/CISREG_LEN;
    /*If which_gene's cisreg_cluster contains multiple genes, there is a chance that the mutation creates a new cluster*/
    cluster_id=genotype->which_cluster[which_gene];
    if(genotype->cisreg_cluster_pool[cluster_id][0][0]>1)
    {
        /*calculate and store the distribution of binding sites before mutation*/
        calc_all_binding_sites_copy(genotype,which_gene,NMIN);
        container=malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
            container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
            container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
        }
    }
    
    /*generate new nucleic acid*/
    nucleotide=set_base_pair(RngStream_RandU01(RS));
    while(nucleotide==Genome[which_nucleotide]) // make sure we get a different nucleic acid
        nucleotide=set_base_pair(RngStream_RandU01(RS));
    Genome[which_nucleotide]=nucleotide;
   
    if(genotype->cisreg_cluster_pool[cluster_id][0][0]>1)   
    {
        /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster_pool*/
        calc_all_binding_sites_copy(genotype,which_gene,NMIN);    
        if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])    
            update_cisreg_cluster_pool(genotype,which_gene,'s');
        else
        {
            /*assuming no change in binding sites*/
            flag_difference=0;
            /*comparing binding sites pair by pair*/
            for(i=0;i<N_BS_bf_mutation;i++)
            {
                if(container[i].BS_pos!=genotype->all_binding_sites[which_gene][i].BS_pos ||
                    container[i].tf_id!=genotype->all_binding_sites[which_gene][i].tf_id ||
                    container[i].mis_match!=genotype->all_binding_sites[which_gene][i].mis_match)
                {    
                    flag_difference=1;
                    break;
                }
            }
            if(flag_difference==1)
               update_cisreg_cluster_pool(genotype,which_gene,'s'); 
        }
        free(container);
    }
    /*record mutation info*/
    mut_record->which_nucleotide=which_nucleotide;
    mut_record->which_gene=which_gene;
    mut_record->nuc_diff[0]=Genome[which_nucleotide];
    genotype->recalc_TFBS[which_gene]=YES;  /*the binding sites on this cis-reg needs to be updated*/    
}

/*
 * Replay substitution. Used when replay the evolution for analysis
 */
void reproduce_substitution(Genotype *genotype, Mutation *mut_record)
{    
    char *Genome;
    int i, N_BS_bf_mutation,which_gene,flag_difference, cluster_id;  
    AllTFBindingSites *container;   
    /*get the mutated gene from record*/
    which_gene=mut_record->which_gene;
    cluster_id=genotype->which_cluster[which_gene];
    if(genotype->cisreg_cluster_pool[cluster_id][0][0]>1)
    {
        /*calculate and store the distribution of binding sites before mutation*/
        calc_all_binding_sites_copy(genotype,which_gene,NMIN);
        container=malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
        N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
        for(i=0;i<N_BS_bf_mutation;i++)
        {
            container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
            container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
            container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
        }
    }
    /*apply the mutation from record*/
    Genome= &genotype->cisreg_seq[0][0];  
    Genome[mut_record->which_nucleotide]=mut_record->nuc_diff[1]; 
    if(genotype->cisreg_cluster_pool[cluster_id][0][0]>1)
    {
        /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster_pool*/
        calc_all_binding_sites_copy(genotype,which_gene,NMIN); 
        if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])
            update_cisreg_cluster_pool(genotype,which_gene,'s'); 
        else
        {
            flag_difference=0;
            for(i=0;i<N_BS_bf_mutation;i++)
            {
                if(container[i].BS_pos!=genotype->all_binding_sites[which_gene][i].BS_pos ||
                    container[i].tf_id!=genotype->all_binding_sites[which_gene][i].tf_id ||
                    container[i].mis_match!=genotype->all_binding_sites[which_gene][i].mis_match)
                {    
                    flag_difference=1;
                    break;
                }
            }
            if(flag_difference==1)
               update_cisreg_cluster_pool(genotype,which_gene,'s'); 
        }
        free(container);   
    }
    genotype->recalc_TFBS[which_gene]=YES;  /*recalc TFBS*/
}

/**
 *Deletion whole cis-reg sequence
 */
void mut_whole_gene_deletion(Genotype *genotype, Mutation *mut_record, RngStream RS) // any gene can be deleted
{
    int i,j;
    int which_gene, protein_id;
    char *temp;        
    /* check which genes can be deleted*/      
    if(genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF>MIN_NON_OUTPUT_GENES)//if there are more than one copies of genes of non-output genes  
    {
        if(genotype->n_output_genes==1)//if the effector have only one copies
        {            
            while(1) //the last copy of the effector cannot be deleted. choose a different gene
            {   
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);                
                if(genotype->is_output[which_gene]==NON_OUTPUT_PROTEIN)
                    break;
            }
        }
        else //otherwise any gene can be deleted
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    }
    else //if there's only one copy of non-output gene
    {        
        while(1)
        {
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);            
            if(genotype->is_output[which_gene]==OUTPUT_PROTEIN)//obviously we can only delete a copy of the effector gene
                break;                           
        }
    }  
    /*record mutation info*/
    mut_record->which_gene=which_gene;     
    /*points the first nucleotide of the cis-reg sequence to be deleted*/
    temp=&genotype->cisreg_seq[which_gene][0];   
    /* shift the cisreg array to overwrite the cis_seq to be deleted */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-1-which_gene);i++) 
    {				
        *temp=*(temp+CISREG_LEN);         
        temp++;				
    }    
    protein_id=genotype->which_protein[which_gene];      
    /*if this tf has only one copy of gene, then we'll delete the binding seq of this tf 
     * and remove the tf from locus_specific_TF_behavior*/ 
    if(genotype->protein_pool[protein_id][0][0]==1) 
    {            
        /* shift the tf_reg array to overwrite the binding sequence to be deleted */
        temp=&genotype->tf_binding_seq[protein_id][0];     
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-1-protein_id);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        }
        /* shift the tf_reg_rc array to overwrite the binding sequence to be deleted */    
        temp=&genotype->tf_binding_seq_rc[protein_id][0];
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-1-protein_id);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        } 
        /*UPDATE locus_specific_TF_behavior*/
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
        {
            for(j=protein_id;j<genotype->nproteins;j++)
                genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i][j+1];
        }    
    }     
    /*Before removing info of which_gene*/
    /*update total loci length*/
    genotype->total_loci_length-=genotype->locus_length[which_gene]; 
    /*update output_gene_ids*/   
    update_output_protein_pool(genotype,which_gene,'w');  
    /*removing info of which_gene*/
    for(i=which_gene;i<genotype->ngenes;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];            
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->protein_syn_rate[i]=genotype->protein_syn_rate[i+1]; 
        genotype->locus_length[i]=genotype->locus_length[i+1];
        genotype->is_output[i]=genotype->is_output[i+1];
        for(j=0;j<MAX_PROTEINS;j++)
            genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i+1][j];
        genotype->recalc_TFBS[i]=1;
    }    
    /*update others*/
    update_node_family_pool(genotype,which_gene,'w'); //be sure to update node_family_pool bf update ngenes and nproteins
    update_protein_pool(genotype,protein_id,which_gene,'w');     
    update_cisreg_cluster_pool(genotype,which_gene,'w');// ngenes is reduced in update_cisreg_cluster_pool 
}

void reproduce_whole_gene_deletion(Genotype *genotype, Mutation *mut_record) // any gene can be deleted
{    
    int which_gene, i,j;
    char *temp;    
    int protein_id; 
    which_gene=mut_record->which_gene;     
    temp = &genotype->cisreg_seq[which_gene][0];	
    for(i=0;i<CISREG_LEN*(genotype->ngenes-which_gene-1);i++) 
    {				
        *temp=*(temp+CISREG_LEN);         
        temp++;				
    }
    protein_id=genotype->which_protein[which_gene];
    if(genotype->protein_pool[protein_id][0][0]==1) 
    {  
        temp=&genotype->tf_binding_seq[protein_id][0];       
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        }            
        temp=&genotype->tf_binding_seq_rc[protein_id][0];
        for(i=0;i<TF_ELEMENT_LEN*(genotype->nproteins-protein_id-1);i++)
        {
            *temp=*(temp+TF_ELEMENT_LEN);
            temp++;
        }         
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
        {
            for(j=protein_id;j<genotype->nproteins;j++)
                genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i][j+1];
        }    
    }
    /*Before removing info of which_gene*/
    /*update total loci length*/
    genotype->total_loci_length-=genotype->locus_length[which_gene]; 
    /*update output_gene_ids*/  
    update_output_protein_pool(genotype,which_gene,'w');  
    /*removing info of which_gene*/
    for(i=which_gene;i<genotype->ngenes;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];            
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->protein_syn_rate[i]=genotype->protein_syn_rate[i+1]; 
        genotype->locus_length[i]=genotype->locus_length[i+1];
        genotype->is_output[i]=genotype->is_output[i+1];
        for(j=0;j<MAX_PROTEINS;j++)
            genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i+1][j];
        genotype->recalc_TFBS[i]=1;
    }    
    /*update others*/
    update_node_family_pool(genotype,which_gene,'w');
    update_protein_pool(genotype,protein_id,which_gene,'w');     
    update_cisreg_cluster_pool(genotype,which_gene,'w');   
}

/*
 * Duplicate a whole cis-reg sequence
 */
void mut_duplication(Genotype *genotype, Mutation *mut_record, RngStream RS) 
{
    int which_gene, i, protein_id, which_node;
    char *temp1, *temp2;   
    
    while(1)
    {
        which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
        which_node=genotype->which_node_family[which_gene];        
        if(genotype->n_output_genes<MAX_OUTPUT_GENES) // not too many effector gene copies
        {
            if(genotype->node_family_pool[which_node][0][0]<MAX_COPIES_PER_NON_OUTPUT_GENE) // not too many copies
                break;            
        }
        else
        {
            if(genotype->node_family_pool[which_node][0][0]<MAX_COPIES_PER_NON_OUTPUT_GENE && genotype->is_output[which_gene]==NON_OUTPUT_PROTEIN) // not too many non-output gene
                break;
        }
    }     
    /*record mutation info*/
    mut_record->which_gene=which_gene; 
    /* copy the promoter*/
    temp1=&genotype->cisreg_seq[which_gene][0]; /* points to the gene to be duplicated*/
    temp2=&genotype->cisreg_seq[genotype->ngenes][0]; /* points to the end of the effector gene */
    /* shift the sequences of the effector gene CISREG_LEN bp */
    for(i=0;i<CISREG_LEN;i++)
        *temp2++=*temp1++;           
    /* copy and paste info to the slot*/
    genotype->min_N_activator_to_transc[genotype->ngenes]=genotype->min_N_activator_to_transc[which_gene];
    genotype->active_to_intermediate_rate[genotype->ngenes]=genotype->active_to_intermediate_rate[which_gene];
    genotype->mRNA_decay_rate[genotype->ngenes]=genotype->mRNA_decay_rate[which_gene];
    genotype->protein_decay_rate[genotype->ngenes]=genotype->protein_decay_rate[which_gene];
    genotype->protein_syn_rate[genotype->ngenes]=genotype->protein_syn_rate[which_gene]; 
    genotype->locus_length[genotype->ngenes]=genotype->locus_length[which_gene];
    genotype->is_output[genotype->ngenes]=genotype->is_output[which_gene];
    for(i=0;i<MAX_PROTEINS;i++)
        genotype->locus_specific_TF_behavior[genotype->ngenes][i]=genotype->locus_specific_TF_behavior[which_gene][i];
    genotype->recalc_TFBS[genotype->ngenes]=YES;
    /*update total loci length*/
    genotype->total_loci_length+=genotype->locus_length[which_gene];
    /*update pools*/
    if(genotype->is_output[genotype->ngenes]==OUTPUT_PROTEIN)       
        update_output_protein_pool(genotype,genotype->ngenes,'d');  
    protein_id=genotype->which_protein[which_gene];   
    update_node_family_pool(genotype,which_gene,'d');
    update_protein_pool(genotype,protein_id,which_gene,'d'); 
    /* update cisreg_cluster_pool*/    
    update_cisreg_cluster_pool(genotype,which_gene,'d');
    /* update gene numbers*/  
    genotype->ngenes++;
}

void reproduce_gene_duplication(Genotype *genotype, Mutation *mut_record) //any gene can be duplicated
{   
    int which_gene, i, protein_id;
    char *temp1, *temp2; 
    /*get the gene to be duplicated from record*/
    which_gene=mut_record->which_gene;      
    temp1=&genotype->cisreg_seq[which_gene][0]; 
    temp2=&genotype->cisreg_seq[genotype->ngenes][0];  
    for(i=0;i<CISREG_LEN;i++) 
        *temp2++=*temp1++;    
    genotype->min_N_activator_to_transc[genotype->ngenes]=genotype->min_N_activator_to_transc[which_gene];
    genotype->active_to_intermediate_rate[genotype->ngenes]=genotype->active_to_intermediate_rate[which_gene];
    genotype->mRNA_decay_rate[genotype->ngenes]=genotype->mRNA_decay_rate[which_gene];
    genotype->protein_decay_rate[genotype->ngenes]=genotype->protein_decay_rate[which_gene];
    genotype->protein_syn_rate[genotype->ngenes]=genotype->protein_syn_rate[which_gene]; 
    genotype->locus_length[genotype->ngenes]=genotype->locus_length[which_gene];
    genotype->is_output[genotype->ngenes]=genotype->is_output[which_gene];   
    genotype->total_loci_length+=genotype->locus_length[which_gene];
    for(i=0;i<MAX_PROTEINS;i++)
        genotype->locus_specific_TF_behavior[genotype->ngenes][i]=genotype->locus_specific_TF_behavior[which_gene][i];    
    genotype->recalc_TFBS[genotype->ngenes]=YES;
    protein_id=genotype->which_protein[which_gene];
    if(genotype->is_output[genotype->ngenes]==OUTPUT_PROTEIN)
        update_output_protein_pool(genotype,which_gene,'d'); 
    update_node_family_pool(genotype,which_gene,'d');    
    update_protein_pool(genotype,protein_id,which_gene,'d'); 
    update_cisreg_cluster_pool(genotype,which_gene,'d');      
    genotype->ngenes++;     
}

/*
 * Mutation to the binding sequence of a TF gene
 */
void mut_binding_sequence(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene, which_nucleotide, protein_id, i;
    char nucleotide;      
    char *tf_binding_seq, *tf_binding_seq_rc,*temp1,*temp2;
    /*get which gene to mutate*/
    if(genotype->flag_effector_is_TF)
        which_gene=RngStream_RandInt(RS,0,genotype->ngenes-1); 
    else
    {
        do 
            which_gene=RngStream_RandInt(RS,0,genotype->ngenes-1);
        while(genotype->is_output[which_gene]==OUTPUT_PROTEIN);
    }
    /*if this TF has more than one copies of gene, then the mutation adds a new tf
     *which requires a new slot in tf_binding_seq and tf_binding_seq_rc to store the new binding seq*/
    protein_id=genotype->which_protein[which_gene];    
    if(genotype->protein_pool[protein_id][0][0]>1)
    {
        /*points to the first nucleotide of the binding sequence*/
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        /*points to an empty slot*/
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        /*copy the binding sequences to empty slots*/
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }
        /*point tf_binding_seq and tf_binding_seq_rc to the new slot so that we can apply mutation later*/
        tf_binding_seq=&genotype->tf_binding_seq[genotype->nproteins][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        /*update node_family_pool*/
        update_node_family_pool(genotype,which_gene,'c');
        /*update protein pool*/
        update_protein_pool(genotype,protein_id,which_gene,'c');    
        /* make Kd for the new protein*/
        genotype->Kd[genotype->nproteins-1]=genotype->Kd[protein_id]; //nproteins has been increased by 1 in update_protein_pool
        /* update activator or repressor numbers, and protein_identity*/
        genotype->protein_identity[genotype->nproteins-1]=genotype->protein_identity[protein_id];
        if(genotype->protein_identity[genotype->nproteins-1]==ACTIVATOR) //mutation to binding seq does not change the identity of a tf
            genotype->N_act++;
        else
            genotype->N_rep++;             
        /* Update locus_specific_TF_behavior: assuming mutation to binding sequence changes nothing to the locus-specific behavior */   
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else /*if this tf has only one copy of gene, no new slot is required*/
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        /*update node_family_pool*/
        update_node_family_pool(genotype,which_gene,'c');
    }
   
    /*mutate the binding sequence, only changing one nucleotide in the binding sequence*/
    which_nucleotide=RngStream_RandInt(RS,0,TF_ELEMENT_LEN-1);
    nucleotide=set_base_pair(RngStream_RandU01(RS));        
    while (nucleotide == tf_binding_seq[which_nucleotide])
        nucleotide=set_base_pair(RngStream_RandU01(RS));
    tf_binding_seq[which_nucleotide]=nucleotide;    
    /*record mutation info*/
    mut_record->which_gene=which_gene;
    mut_record->which_nucleotide=which_nucleotide;
    mut_record->nuc_diff[0]=nucleotide;    
    /* update the reverse complement sequence*/
    switch (nucleotide)
    {
        case 'g':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; 
    }  
    /*decide whether to update cisreg clusters. Mutation to binding seq may differentiate bs distributions among genes in a cluster*/     
    update_cisreg_cluster_pool(genotype,which_gene,'c');  
    /* Remember to calculate TFBSs if this mutant is accepted */
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;  
}

void reproduce_mut_binding_sequence(Genotype *genotype, Mutation *mut_record)
{    
    int which_gene, which_nucleotide, protein_id, i;    
    char *tf_binding_seq, *tf_binding_seq_rc, *temp1,*temp2;     
    which_gene=mut_record->which_gene;
    protein_id=genotype->which_protein[which_gene];  
    if(genotype->protein_pool[protein_id][0][0]>1)
    {        
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }      
        tf_binding_seq=&genotype->tf_binding_seq[genotype->nproteins][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[genotype->nproteins][0];     
        update_node_family_pool(genotype,which_gene,'c');
        update_protein_pool(genotype,protein_id,which_gene,'c');  
        genotype->Kd[genotype->nproteins-1]=genotype->Kd[protein_id];
        genotype->protein_identity[genotype->nproteins-1]=genotype->protein_identity[protein_id];
        if(genotype->protein_identity[genotype->nproteins-1]==ACTIVATOR) //mutation to binding seq does not change the identity of a tf
            genotype->N_act++;
        else
            genotype->N_rep++;     
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else 
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        update_node_family_pool(genotype,which_gene,'c');
    }     
    which_nucleotide=mut_record->which_nucleotide;        
    tf_binding_seq[which_nucleotide]=mut_record->nuc_diff[1];
    switch (tf_binding_seq[which_nucleotide])
    {
        case 'g':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='c'; break;
        case 'c':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='g'; break;
        case 'a':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='t'; break;
        case 't':
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
    }     
    update_cisreg_cluster_pool(genotype,which_gene,'c');    
    for(i=0;i<genotype->ngenes;i++)
        genotype->recalc_TFBS[i]=YES;   
}

/* Mutations to the rate of mRNA_decay, translation, protein_decay, and pic_disassembly 
 * will be mutated. 
 */
void mut_kinetic_constant(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    float random;
    int random2;    
    int which_gene;
    float total_mut_flux; 
    float total_mut_ACT_to_INT; 
    float total_mut_mRNA_decay;
    float total_mut_protein_syn_rate;
    float total_mut_protein_decay;
    float total_mut_cooperation;
    
    total_mut_ACT_to_INT=(genotype->ngenes-N_SIGNAL_TF)*MUT_ACT_to_INT;
    total_mut_mRNA_decay=genotype->total_loci_length*MUT_mRNA_decay;
    total_mut_protein_decay=genotype->total_loci_length*MUT_protein_decay;
    total_mut_protein_syn_rate=genotype->total_loci_length*MUT_protein_syn_rate;
    total_mut_cooperation=(genotype->ngenes-genotype->ngenes)*MUT_cooperation;
    total_mut_flux=total_mut_ACT_to_INT+total_mut_mRNA_decay+total_mut_protein_decay+total_mut_protein_syn_rate+total_mut_cooperation;
   
    while(1)
    {
        random=RngStream_RandU01(RS)*total_mut_flux;
        /********************************** mut kdis *******************************/
        if(random<=total_mut_ACT_to_INT) 
        {  
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);            
            genotype->active_to_intermediate_rate[which_gene]=mut_make_new_value(genotype->active_to_intermediate_rate[which_gene],
                                                                                miu_ACT_TO_INT_RATE,
                                                                                sigma_ACT_TO_INT_RATE,
                                                                                MAX_ACT_TO_INT_RATE,
                                                                                MIN_ACT_TO_INT_RATE,
                                                                                RS,
                                                                                mut_record,
                                                                                BOUND_INCLUDED);
            /*record mutation info*/
            mut_record->which_gene=which_gene;
            mut_record->kinetic_type=0;
            mut_record->kinetic_diff=genotype->active_to_intermediate_rate[which_gene];
            break;
        }
        else 
        {
            random-=total_mut_ACT_to_INT;
            /****************** mut cooperation ***********************/
            if(random<=total_mut_cooperation && MUT_cooperation!=0.0) /// this part is not correct
            {   
                which_gene=0;
                while(genotype->which_protein[which_gene]!=genotype->nproteins-1)
                    which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                mut_record->which_gene=which_gene;
                if(genotype->min_N_activator_to_transc[which_gene]==1)
                    genotype->min_N_activator_to_transc[which_gene]=2;
                else
                    genotype->min_N_activator_to_transc[which_gene]=1;
                mut_record->kinetic_type=4;
                mut_record->kinetic_diff=(float)genotype->min_N_activator_to_transc[which_gene];
                if(genotype->cisreg_cluster_pool[genotype->which_cluster[which_gene]][0][0]>1)
                    update_cisreg_cluster_pool(genotype,which_gene,'s'); 
                break;
            }
            else
            { 
                /* The rest mutation rates are proportional to gene length, so we can pick which genes first*/
                random2=RngStream_RandInt(RS,1,genotype->total_loci_length);
                which_gene=N_SIGNAL_TF-1;
                while(which_gene<genotype->ngenes-1 && random2>0)
                {
                    which_gene++;
                    random2-=genotype->locus_length[which_gene];
                } 
                mut_record->which_gene=which_gene;
                /* now pick which type of mutation*/
                random-=total_mut_cooperation;
                /******************************** mut mRNAdecay **********************/
                if(random<=total_mut_mRNA_decay) 
                {       
                    genotype->mRNA_decay_rate[which_gene]=mut_make_new_value(genotype->mRNA_decay_rate[which_gene],miu_mRNA_decay,sigma_mRNA_decay, MAX_MRNA_DECAY, MIN_MRNA_DECAY,RS,mut_record,BOUND_INCLUDED);
                    /*record mutation info*/
                    mut_record->kinetic_type=1;
                    mut_record->kinetic_diff=genotype->mRNA_decay_rate[which_gene];
                    break;
                }
                else 
                {
                    random-=total_mut_mRNA_decay;
                    /*************************** mut translation *********************/
                    if(random<=total_mut_protein_syn_rate) 
                    {       
                        genotype->protein_syn_rate[which_gene]=mut_make_new_value(genotype->protein_syn_rate[which_gene],miu_protein_syn_rate,sigma_protein_syn_rate, MAX_PROTEIN_SYN_RATE, MIN_PROTEIN_SYN_RATE,RS,mut_record,BOUND_INCLUDED);        
                        mut_record->kinetic_type=2;
                        mut_record->kinetic_diff=genotype->protein_syn_rate[which_gene];
                        break;
                    }
                    else 
                    {  
                        random-=total_mut_protein_syn_rate;
                        /********************* mut protein decay **********************/
                        if(random<=total_mut_protein_decay)
                        {                   
                            genotype->protein_decay_rate[which_gene]=mut_make_new_value(genotype->protein_decay_rate[which_gene],miu_protein_decay,sigma_protein_decay,MAX_PROTEIN_DECAY,MIN_PROTEIN_DECAY,RS,mut_record,BOUND_INCLUDED);                  
                            mut_record->kinetic_type=3;
                            mut_record->kinetic_diff=genotype->protein_decay_rate[which_gene];
                            break;
                        }                    
                    } 
                }
            }
        }
    }
}

void reproduce_mut_kinetic_constant(Genotype *genotype, Mutation *mut_record)
{    
    int which_gene;            
    which_gene=mut_record->which_gene; /* which gene */ 
    switch (mut_record->kinetic_type)
    {
        case 0: /* mut r_act_to_int */
            genotype->active_to_intermediate_rate[which_gene]=mut_record->kinetic_diff;
            break;
        case 1: /* mut mRNAdecay */ 
            genotype->mRNA_decay_rate[which_gene]=mut_record->kinetic_diff;
            break;
        case 2: /* mut translation */                      
            genotype->protein_syn_rate[which_gene]=mut_record->kinetic_diff;           
            break;        
        case 3: /* mut protein decay */
            genotype->protein_decay_rate[which_gene]=mut_record->kinetic_diff;            
            break;
        case 4: /* mut cooperation */
            genotype->min_N_activator_to_transc[which_gene]=(int)mut_record->kinetic_diff;
            if(genotype->cisreg_cluster_pool[genotype->which_cluster[which_gene]][0][0]>1)
                update_cisreg_cluster_pool(genotype,which_gene,'s'); 
    }
}

void mut_locus_length(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene;
    int random;    
   
    random=RngStream_RandInt(RS,1,genotype->total_loci_length);
    which_gene=N_SIGNAL_TF-1;
    while(which_gene<genotype->ngenes-1 && random>0)
    {
        which_gene++;
        random-=genotype->locus_length[which_gene];
    }        
    
    genotype->total_loci_length-=genotype->locus_length[which_gene];
    
    if(genotype->locus_length[which_gene]==MAX_GENE_LENGTH)
        genotype->locus_length[which_gene]-=1;
    else if(genotype->locus_length[which_gene]==MIN_GENE_LENGTH)
        genotype->locus_length[which_gene]+=1;
    else                            
        genotype->locus_length[which_gene]+=(RngStream_RandU01(RS)<0.5)?1:-1;  
    
    genotype->total_loci_length+=genotype->locus_length[which_gene];  
    mut_record->mut_type='l';   
    mut_record->which_gene=which_gene;
    mut_record->kinetic_diff=(float)genotype->locus_length[which_gene];
}

void reproduce_mut_locus_length(Genotype *genotype, Mutation *mut_record)
{
    genotype->total_loci_length-=genotype->locus_length[mut_record->which_gene];   
    genotype->locus_length[mut_record->which_gene]=(int)mut_record->kinetic_diff;
    genotype->total_loci_length+=(int)genotype->locus_length[mut_record->which_gene]; 
}

/*activator to repressor or vice versa*/
void mut_identity(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int tf_id,protein_id,i;   
    char *tf_binding_seq,*tf_binding_seq_rc,*temp1,*temp2; 
    /*which tf gene to mutate*/
    if(genotype->flag_effector_is_TF)
        tf_id = RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1); // the first sensor tf must be an activator, therefore is not subject to mutation
    else
    {
        do 
            tf_id=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
        while(genotype->is_output[tf_id]==OUTPUT_PROTEIN);
    }
    protein_id=genotype->which_protein[tf_id];  
    /*save record*/
    mut_record->which_gene=tf_id;       
    /* if this tf gene has more than one copies, the mutation adds a new protein*/
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        /*give the new protein its own binding sequence*/
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0]; 
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        } 
        /*update node_family_pool*/
        update_node_family_pool(genotype,tf_id,'e');
        update_protein_pool(genotype,protein_id,tf_id,'e');  
        /* update Kd*/
        genotype->Kd[genotype->nproteins-1]=genotype->Kd[protein_id]; //nproteins has been updated in update_protein_pool
        /* update_protein_pool put the new protein at nproteins and then increases nproteins by 1, 
        * so the new protein is at nproteins-1 now. Note that N_act and N_rep is updated in update_protein_pool*/ 
        genotype->protein_identity[genotype->nproteins-1]=(genotype->protein_identity[protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR;   
        /* update protein_identity*/
        if(genotype->protein_identity[genotype->nproteins-1]==REPRESSOR) 
            genotype->N_rep++;  /* an activator turns into a repressor */
        else
            genotype->N_act++;
        /*the locus specific behavior is the identity of the mutated tf*/
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->protein_identity[genotype->nproteins-1];
    }
    else //otherwise we just flip the property of an existing TF
    {
        genotype->protein_identity[protein_id]=(genotype->protein_identity[protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id]==ACTIVATOR)//mutate to activator
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else //otherwise to repressor
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][protein_id]=genotype->protein_identity[protein_id];
        /*update node_family_pool*/
        update_node_family_pool(genotype,tf_id,'e');
    } 
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES; /* recalculate binding sites on every promoter */        
}

void reproduce_mut_identity(Genotype *genotype, Mutation *mut_record)
{
    int tf_id, protein_id,i;
    char *tf_binding_seq,*tf_binding_seq_rc,*temp1,*temp2;    
    tf_id = mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }   
        update_node_family_pool(genotype,tf_id,'e');
        update_protein_pool(genotype,protein_id,tf_id,'e');         
        genotype->Kd[genotype->nproteins-1]=genotype->Kd[protein_id];       
        genotype->protein_identity[genotype->nproteins-1]=(genotype->protein_identity[protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[genotype->nproteins-1]==REPRESSOR) 
            genotype->N_rep++;
        else
            genotype->N_act++;
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->protein_identity[genotype->nproteins-1];
    }
    else
    {
        genotype->protein_identity[protein_id]=(genotype->protein_identity[protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id]==ACTIVATOR)
        {
            genotype->N_act++;
            genotype->N_rep--;
        }
        else
        { 
            genotype->N_rep++;
            genotype->N_act--;
        }
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][protein_id]=genotype->protein_identity[protein_id];
        update_node_family_pool(genotype,tf_id,'e');
    } 
    
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=YES;   
}

/*
 * mutate affinity of TF
 */
void mut_Kd(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int i;   
    float new_Kd; 
    int tf_id,protein_id;
    char *tf_binding_seq,*tf_binding_seq_rc,*temp1,*temp2;
    /*which TF to mutate*/
    if(genotype->flag_effector_is_TF)
        tf_id=RngStream_RandInt(RS,0,genotype->ngenes-1); 
    else
    {
        do 
            tf_id=RngStream_RandInt(RS,0,genotype->ngenes-1);
        while(genotype->is_output[tf_id]==OUTPUT_PROTEIN);
    }
    protein_id=genotype->which_protein[tf_id];   
    /*generate a new koff */    
    new_Kd=mut_make_new_value(genotype->Kd[protein_id],miu_Kd,sigma_Kd,MAX_KD,MIN_KD,RS,mut_record,BOUND_EXCLUDED);   
    /* if this tf gene has more than one copies, the mutation adds a new protein*/    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }
//        /*update node_family_pool*/
//        update_node_family_pool(genotype,tf_id,'f');
        /*update protein pool*/
        update_protein_pool(genotype,protein_id,tf_id,'f');
        /* update protein_identity. update_protein_pool put the new protein at nproteins 
         * and then increases nproteins by 1, so the new protein is at nproteins-1 now.*/
        genotype->protein_identity[genotype->nproteins-1]=genotype->protein_identity[protein_id];       
        if(genotype->protein_identity[genotype->nproteins-1]==ACTIVATOR) 
            genotype->N_act++;  
        else
            genotype->N_rep++;   
        genotype->Kd[genotype->nproteins-1]=new_Kd;  
        /* Update locus_specific_TF_behavior: assuming mutation to Kd changes nothing to the locus-specific behavior */   
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }                                                                                                           
    else
    {
        genotype->Kd[protein_id]=new_Kd;  
//        update_node_family_pool(genotype,tf_id,'f');
    }
    /*record mutation*/
    mut_record->which_gene=tf_id;
    mut_record->kinetic_diff=new_Kd; 
    /* recalculate binding sites on every promoter */
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;
}

void reproduce_mut_Kd(Genotype *genotype, Mutation *mut_record)
{   
    int tf_id,protein_id,i;
    char *tf_binding_seq, *tf_binding_seq_rc, *temp1, *temp2;    
    tf_id=mut_record->which_gene;
    protein_id=genotype->which_protein[tf_id];    
    if(genotype->protein_pool[protein_id][0][0]!=1)
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
        temp1=&genotype->tf_binding_seq[genotype->nproteins][0];
        temp2=&genotype->tf_binding_seq_rc[genotype->nproteins][0];
        for(i=0;i<TF_ELEMENT_LEN;i++)
        {
            *temp1++=*tf_binding_seq++;
            *temp2++=*tf_binding_seq_rc++;
        }
        update_protein_pool(genotype,protein_id,tf_id,'f');  
        genotype->protein_identity[genotype->nproteins-1]=genotype->protein_identity[protein_id];       
        if(genotype->protein_identity[genotype->nproteins-1]==ACTIVATOR) 
            genotype->N_act++;  
        else
            genotype->N_rep++;  
        genotype->Kd[genotype->nproteins-1]=mut_record->kinetic_diff;  
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else    
        genotype->Kd[protein_id]=mut_record->kinetic_diff; 
    for(i=0;i<genotype->ngenes;i++) 
        genotype->recalc_TFBS[i]=YES;   
}

/*mutate an effector gene to normal tf gene*/
void mut_effector2TF(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_gene;
    while(1)
    {
        which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
        if(genotype->is_output[which_gene]==OUTPUT_PROTEIN)
            break;
    }
    mut_record->which_gene=which_gene; 
    update_output_protein_pool(genotype,which_gene,'t');       
    update_node_family_pool(genotype,which_gene,'t');
}

void reproduce_effector2TF(Genotype *genotype, Mutation *mut_record)
{         
    update_output_protein_pool(genotype,mut_record->which_gene,'t');  
    update_node_family_pool(genotype,mut_record->which_gene,'t');
}

/* mut gene specific tf behavior*/
void mut_locus_specific_tf_behavior(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int i, gene_id, protein_ids[genotype->nproteins], which_cluster, which_node;
    /*which locus to mutate*/
    gene_id=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    /*which node is flipped*/
    if(genotype->flag_effector_is_TF)
        which_node=RngStream_RandInt(RS,0,genotype->N_node_families-1);
    else
    {
        do 
            which_node=RngStream_RandInt(RS,0,genotype->N_node_families-1);
        while(genotype->is_output[genotype->node_family_pool[which_node][1][0]]==OUTPUT_PROTEIN);        
    }
    for(i=0;i<genotype->nproteins;i++)
        protein_ids[i]=NA;
    for(i=0;i<genotype->node_family_pool[which_node][0][0];i++)
        protein_ids[genotype->which_protein[genotype->node_family_pool[which_node][1][i]]]=1;
    /*flip the TFs in the node*/
    for(i=0;i<genotype->nproteins;i++)
    {
        if(protein_ids[i]==1)
            genotype->locus_specific_TF_behavior[gene_id][i]=(genotype->locus_specific_TF_behavior[gene_id][i]==ACTIVATOR)?REPRESSOR:ACTIVATOR;
    }
    
    which_cluster=genotype->which_cluster[gene_id];    
    /*if the cis-regulatory sequence of gene_id is not unique, the mutation kicks gene_id out of its cisreg_cluster_pool*/
    if(genotype->cisreg_cluster_pool[which_cluster][0][0]>1)
        update_cisreg_cluster_pool(genotype,gene_id,'s');    
    mut_record->which_gene=gene_id;
    mut_record->which_protein=which_node;   
    genotype->recalc_TFBS[gene_id]=YES;    
}

void reproduce_mut_locus_specific_tf_behavior(Genotype *genotype, Mutation *mut_record)
{
    int i, protein_ids[genotype->nproteins],which_node;
    which_node=mut_record->which_protein;
    for(i=0;i<genotype->nproteins;i++)
        protein_ids[i]=NA;
    for(i=0;i<genotype->node_family_pool[which_node][0][0];i++)
        protein_ids[genotype->which_protein[genotype->node_family_pool[which_node][1][i]]]=1;
    for(i=0;i<genotype->nproteins;i++)
    {
        if(protein_ids[i]==1)
            genotype->locus_specific_TF_behavior[mut_record->which_gene][i]=(genotype->locus_specific_TF_behavior[mut_record->which_gene][i]==ACTIVATOR)?REPRESSOR:ACTIVATOR;
    }
    if(genotype->cisreg_cluster_pool[genotype->which_cluster[mut_record->which_gene]][0][0]>1)
        update_cisreg_cluster_pool(genotype,mut_record->which_gene,'s');    
}

/*
 *Wrapper of all mutation functions
 */
void mutate(Genotype *genotype, RngStream RS, Mutation *mut_record)
{
    int i;   
    draw_mutation(genotype, &(mut_record->mut_type),RS);
    /*reset records of nucleotides*/
    for(i=0;i<3;i++)
        mut_record->nuc_diff[i]='\0';
    switch (mut_record->mut_type)
    {
        case 's': //substitution in cis-reg       		
            mut_substitution(genotype,mut_record,RS);
            break;      
        case 'w': // whole gene deletion.          
            mut_whole_gene_deletion(genotype,mut_record,RS);            
            break;
        case 'd': // Whole gene duplication                  
            mut_duplication(genotype,mut_record,RS);            
            break;
        case 'c': //binding sequence 
            mut_binding_sequence(genotype,mut_record,RS);            
            break;
        case 'k': //mutations to kinetic constants        
            mut_kinetic_constant(genotype, mut_record,RS);           
            break;
        case 'e': //activator to repressor or the reverse
            mut_identity(genotype, mut_record, RS);            
            break;
        case 'f': //mutations to the Kd of a tf
            mut_Kd(genotype,mut_record,RS);            
            break;
        case 'l': //mutate locus length
            mut_locus_length(genotype,mut_record, RS);
            break;
        case 't': //effector to regular TF
            mut_effector2TF(genotype,mut_record,RS);
            break;
        case 'o':
            mut_locus_specific_tf_behavior(genotype, mut_record, RS);
            break;
    }
}

/* this function perform mutation indicated by input. Used to reproduce the genotype following a serial of mutations*/
void reproduce_mutate(Genotype *genotype, Mutation *mut_record)
{ 
    switch (mut_record->mut_type)
    {
        case 's': //substitution  
            reproduce_substitution(genotype,mut_record);
            break;        
        case 'w': // whole gene deletion.  
            reproduce_whole_gene_deletion(genotype,mut_record);
            break;
        case 'd': // Whole gene duplication                 
            reproduce_gene_duplication(genotype,mut_record);
            break;        
        case 'c': //binding sequence        
            reproduce_mut_binding_sequence(genotype,mut_record);
            break;
        case 'k': //mutations in kinetic constants        
            reproduce_mut_kinetic_constant(genotype, mut_record);            
            break;
        case 'e': //changing the identity of a TF
            reproduce_mut_identity(genotype, mut_record);
            break;
        case 'f':
            reproduce_mut_Kd(genotype,mut_record);
            break;  
        case 'l':
            reproduce_mut_locus_length(genotype, mut_record);
            break;
        case 't': //effector to regular TF
            reproduce_effector2TF(genotype,mut_record);
            break;
        case 'o':
            reproduce_mut_locus_specific_tf_behavior(genotype, mut_record);
            break;
    }   
}

/* this function calculates the probability of different mutations based on
 * the current genotype, and draws a mutation accordingly.
 */
void draw_mutation(Genotype *genotype, char *mut_type, RngStream RS)
{
    float random;
    float tot_mut_rate=0.0;
    float tot_subs_rate, tot_dup_rate, tot_sil_rate, tot_mut_kin_rate, tot_mut_identity_rate, tot_mut_binding_seq_rate, tot_mut_koff_rate; 
    float tot_effector2TF, tot_mut_locus_length_rate, tot_mut_locus_specific_tf_behavior;
    int N_target_genes,N_effector_protein;
    int i,j;    
    /* duplication rate*/  
    tot_dup_rate=0.0;
    N_target_genes=0;
    if(genotype->ngenes<MAX_GENES-1) // not too many genes
    {
        if(genotype->n_output_genes<MAX_OUTPUT_GENES) // not too many output genes
        {            
            for(i=N_SIGNAL_TF;i<genotype->N_node_families;i++)
            {
                if(genotype->node_family_pool[i][0][0]<MAX_COPIES_PER_NON_OUTPUT_GENE)//not too many non-output genes
                    N_target_genes+=genotype->node_family_pool[i][0][0];
            }
        }
        else
        {
            for(i=N_SIGNAL_TF;i<genotype->N_node_families;i++)
            {
                if(genotype->node_family_pool[i][0][0]<MAX_COPIES_PER_NON_OUTPUT_GENE && genotype->is_output[genotype->node_family_pool[i][1][0]]!=OUTPUT_PROTEIN)//not too many non-output genes   
                    N_target_genes+=genotype->node_family_pool[i][0][0];
            }           
        }
    }         
    tot_dup_rate=N_target_genes*DUPLICATION;
    tot_mut_rate+=tot_dup_rate;     
    /* silencing rate*/ 
    N_target_genes=0;
    if(genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF>MIN_NON_OUTPUT_GENES)//if there's more than MIN_NON_OUTPUT_GENES
        N_target_genes+=genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF;
    if(genotype->n_output_genes>1)
        N_target_genes+=genotype->n_output_genes; 
    tot_sil_rate=N_target_genes*SILENCING; 
    tot_mut_rate+=tot_sil_rate;    
    
    /* calc total susbtitution rate*/
    tot_subs_rate=(genotype->ngenes-N_SIGNAL_TF)*CISREG_LEN*SUBSTITUTION; //NA to the sensor TF gene
    tot_mut_rate+=tot_subs_rate;  
    
    /* mut in kinetic constants */    
    tot_mut_kin_rate=(genotype->ngenes-N_SIGNAL_TF)*MUT_ACT_to_INT; // NA to the sensor TF gene
    tot_mut_kin_rate+=genotype->total_loci_length*(MUT_protein_syn_rate+MUT_mRNA_decay+MUT_protein_decay);
    tot_mut_kin_rate+=(genotype->ngenes-N_SIGNAL_TF)*MUT_cooperation;
    tot_mut_rate+=tot_mut_kin_rate;
    
    /* mut in binding seq*/
    tot_mut_binding_seq_rate=genotype->ngenes*MUT_binding_seq; // NA to the effector genes
    if(!(genotype->flag_effector_is_TF))
        tot_mut_binding_seq_rate-=genotype->n_output_genes*MUT_binding_seq;
    tot_mut_rate+=tot_mut_binding_seq_rate;    
    
    /* mut in identity*/
    tot_mut_identity_rate=(genotype->ngenes-N_SIGNAL_TF)*MUT_identity; // NA to the sensor TF gene
    if(!(genotype->flag_effector_is_TF))
        tot_mut_identity_rate-=genotype->n_output_genes*MUT_identity;
    tot_mut_rate+=tot_mut_identity_rate;
    
    /* mut in Kd*/
    tot_mut_koff_rate=genotype->ngenes*MUT_Kd;  
    if(!(genotype->flag_effector_is_TF))
        tot_mut_koff_rate-=genotype->n_output_genes*MUT_Kd;
    tot_mut_rate+=tot_mut_koff_rate;    
    
    /* mut locus length*/
    tot_mut_locus_length_rate=MUT_GENE_LENGTH_RATE*genotype->total_loci_length;
    tot_mut_rate+=tot_mut_locus_length_rate;
    
    /*effector to regular TF*/
    if(genotype->flag_effector_is_TF)
    {
        if(genotype->n_output_genes>1)
            tot_effector2TF=genotype->n_output_genes*MUT_effector_to_TF;
        else
            tot_effector2TF=0.0;
    }
    else
        tot_effector2TF=0.0;
    tot_mut_rate+=tot_effector2TF;
    
    /*gene-specific behavior of TF*/
    tot_mut_locus_specific_tf_behavior=(genotype->ngenes-N_SIGNAL_TF)*genotype->nproteins*MUT_locus_specific_tf_behavior;
    if(!(genotype->flag_effector_is_TF))
    {
        N_effector_protein=0;
        for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
        {
            if(genotype->is_output[genotype->protein_pool[i][1][0]]==OUTPUT_PROTEIN)
                N_effector_protein++;
        }
        tot_mut_locus_specific_tf_behavior-=N_effector_protein*(genotype->ngenes-N_SIGNAL_TF)*MUT_locus_specific_tf_behavior;
    }
    tot_mut_rate+=tot_mut_locus_specific_tf_behavior;        
        
    /*Draw a mutation based on the above rates*/
    while(1)
    {
        random=RngStream_RandU01(RS)*tot_mut_rate;    
        if(random<=tot_subs_rate) 
        {
            *mut_type='s';                      /* substitution*/
            break;
        }
        else
        {   
            random-=tot_subs_rate;
            if(random<= tot_dup_rate)
            {
                *mut_type='d';                          /* gene duplication */
                break;
            }
            else
            {
                random-=tot_dup_rate;
                if(random<=tot_sil_rate)
                {
                    *mut_type='w';                              /* gene deletion*/  
                    break;
                }
                else
                {                  
                    random-=tot_sil_rate; 
                    if(random<=tot_mut_locus_length_rate)
                    {
                        *mut_type='l';      /* mut locus length */
                        break;
                    }                    
                    else
                    {                        
                        random-=tot_mut_locus_length_rate;
                        if(random<=tot_mut_kin_rate)
                        {
                            *mut_type='k';                  /* mut kinetic const*/
                            break;
                        }                        
                        else
                        {
                            random-=tot_mut_kin_rate;
                            if(random<=tot_mut_binding_seq_rate)
                            {
                                *mut_type='c';              /* mut binding seq*/
                                break;
                            }                            
                            else
                            {
                                random-=tot_mut_binding_seq_rate;
                                if(random<=tot_mut_koff_rate)
                                {                                
                                    *mut_type='f';          /* mut koff*/
                                    break;
                                }                               
                                else
                                {                                    
                                    random-=tot_mut_koff_rate;
                                    if(random<=tot_mut_identity_rate)
                                    {
                                        *mut_type='e';           /* mut identity of a TF */
                                        break;
                                    }  
                                    else
                                    {
                                        random-=tot_mut_identity_rate;
                                        if(random<=tot_effector2TF)
                                        {
                                            *mut_type='t';      /* mut effector to regulator TF */
                                            break;
                                        }
                                        else
                                        {
                                            random-=tot_effector2TF;
                                            if(random<=tot_mut_locus_specific_tf_behavior)
                                            {
                                                *mut_type='o';
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }                
                }
            }
        }  
    }
}


/*******************************************************************************
 *
 *                              Private functions 
 * 
 ******************************************************************************/
static float mut_make_new_value(float old_val, float miu, float sigma, float upper_bound, float lower_bound, RngStream RS, Mutation *mut_record, int boudary_condition)
{
    float new_val;
    new_val=old_val;
    while(new_val==old_val) 
    {
        new_val=old_val*pow(10.0,sigma*gasdev(RS)+mutational_regression_rate*(miu-log10(old_val)));
        if(boudary_condition==BOUND_INCLUDED)
        {
            if(new_val>upper_bound)
            {
                new_val=upper_bound;
                mut_record->N_hit_bound++;
            }
            if(new_val<lower_bound)
            {
                new_val=lower_bound;
                mut_record->N_hit_bound++;
            }
        }
        else
        {
            if(new_val>=upper_bound || new_val<=lower_bound)
            {
                new_val=old_val;
                mut_record->N_hit_bound++;
            }
        }      
    }
    return new_val;
}


/*
 *We use which_protein to look for the protein encoded by a given gene copy,
 *and protein_pool to look for genes encoding a given protein. These two tables
 *are updated upon mutations.
 */
static void update_protein_pool(Genotype *genotype, int which_protein, int which_gene, char mut_type)
{
    int i, j; 
    /*Protein_pool stores the numbers of gene copies that encode a given protein, and the ids of these gene copies.
     *One important thing is that the genes encoding a given protein are not stored by the order of their ids in protein_pool.
     *To delete a gene, which might be the only gene encoding a given protein, we shift protein_pool to overwrite the to-be-deleted gene
     *We need to update the ids of the remaining genes and proteins
     *For gene duplication, the new gene is always add to the end of the list of genes encoding a given protein.
     *A new protein is also add to the end of protein_pool
     *which_protein can be updated easily; changing which protein a gene encodes is always easy. For deletion, 
     *we just shift the array to overwrite the to-be-deleted gene and update the ids of the remaining genes.*/ 
    switch (mut_type)
    {
        case 'w':/*a whole gene deletion*/        
            if(genotype->protein_pool[which_protein][0][0]==1) /* if this is the only gene copy,we also need to delete a protein*/
            { 
                /* for gene<which_gene, reduce which_protein by 1 */
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i]<which_protein)?genotype->which_protein[i]:genotype->which_protein[i]-1;//the deletion also changes the ids of proteins
                /* for gene>=which_gene, shift and reduce which_protein by 1 */                
                for(i=which_gene;i<genotype->ngenes;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i+1]>which_protein)?genotype->which_protein[i+1]-1:genotype->which_protein[i+1];  
                /* UPDATE the number of activators or that of repressors*/
                if(genotype->protein_identity[which_protein]==ACTIVATOR) //protein identity is updated in the next step
                    genotype->N_act--;
                else
                    genotype->N_rep--;
                /* remove which_protein from protein_identity and Kd */              
                for(i=which_protein;i<genotype->nproteins;i++)
                {
                    genotype->protein_identity[i]=genotype->protein_identity[i+1];                    
                    genotype->Kd[i]=genotype->Kd[i+1];                    
                }  
                /* since one protein is deleted, all genes need to recalc binding sites*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)                    
                    genotype->recalc_TFBS[i]=YES; /* recalc BS */                 
            }  
            else /*if the protein has more than one genes, then no protein is deleted*/
            {
                /*shift which_protein to delete which_gene*/
                for(i=which_gene;i<genotype->ngenes;i++)
                    genotype->which_protein[i]=genotype->which_protein[i+1];                
            }   
            /*rebuild protein_pool from which_proteins*/
            rebuild_protein_pool(genotype);            
            break;
        case 'd': /*a gene duplication*/
            /* add the duplicate to protein_pool, but do not change nproteins*/    
            genotype->protein_pool[which_protein][1][genotype->protein_pool[which_protein][0][0]]=genotype->ngenes; //append newly duplicated gene to the end
            genotype->protein_pool[which_protein][0][0]++;           
            /*update which_protein*/          
            genotype->which_protein[genotype->ngenes]=which_protein; 
            break;    
        default: /*mutation in tf binding seq, in the identity of a TF, and creating a new tf and mutation in tf Kd*/
            /* remove this copy of gene from the original protein_pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            /*shift protein_pool to delete which_gene*/
            for(;i<genotype->protein_pool[which_protein][0][0]-1;i++) 
                genotype->protein_pool[which_protein][1][i]=genotype->protein_pool[which_protein][1][i+1]; 
            /*one less gene copy to encoding which_protein*/
            genotype->protein_pool[which_protein][0][0]--;                               
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins; //put the new protein to the end
            genotype->protein_pool[genotype->nproteins][0][0]=1;
            genotype->protein_pool[genotype->nproteins][1][0]=which_gene; 
            /* finally, update protein numbers*/
            genotype->nproteins++; 
    }
}

/*To reduce the amount of calculation on the probability of binding distributions, we group gene copies
 *that are created by whole gene duplication. We call such a group a cis-reg cluster because gene copies 
 *in the group should have the same cis-reg sequence. For each cis-reg cluster we only need to calculate
 *the probability of binding distributions once. However, substitutions is cis-reg sequence can create/remove
 *binding sites, therefore we need to check whether a gene copy is still in the original cis-reg cluster 
 *after mutation.We use cisreg_cluster_pool and which_cluster to track the bi-way relation between a gene and 
 *a cis-reg cluster.*/
void update_cisreg_cluster_pool(Genotype *genotype, int which_gene, char mut_type)
{
    /*In a cis-reg cluster, gene copies are ordered ascendingly by their ids. There are no empty slots in the list 
     *of gene copies. Empty slots after the list are marked by -1. We do not track the number of gene copies in a
     *cluster. In order to count the number of gene copies correctly, marking the empty slots accurately become important.
     *Therefore, we always set all the slots in a cluster to -1, when deleting or overwrting the cluster. 
     */
    int cisreg_seq_cluster_id, i; 
    
    switch(mut_type)
    {
        case 'c': //a mutation to the binding sequence of a TF happened. this can create new clusters. 
            regroup_cisreg_clusters(genotype);
            rebuild_cisreg_cluster_pool(genotype);
            break;
            
        case 's': //the cluster contains more than 1 gene, and a new cluster is created
            genotype->which_cluster[which_gene]=genotype->N_cisreg_clusters;
            rebuild_cisreg_cluster_pool(genotype);
            break;
            
        case 'w': /*gene deletion*/          
            cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
            for(i=which_gene;i<genotype->ngenes;i++)
                genotype->which_cluster[i]=genotype->which_cluster[i+1];
            /*if a cluster contains only 1 gene, the cluster is deleted*/
            if(genotype->cisreg_cluster_pool[cisreg_seq_cluster_id][0][0]==1) 
            { 
                for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
                    genotype->which_cluster[i]=(genotype->which_cluster[i]>cisreg_seq_cluster_id)?genotype->which_cluster[i]-1:genotype->which_cluster[i];  //update cluster_id              
            }  
            genotype->ngenes--;
            rebuild_cisreg_cluster_pool(genotype);
            break;
            
        case 'd': /*gene duplication*/     
            cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
            /*the duplicate is always appended to the end of gene list. Note that 
             *ngenes has not been updated when update_cisreg_cluster is called.*/
            genotype->cisreg_cluster_pool[cisreg_seq_cluster_id][1][genotype->cisreg_cluster_pool[cisreg_seq_cluster_id][0][0]]=genotype->ngenes;
            genotype->cisreg_cluster_pool[cisreg_seq_cluster_id][0][0]++;
            /*update which_cluster*/
            genotype->which_cluster[genotype->ngenes]=cisreg_seq_cluster_id;            
            break;  
    }
}

static void update_output_protein_pool(Genotype *genotype, int which_gene, char mut_type)
{    
    int i;
    switch(mut_type)
    {
        case 'w'://whole-gene deletion
            /*first reduce gene_ids that are greater than which_gene by 1*/
            for(i=0;i<genotype->n_output_genes;i++)            
                genotype->output_gene_ids[i]=(genotype->output_gene_ids[i]>which_gene)?genotype->output_gene_ids[i]-1:genotype->output_gene_ids[i];
            /*if an output gene is deleted, shift output_protein_pool to overwrite which_gene*/
            if(genotype->is_output[which_gene]==OUTPUT_PROTEIN)
            {
                i=0;
                while(genotype->output_gene_ids[i]!=which_gene) i++; //find which_gene in output_gene_ids
                for(;i<genotype->n_output_genes;i++)
                    genotype->output_gene_ids[i]=genotype->output_gene_ids[i+1];
                genotype->n_output_genes--;
            }            
            break;    
        case 'd': //gene duplication
            /*add the duplicate to output_gene_ids*/
            genotype->output_gene_ids[genotype->n_output_genes]=genotype->ngenes;
            genotype->n_output_genes++; 
            break;
        case 't': //effector2TF
            /*remove which_gene from output_gene_ids*/
            i=0;
            while(genotype->output_gene_ids[i]!=which_gene) i++;
            for(;i<genotype->n_output_genes;i++)
                genotype->output_gene_ids[i]=genotype->output_gene_ids[i+1];
            genotype->n_output_genes--;
            genotype->is_output[which_gene]=NON_OUTPUT_PROTEIN;     
    }    
}

static void rebuild_protein_pool(Genotype *genotype)
{
    int i,j,protein_id,n_copies,n_proteins;    
    /*reset protein_pool*/
    for(i=N_SIGNAL_TF;i<genotype->nproteins;i++)
    {
        genotype->protein_pool[i][0][0]=0;
        for(j=0;j<genotype->ngenes;j++)
            genotype->protein_pool[i][1][j]=NA;
    }         
    /*rebuild protein_pool using which_protein*/
    n_proteins=0;
    for(i=N_SIGNAL_TF;i<genotype->ngenes-1;i++) //note that one gene is deleted but ngenes has not been updated yet
    {
        protein_id=genotype->which_protein[i];
        n_proteins=(protein_id>n_proteins)?protein_id:n_proteins;
        n_copies=genotype->protein_pool[protein_id][0][0];
        genotype->protein_pool[protein_id][1][n_copies]=i;
        genotype->protein_pool[protein_id][0][0]++;
    }
    genotype->nproteins=n_proteins+N_SIGNAL_TF;
}

static void update_node_family_pool(Genotype *genotype, int which_gene, char mut_type)
{
    int i, node_id;    
    node_id=genotype->which_node_family[which_gene];
    switch(mut_type)
    {
        case 'w':   //gene deletion
            for(i=which_gene;i<genotype->ngenes;i++)
                genotype->which_node_family[i]=genotype->which_node_family[i+1]; 
            /*if the node family has only one member, then delete the node family*/
            if(genotype->node_family_pool[node_id][0][0]==1)
            {
                for(i=0;i<genotype->ngenes;i++)
                    genotype->which_node_family[i]=(genotype->which_node_family[i]>node_id)?genotype->which_node_family[i]-1:genotype->which_node_family[i];                     
            }
            //using which_node_family to rebuild the pool
            rebuild_node_family_pool(genotype);
            break;
        case 'd':   //gene duplication   
            /*add the duplicate to the node family node_id*/
            genotype->node_family_pool[node_id][1][genotype->node_family_pool[node_id][0][0]]=genotype->ngenes;
            genotype->node_family_pool[node_id][0][0]++;
            genotype->which_node_family[genotype->ngenes]=node_id;
            break;
        default:
            /*if the node family has multiple member, then a new node family is produced*/
            if(genotype->node_family_pool[node_id][0][0]>1)
            {                
                /*remove the gene from the old node family*/
                i=0;
                while(genotype->node_family_pool[node_id][1][i]!=which_gene) i++;
                for(;i<genotype->node_family_pool[node_id][0][0];i++)
                    genotype->node_family_pool[node_id][1][i]=genotype->node_family_pool[node_id][1][i+1]; //shift node_family_pool to overwrite which_gene
                genotype->node_family_pool[node_id][0][0]--; //the node_family has one less member
                /*the gene is now a new node*/
                genotype->which_node_family[which_gene]=genotype->N_node_families; //the id of the new node is N_node_families
                genotype->node_family_pool[genotype->N_node_families][0][0]=1;
                genotype->node_family_pool[genotype->N_node_families][1][0]=which_gene;
                genotype->N_node_families++;
           }
    }    
}

static void rebuild_node_family_pool(Genotype *genotype)
{
    int i, j, n_members, node_id, n_node_families;
    /*reset node_family_pool*/
    for(i=N_SIGNAL_TF;i<genotype->N_node_families;i++)
    {
        genotype->node_family_pool[i][0][0]=0;
        for(j=0;j<genotype->ngenes;j++)
            genotype->node_family_pool[i][1][j]=NA;
    }
    /*reubuild from which_node_family*/
    n_node_families=0;   
    for(i=N_SIGNAL_TF;i<genotype->ngenes-1;i++) //NOTE ngenes has not been updated
    {
        node_id=genotype->which_node_family[i];
        n_node_families=(n_node_families>node_id)?n_node_families:node_id; //count node families
        n_members=genotype->node_family_pool[node_id][0][0];
        genotype->node_family_pool[node_id][1][n_members]=i;
        genotype->node_family_pool[node_id][0][0]++;
    }
    genotype->N_node_families=n_node_families+N_SIGNAL_TF; 
} 

static void regroup_cisreg_clusters(Genotype *genotype)
{    
    int i,j;
    int N_genes_to_be_sorted_this_round;
    int N_genes_to_be_sorted_next_round;
    int reference_gene,gene_to_be_sorted;
    int genes_to_be_sorted[genotype->ngenes];
    int new_cluster_id;   
    int no_difference;  
    
    /* recalculate the binding sites on every promoter*/
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=YES; 
    calc_all_binding_sites(genotype,NMIN); 
    
    /*sort genes into new cisreg clusters based on TFBSs distribution*/
    for(i=0;i<genotype->ngenes-1;i++)
        genes_to_be_sorted[i]=i+1;
    new_cluster_id=1;
    N_genes_to_be_sorted_next_round=genotype->ngenes-1;     
    while(N_genes_to_be_sorted_next_round>0)  
    {
        N_genes_to_be_sorted_this_round=N_genes_to_be_sorted_next_round;
        N_genes_to_be_sorted_next_round=0;
        reference_gene=genes_to_be_sorted[0]; //compare every gene to be sorted to this gene
        genotype->which_cluster[reference_gene]=new_cluster_id;
        for(i=1;i<N_genes_to_be_sorted_this_round;i++)
        {
            no_difference=1;
            gene_to_be_sorted=genes_to_be_sorted[i];
            if(genotype->binding_sites_num[gene_to_be_sorted]==genotype->binding_sites_num[reference_gene])
            {
                for(j=0;j<genotype->binding_sites_num[reference_gene];j++)
                {
                    if(genotype->all_binding_sites[reference_gene][j].BS_pos!=genotype->all_binding_sites[gene_to_be_sorted][j].BS_pos ||
                        genotype->all_binding_sites[reference_gene][j].tf_id!=genotype->all_binding_sites[gene_to_be_sorted][j].tf_id ||
                        genotype->all_binding_sites[reference_gene][j].mis_match!=genotype->all_binding_sites[gene_to_be_sorted][j].mis_match)
                    {
                        no_difference=0;
                        break;
                    }
                }
            }
            else
                no_difference=0;    
            /*if to_be_sorted and reference have identical TFBS distribution*/     
            if(no_difference)       
                genotype->which_cluster[gene_to_be_sorted]=new_cluster_id;
            else
            {
                genes_to_be_sorted[N_genes_to_be_sorted_next_round]=gene_to_be_sorted;
                N_genes_to_be_sorted_next_round++;
            }
        }
        new_cluster_id++;
    }
}

static void rebuild_cisreg_cluster_pool(Genotype *genotype)
{
    int i,j, n_members, cluster_id, n_clusters;
    /*reset cisreg_cluster_pool*/    
    for(i=N_SIGNAL_TF;i<genotype->N_cisreg_clusters;i++)
    {        
        for(j=0;j<genotype->cisreg_cluster_pool[i][0][0];j++)
            genotype->cisreg_cluster_pool[i][1][j]=NA;
        genotype->cisreg_cluster_pool[i][0][0]=0;
    }
    /*use genotype->which_cluster to rebuild cisreg_cluster_pool*/
    n_clusters=0;
    for(i=N_SIGNAL_TF;i<genotype->ngenes;i++) 
    {
        cluster_id=genotype->which_cluster[i];
        n_clusters=(n_clusters>cluster_id)?n_clusters:cluster_id; //count clusters
        n_members=genotype->cisreg_cluster_pool[cluster_id][0][0];
        genotype->cisreg_cluster_pool[cluster_id][1][n_members]=i;
        genotype->cisreg_cluster_pool[cluster_id][0][0]++;
    }   
    genotype->N_cisreg_clusters=n_clusters+N_SIGNAL_TF;
}
