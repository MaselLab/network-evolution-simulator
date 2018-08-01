/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RngStream.h"
#include "netsim.h"
#include "numerical.h"
#include "mutation.h"

const int BOUND_INCLUDED=1;
const int BOUND_EXCLUDED=0;

/*
 ************ begin of mutation functions **************
 */
/*single nucleic acid substitution in cis-reg*/
void mut_substitution(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int which_nucleotide, which_gene,i,N_BS_bf_mutation,flag_difference;
    char *Genome,nucleotide; 
    AllTFBindingSites *container; // used to store binding sites before substitution
                                  // we compare whether substituion changes any binding sites,
                                  // in order to determine whether mutation creates a unique cis-reg     
                                  // sequence whose binding configuration always needs computation     
    /*points to the current cis-reg*/
    Genome= &genotype->cisreg_seq[0][0]; 
    /*whether to simulate bias in the frequency of different substitutions*/
#if !SIMPLE_SUBSTITUTION 
    float random;
    while(1)
    {
        /*this is the nucleic acid that's going to be mutated*/ 
        which_nucleotide=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
        random=RngStream_RandU01(RS);
        if (Genome[which_nucleotide]=='a'||Genome[which_nucleotide]=='t')
        {
             if(random<=0.2)break;  // 20% of mutations hit A or T
        }	            
        else
        { 
             if(random>0.2)break; //80% hit C or G
        }
    }  
    /*this is gene hit by the mutation*/
    which_gene=which_nucleotide/CISREG_LEN;
    /*calculate and store the distribution of binding sites before mutation*/
    calc_all_binding_sites_copy(genotype,which_gene);
    container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
    N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
    for(i=0;i<N_BS_bf_mutation;i++)
    {
        container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
        container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
        container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
    }
    /*generate new nucleic acid*/
    random=RngStream_RandU01(RS);
    switch (Genome[which_nucleotide])
    {
        case 'a':
            if(random<=0.2)
                Genome[which_nucleotide]='c'; // 20% of mutation to A changes it to C
            else
            {
                if((random-0.2)<=0.3)
                    Genome[which_nucleotide]='t'; // 30% changes it to T
                else
                    Genome[which_nucleotide]='g'; // 50% changes it to G
            }
            break;
        case 't':
            if(random<=0.2)
                Genome[which_nucleotide]='g';
            else
            {
                if((random-0.2)<=0.3)
                    Genome[which_nucleotide]='a';
                else
                    Genome[which_nucleotide]='c';
            }
            break;        
        case 'c':
            if(random<=0.2)
                Genome[which_nucleotide]='g';
            else
            {
                if((random-0.2)<=0.3)
                    Genome[which_nucleotide]='t';
                else
                    Genome[which_nucleotide]='a';
            }
            break;           
        case 'g':  
            if(random<=0.2)
                Genome[which_nucleotide]='c';
            else
            {
                if((random-0.2)<=0.3)
                    Genome[which_nucleotide]='a';
                else
                    Genome[which_nucleotide]='t';
            }
            break;
    }   
#else // if do not simulate bias in substitution
    /*this is the nucleic acid that's going to be mutated*/        
    which_nucleotide=RngStream_RandInt(RS,N_SIGNAL_TF*CISREG_LEN,genotype->ngenes*CISREG_LEN-1);
    /*this is gene hit by the mutation*/
    which_gene=which_nucleotide/CISREG_LEN;
    /*calculate and store the distribution of binding sites before mutation*/
    calc_all_binding_sites_copy(genotype,which_gene);
    container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
    N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
    for(i=0;i<N_BS_bf_mutation;i++)
    {
        container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
        container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
        container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
    }
    /*generate new nucleic acid*/
    nucleotide=set_base_pair(RngStream_RandU01(RS));
    while(nucleotide==Genome[which_nucleotide]) // make sure we get a different nucleic acid
        nucleotide=set_base_pair(RngStream_RandU01(RS));
    Genome[which_nucleotide]=nucleotide;
#endif    
    /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster*/
    calc_all_binding_sites_copy(genotype,which_gene);    
    if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])
    {
        update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1);  
    }
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
           update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
    free(container);
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
    int i, N_BS_bf_mutation,which_gene,flag_difference;    
    AllTFBindingSites *container;   
    /*get the mutated gene from record*/
    which_gene=mut_record->which_gene;
    /*calculate and store the distribution of binding sites before mutation*/
    calc_all_binding_sites_copy(genotype,which_gene);
    container = malloc(genotype->N_allocated_elements*sizeof(AllTFBindingSites));
    N_BS_bf_mutation=genotype->binding_sites_num[which_gene];
    for(i=0;i<N_BS_bf_mutation;i++)
    {
        container[i].BS_pos=genotype->all_binding_sites[which_gene][i].BS_pos;
        container[i].tf_id=genotype->all_binding_sites[which_gene][i].tf_id;
        container[i].mis_match=genotype->all_binding_sites[which_gene][i].mis_match;
    }
    /*apply the mutation from record*/
    Genome= &genotype->cisreg_seq[0][0];    
    Genome[mut_record->which_nucleotide]=mut_record->nuc_diff[1]; 
    /*compare binding site bf and aft substitution to decide whether to update cisreg_cluster*/
    calc_all_binding_sites_copy(genotype,which_gene);    
    if(N_BS_bf_mutation!=genotype->binding_sites_num[which_gene])
    {
        update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1);  
    }
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
           update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
    }
    free(container);    
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
    if(genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF>1)//if there are more than one copies of genes of non-output genes  
    {
        if(genotype->n_output_genes==1)//if the effector have only one copies
        {            
            while(1) //the last copy of the effector cannot be deleted. choose a different gene
            {   
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                protein_id=genotype->which_protein[which_gene];
                if(genotype->protein_identity[protein_id][1]==NON_OUTPUT_PROTEIN)
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
            protein_id=genotype->which_protein[which_gene];
            if(genotype->protein_identity[protein_id][1]!=NON_OUTPUT_PROTEIN)//obviously we can only delete a copy of the effector gene
                break;                           
        }
    }
    /*record mutation info*/
    mut_record->which_gene=which_gene;     
    /*points the first nucleotide of the cis-reg sequence to be deleted*/
    temp = &genotype->cisreg_seq[which_gene][0];   
    /* shift the cisreg array to overwrite the cis_seq to be deleted */
    for(i=0;i<CISREG_LEN*(genotype->ngenes-1-which_gene);i++) 
    {				
        *temp=*(temp+CISREG_LEN);         
        temp++;				
    }    
    /* if the to-be-deleted gene is a tf gene*/
    protein_id=genotype->which_protein[which_gene];      
    /*if this tf has only one copy of gene, then we'll delete the binding seq of this tf and remove the tf from locus_specific_TF_behavior*/ 
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
            for(j=N_SIGNAL_TF;j<genotype->nproteins;j++)
                genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i][j+1];
        }    
    }     
    /*update total loci length before removing info of which_gene*/
    genotype->total_loci_length-=genotype->locus_length[which_gene]; 
    /* remove it from PIC_assembly, mRNAdecay, proteinDecay, translation and re_calc*/
    for(i=which_gene;i<genotype->ngenes-1;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];            
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->protein_syn_rate[i]=genotype->protein_syn_rate[i+1];  
        genotype->locus_length[i]=genotype->locus_length[i+1];
        for(j=0;j<NPROTEINS;j++)
            genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i+1][j];
        genotype->recalc_TFBS[i]=1;        
    }
    /* if the to-be-deleted gene is an effector gene, change n_output_genes as well*/
    if(genotype->protein_identity[protein_id][1]!=NON_OUTPUT_PROTEIN)
        genotype->n_output_genes--;    
    /* now change protein_pool and cisreg_cluster*/   
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);  
    genotype->ngenes--;
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
    }
    /*update total loci length before removing info of which_gene*/
    genotype->total_loci_length-=genotype->locus_length[which_gene]; 
    for(i=which_gene;i<genotype->ngenes-1;i++)
    {
        genotype->min_N_activator_to_transc[i]=genotype->min_N_activator_to_transc[i+1];
        genotype->active_to_intermediate_rate[i]=genotype->active_to_intermediate_rate[i+1];            
        genotype->mRNA_decay_rate[i]=genotype->mRNA_decay_rate[i+1];
        genotype->protein_decay_rate[i]=genotype->protein_decay_rate[i+1];
        genotype->protein_syn_rate[i]=genotype->protein_syn_rate[i+1]; 
        genotype->locus_length[i]=genotype->locus_length[i+1];
        for(j=0;j<NPROTEINS;j++)
            genotype->locus_specific_TF_behavior[i][j]=genotype->locus_specific_TF_behavior[i+1][j];
        genotype->recalc_TFBS[i]=1;
    }
    if(genotype->protein_identity[protein_id][1]!=NON_OUTPUT_PROTEIN)
        genotype->n_output_genes--; 
    update_protein_pool(genotype,protein_id,which_gene,'w'); 
    update_cisreg_cluster(genotype,which_gene,'w',NULL,-1,-1);
    genotype->ngenes--;   
}

/*
 * Duplicate a whole cis-reg sequence
 */
void mut_duplication(Genotype *genotype, Mutation *mut_record, RngStream RS) 
{
    int which_gene, i, protein_id;
    char *temp1, *temp2;   
    
    if(genotype->n_output_genes>=MAX_OUTPUT_GENES) // too many effector gene copies
    {     
        //note that it's not possible to have too many effector gene copies and too many tf gene copies at the same time
        //because that will make duplication rate 0.
        while(1)
        {
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
            protein_id=genotype->which_protein[which_gene];
            if(genotype->protein_identity[protein_id][1]==-1)//if which_gene is non-output
                break;
        }        
    }
    else
    {
        if(genotype->ngenes-genotype->n_output_genes>=MAX_NON_OUTPUT_GENES) // too many non-output gene copies
        {
            while(1)
            {
                which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
                protein_id=genotype->which_protein[which_gene];
                if(genotype->protein_identity[protein_id][1]!=-1)
                    break;
            }    
        }            
        else //any gene can be duplicated
            which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1); 
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
    for(i=0;i<NPROTEINS;i++)
        genotype->locus_specific_TF_behavior[genotype->ngenes][i]=genotype->locus_specific_TF_behavior[which_gene][i];
    genotype->recalc_TFBS[genotype->ngenes]=1;
    /*update total loci length*/
    genotype->total_loci_length+=genotype->locus_length[which_gene];
    /* update protein_pool*/
    protein_id=genotype->which_protein[which_gene];    
    update_protein_pool(genotype,protein_id,which_gene,'d'); 
    /* update cisreg_cluster*/    
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    /* update gene numbers*/  
    if(genotype->protein_identity[protein_id][1]!=-1)//note duplication do not change nproteins           
        genotype->n_output_genes++;     
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
    /*update total loci length*/
    genotype->total_loci_length+=genotype->locus_length[which_gene];
    for(i=0;i<NPROTEINS;i++)
        genotype->locus_specific_TF_behavior[genotype->ngenes][i]=genotype->locus_specific_TF_behavior[which_gene][i];    
    genotype->recalc_TFBS[genotype->ngenes]=1;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'d');     
    update_cisreg_cluster(genotype,which_gene,'d',NULL,-1,-1);
    if(genotype->protein_identity[protein_id][1]!=-1)
        genotype->n_output_genes++;       
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
    which_gene=RngStream_RandInt(RS,0,genotype->ngenes-1);  
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
        /*update protein pool*/
        update_protein_pool(genotype,protein_id,which_gene,'c');          
        /* Update locus_specific_TF_behavior: assuming mutation to binding sequence changes nothing to the locus-specific behavior */   
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else /*if this tf has only one copy of gene, no new slot is required*/
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
    }
    /*mutate the binding sequence*/
    /*mutation only changes one nucleotide in the binding sequence*/
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
            tf_binding_seq_rc[TF_ELEMENT_LEN-which_nucleotide-1]='a'; break;
    }  
    /* The binding sites on every promoter needs recalculation */
    for(i=0;i<genotype->ngenes;i++)    
        genotype->recalc_TFBS[i]=YES;
    /*decide whether to update cisreg clusters. Mutation to binding seq may differentiate bs distributions among genes in a cluster*/       
    int new_clusters[NGENES][NGENES];  
    int genes_in_cluster[NGENES];
    int N_genes_in_cluster,no_difference,reference_gene,gene_to_be_sorted;
    int N_new_clusters,N_genes_in_new_cluster,j,k;
    calc_all_binding_sites(genotype); 
    i=N_SIGNAL_TF;
    while(genotype->cisreg_cluster[i][0]!=NA)/*check each cisreg cluster*/
    {       
        for(j=0;j<NGENES;j++)
        {
            for(k=0;k<NGENES;k++)
                new_clusters[j][k]=NA;
        }  
        /*copy an original cluster and count genes in the cluster*/
        N_genes_in_cluster=0;
        while(genotype->cisreg_cluster[i][N_genes_in_cluster]!=NA)
        {
            genes_in_cluster[N_genes_in_cluster]=genotype->cisreg_cluster[i][N_genes_in_cluster];
            N_genes_in_cluster++;
        }
        N_new_clusters=0; 
        /*the while loop below sort genes in a cluster into groups based on whether they have the same BS distributions*/ 
        /*We use one gene in the original cluster as a reference, and sort all genes 
         *that have the same binding sites as the reference gene into a new cluster. 
         *Genes that are different from the reference gene are sorted again similarly 
         *through iterations. After all genes are sorted, we check if the original 
         *cluster turns into multiple new cluster.*/
        while(N_genes_in_cluster>0)//run until every gene in the original cluster has been sorted into new clusters
        {
            reference_gene=genes_in_cluster[0];//compare each gene to the first gene in cluster            
            N_genes_in_new_cluster=0;
            which_gene=0;//start comparison from the first gene in the cluster
            while(which_gene<N_genes_in_cluster)
            {
                no_difference=1;
                gene_to_be_sorted=genes_in_cluster[which_gene];                
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
                if(no_difference)//if the gene has the same binding sites as the reference gene
                {
                    /*put the gene into the new cluster*/
                    new_clusters[N_new_clusters][N_genes_in_new_cluster]=genes_in_cluster[which_gene];
                    N_genes_in_new_cluster++;
                    /*shift to remove the gene from the copy of the original cluster*/
                    for(j=which_gene;j<N_genes_in_cluster-1;j++)
                        genes_in_cluster[j]=genes_in_cluster[j+1];
                    N_genes_in_cluster--;                
                }
                else //otherwise, compare next gene with the reference gene                
                    which_gene++;                
            }
            N_new_clusters++;
        }
        if(N_new_clusters!=1)//if the original cluster turns into multiple clusters
            update_cisreg_cluster(genotype,-1,'c',new_clusters,N_new_clusters,i);
        i++;
    }
    /* Calling calc_all_binding_sites reset recalc_TFBS to 0, so we need to turn them on again */
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
        update_protein_pool(genotype,protein_id,which_gene,'c');  
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }    
    else 
    {    
        tf_binding_seq=&genotype->tf_binding_seq[protein_id][0];
        tf_binding_seq_rc=&genotype->tf_binding_seq_rc[protein_id][0];
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
    for(i=0;i<genotype->ngenes;i++)
        genotype->recalc_TFBS[i]=YES;
    calc_all_binding_sites(genotype);
    int new_clusters[NGENES][NGENES],genes_in_cluster[NGENES];
    int N_genes_in_cluster,no_difference,reference_gene,gene_to_be_sorted;
    int N_new_clusters,N_genes_in_new_cluster,j,k;
    i=N_SIGNAL_TF;
    while(genotype->cisreg_cluster[i][0]!=NA)
    {        
        N_new_clusters=0;
        N_genes_in_cluster=0;
        for(j=0;j<NGENES;j++)
        {
            for(k=0;k<NGENES;k++)
                new_clusters[j][k]=NA;
        }  
        while(genotype->cisreg_cluster[i][N_genes_in_cluster]!=NA)
        {            
            genes_in_cluster[N_genes_in_cluster]=genotype->cisreg_cluster[i][N_genes_in_cluster];
            N_genes_in_cluster++;
        }
        while(N_genes_in_cluster>0)
        {
            reference_gene=genes_in_cluster[0];
            N_genes_in_new_cluster=0;
            which_gene=0;
            while(which_gene<N_genes_in_cluster)
            {
                no_difference=YES;
                gene_to_be_sorted=genes_in_cluster[which_gene];                
                if(genotype->binding_sites_num[gene_to_be_sorted]==genotype->binding_sites_num[reference_gene])
                {
                    for(j=0;j<genotype->binding_sites_num[reference_gene];j++)
                    {
                        if(genotype->all_binding_sites[reference_gene][j].BS_pos!=genotype->all_binding_sites[gene_to_be_sorted][j].BS_pos ||
                            genotype->all_binding_sites[reference_gene][j].tf_id!=genotype->all_binding_sites[gene_to_be_sorted][j].tf_id ||
                            genotype->all_binding_sites[reference_gene][j].mis_match!=genotype->all_binding_sites[gene_to_be_sorted][j].mis_match)
                        {
                            no_difference=NO;
                            break;
                        }
                    }
                }
                else
                    no_difference=NO;                
                if(no_difference)
                {
                    new_clusters[N_new_clusters][N_genes_in_new_cluster]=genes_in_cluster[which_gene];
                    N_genes_in_new_cluster++;
                    for(j=which_gene;j<N_genes_in_cluster-1;j++)
                        genes_in_cluster[j]=genes_in_cluster[j+1];
                    N_genes_in_cluster--;                
                }
                else 
                    which_gene++;
                
            }
            N_new_clusters++;
        }
        if(N_new_clusters!=1)
            update_cisreg_cluster(genotype,-1,'c',new_clusters,N_new_clusters,i);
        i++;
    }    
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
                update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
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

float mut_make_new_value(float old_val, float miu, float sigma, float upper_bound, float lower_bound, RngStream RS, Mutation *mut_record, int boudary_condition)
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
            update_cisreg_cluster(genotype,which_gene,'s',NULL,-1,-1); 
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
    
    if(genotype->locus_length[which_gene]>=MAX_GENE_LENGTH)
        genotype->locus_length[which_gene]-=1;
    else if(genotype->locus_length[which_gene]<=MAX_GENE_LENGTH)
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
    tf_id = RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1); // the first sensor tf must be an activator, therefore is not subject to mutation
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
        update_protein_pool(genotype,protein_id,tf_id,'e');
        /* update_protein_pool put the new protein at nproteins and then increases nproteins by 1, 
        * so the new protein is at nproteins-1 now. Note that N_act and N_rep is updated in update_protein_pool*/ 
        genotype->protein_identity[genotype->nproteins-1][0]=(genotype->protein_identity[genotype->nproteins-1][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR;    
        /*the locus specific behavior is not flipped by this mutation*/
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->protein_identity[genotype->nproteins-1][0];
    }
    else //otherwise we just flip the property of an exisiting TF
    {
        genotype->protein_identity[protein_id][0]=(genotype->protein_identity[protein_id][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id][0]==ACTIVATOR)//mutate to activator
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
            genotype->locus_specific_TF_behavior[i][protein_id]=genotype->protein_identity[protein_id][0];
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
        update_protein_pool(genotype,protein_id,tf_id,'e');  
        genotype->protein_identity[genotype->nproteins-1][0]=(genotype->protein_identity[genotype->nproteins-1][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->protein_identity[genotype->nproteins-1][0];
    }
    else
    {
        genotype->protein_identity[protein_id][0]=(genotype->protein_identity[protein_id][0]==ACTIVATOR)?REPRESSOR:ACTIVATOR; 
        if(genotype->protein_identity[protein_id][0]==ACTIVATOR)
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
            genotype->locus_specific_TF_behavior[i][protein_id]=genotype->protein_identity[protein_id][0];
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
    tf_id=RngStream_RandInt(RS,0,genotype->ngenes-1);
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
        update_protein_pool(genotype,protein_id,tf_id,'f');
        /* update_protein_pool put the new protein at nproteins and then increases nproteins by 1, 
         * so the new protein is at nproteins-1 now.*/
        genotype->Kd[genotype->nproteins-1]=new_Kd;  
        /* Update locus_specific_TF_behavior: assuming mutation to Kd changes nothing to the locus-specific behavior */   
        for(i=N_SIGNAL_TF;i<genotype->ngenes;i++)                            
            genotype->locus_specific_TF_behavior[i][genotype->nproteins-1]=genotype->locus_specific_TF_behavior[i][protein_id];  
    }                                                                                                           
    else
        genotype->Kd[protein_id]=new_Kd;     
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
    int which_gene,protein_id;
    while(1)
    {
        which_gene=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
        if(genotype->protein_identity[genotype->which_protein[which_gene]][1]!=-1)
            break;
    }
    mut_record->which_gene=which_gene;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'t');
    genotype->n_output_genes--;
}

void reproduce_effector2TF(Genotype *genotype, Mutation *mut_record)
{
    int which_gene, protein_id;
    which_gene=mut_record->which_gene;
    protein_id=genotype->which_protein[which_gene];
    update_protein_pool(genotype,protein_id,which_gene,'t');
    genotype->n_output_genes--;
}

/* mut gene specific tf behavior*/
void mut_locus_specific_tf_behavior(Genotype *genotype, Mutation *mut_record, RngStream RS)
{
    int gene_id, protein_id, which_cluster;
    /*which locus to mutate*/
    gene_id=RngStream_RandInt(RS,N_SIGNAL_TF,genotype->ngenes-1);
    /*which TF is flipped*/
    protein_id=RngStream_RandInt(RS,0,genotype->nproteins-1);
    /*flip the chosen TF*/
    genotype->locus_specific_TF_behavior[gene_id][protein_id]=(genotype->locus_specific_TF_behavior[gene_id][protein_id]==ACTIVATOR)?REPRESSOR:ACTIVATOR;
    /*if gene_id does not have a unique cis-reg, this mutation separate gene_id from its original cisreg cluster*/
    which_cluster=genotype->which_cluster[gene_id];    
    /*if the cis-regulatory sequence of gene_id is not unique, the mutation kicks gene_id out of its cisreg_cluster*/
    if(genotype->cisreg_cluster[which_cluster][1]!=NA)
        update_cisreg_cluster(genotype,gene_id,'s',NULL,-1,-1);    
    mut_record->which_gene=gene_id;
    mut_record->which_protein=protein_id;   
    genotype->recalc_TFBS[gene_id]=YES;    
}

void reproduce_mut_locus_specific_tf_behavior(Genotype *genotype, Mutation *mut_record)
{
    genotype->locus_specific_TF_behavior[mut_record->which_gene][mut_record->which_protein]=(genotype->locus_specific_TF_behavior[mut_record->which_gene][mut_record->which_protein]==ACTIVATOR)?REPRESSOR:ACTIVATOR;
    if(genotype->cisreg_cluster[genotype->which_cluster[mut_record->which_gene]][1]!=NA)
        update_cisreg_cluster(genotype,mut_record->which_gene,'s',NULL,-1,-1);    
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
void reproduce_mutate(Genotype *genotype, Mutation *mut_record,RngStream RS)
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
    int N_target_genes;
    
    /* duplication rate*/  
    tot_dup_rate=0.0;
    N_target_genes=genotype->ngenes-N_SIGNAL_TF;//NA to sensor TF 
    if(genotype->ngenes-genotype->n_output_genes>=MAX_NON_OUTPUT_GENES)//too many non-output genes
        N_target_genes-=genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF; //do not duplicate non-output gene anymore
    if(genotype->n_output_genes>=MAX_OUTPUT_GENES)//too many effector gene
        N_target_genes-=genotype->n_output_genes;//do not duplicate effector gene anymore
    tot_dup_rate=N_target_genes*DUPLICATION;
    tot_mut_rate+=tot_dup_rate; 
    
    /* silencing rate*/ 
    N_target_genes=0;
    if(genotype->ngenes-genotype->n_output_genes-N_SIGNAL_TF>1)//if there's more than one copy of non-output genes
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
    tot_mut_rate+=tot_mut_binding_seq_rate;    
    
    /* mut in identity*/
    tot_mut_identity_rate=(genotype->ngenes-N_SIGNAL_TF)*MUT_identity; // NA to the sensor TF gene
    tot_mut_rate+=tot_mut_identity_rate;
    
    /* mut in Kd*/
    tot_mut_koff_rate=genotype->ngenes*MUT_Kd;  
    tot_mut_rate+=tot_mut_koff_rate;    
    
    /* mut locus length*/
    tot_mut_locus_length_rate=MUT_LOCUS_LENGTH_RATE*genotype->total_loci_length;
    tot_mut_rate+=tot_mut_locus_length_rate;
    
    /*effector to regular TF*/
    if(genotype->n_output_genes>1)
        tot_effector2TF=genotype->n_output_genes*MUT_effector_to_TF;
    else
        tot_effector2TF=0.0;
    tot_mut_rate+=tot_effector2TF;
    
    /*gene-specific behavior of TF*/
    tot_mut_locus_specific_tf_behavior=(genotype->ngenes-N_SIGNAL_TF)*MUT_locus_specific_tf_behavior;
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

/*
 *We use which_protein to look for the protein encoded by a given gene copy,
 *and protein_pool to look for genes encoding a given protein. These two tables
 *are updated upon mutations.
 */
void update_protein_pool(Genotype *genotype, int which_protein, int which_gene, char mut_type)
{
    int i, j, gene_id, which_output, tf_family_id; 
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
                /*
                 * UPDATE protein_pool for protein>=which_protein
                 */
                /*shift protein>which_protein to overwrite the to-be-deleted protein*/
                for(i=which_protein;i<genotype->nproteins-1;i++)  
                {   
                    /*reset the portion of protein_pool to be overwritten*/
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                        genotype->protein_pool[i][1][j]=NA;
                    /*overwrite*/
                    genotype->protein_pool[i][0][0]=genotype->protein_pool[i+1][0][0]; //this is number of gene copies encoding a protein
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                    {
                        gene_id=genotype->protein_pool[i+1][1][j];//these are the gene copies encoding a protein
                        /*note that deletion changes the ids of the remaining genes!!! Any gene that is greater than which_gene is reduced by one*/
                        genotype->protein_pool[i][1][j]=(gene_id>which_gene)?gene_id-1:gene_id;
                    }   
                }
                /* 
                 * update TF_family_pool
                 */
                tf_family_id=genotype->which_TF_family[which_protein];
                if(genotype->TF_family_pool[tf_family_id][0][0]==1) /*there is one member in the family, so delete the family*/
                {
                    /*shift tf_family_id to overwrite*/
                    for(i=tf_family_id;i<genotype->nTF_families;i++)
                    {
                        /*reset entries of family i*/
                        for(j=0;j<genotype->TF_family_pool[i][0][0];j++)
                            genotype->TF_family_pool[i][1][j]=NA;
                        /*copy entries of family i+1 to family i. Proteins with id > protein_id needs a new id */
                        for(j=0;j<genotype->TF_family_pool[i+1][0][0];j++)
                            genotype->TF_family_pool[i][1][j]=(genotype->TF_family_pool[i+1][1][j]>which_protein)?genotype->TF_family_pool[i+1][1][j]-1:genotype->TF_family_pool[i+1][1][j];
                        genotype->TF_family_pool[i][0][0]=genotype->TF_family_pool[i+1][0][0];
                    }
                    /* update the protein id for tf family < tf_family_id*/
                    for(i=0;i<tf_family_id;i++)
                    {
                        for(j=0;j<genotype->TF_family_pool[i][0][0];j++)
                            genotype->TF_family_pool[i][1][j]=(genotype->TF_family_pool[i][1][j]>which_protein)?genotype->TF_family_pool[i][1][j]-1:genotype->TF_family_pool[i][1][j];
                    }
                    /*update nTF_family*/
                    genotype->nTF_families--;
                    /*shift which_tf_family and update the family id of protein with family id > tf_family_id*/
                    for(i=which_protein;i<genotype->nproteins-1;i++)
                        genotype->which_TF_family[i]=(genotype->which_TF_family[i+1]>tf_family_id)?genotype->which_TF_family[i+1]-1:genotype->which_TF_family[i+1];
                    /*update the family id for proteins < which_protein*/
                    for(i=0;i<which_protein;i++)
                        genotype->which_TF_family[i]=(genotype->which_TF_family[i]>tf_family_id)?genotype->which_TF_family[i]-1:genotype->which_TF_family[i];
                }
                else /* the tf family has multiple proteins, no need to delete the family*/
                {
                    /*find where is protein_id in its family*/
                    for(i=0;i<genotype->TF_family_pool[tf_family_id][0][0];i++)
                    {
                        if(genotype->TF_family_pool[tf_family_id][1][i]==which_protein)
                            break;
                    }
                    /*shift entries in its family*/
                    for(;i<genotype->TF_family_pool[tf_family_id][0][0];i++)
                        genotype->TF_family_pool[tf_family_id][1][i]=genotype->TF_family_pool[tf_family_id][1][i+1];
                    genotype->TF_family_pool[tf_family_id][0][0]--;
                    /*update protein ids in all families*/
                    for(i=0;i<genotype->nTF_families;i++)
                    {
                        for(j=0;j<genotype->TF_family_pool[i][0][0];j++)
                            genotype->TF_family_pool[i][1][j]=(genotype->TF_family_pool[i][1][j]>which_protein)?genotype->TF_family_pool[i][1][j]-1:genotype->TF_family_pool[i][1][j];
                    } 
                    /*shift which_tf_family*/
                    for(i=which_protein;i<genotype->nproteins-1;i++)
                        genotype->which_TF_family[i]=genotype->which_TF_family[i+1];
                }        
                /*
                 * UPDATE which_protein
                 */
                /* update which_protein for gene<which_gene in which_protein*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i]<which_protein)?genotype->which_protein[i]:genotype->which_protein[i]-1;//the deletion also changes the ids of proteins
                /* shift and update which_protein for gene>=which_gene in which_protein*/                
                for(i=which_gene;i<genotype->ngenes-1;i++)
                    genotype->which_protein[i]=(genotype->which_protein[i+1]>which_protein)?genotype->which_protein[i+1]-1:genotype->which_protein[i+1];  
                /*
                 * UPDATE the number of activators or that of repressors
                 */
                if(genotype->protein_identity[which_protein][0]==ACTIVATOR) 
                    genotype->N_act--;
                else
                    genotype->N_rep--;
                /* remove which_protein from protein_identity and Kd */              
                for(i=which_protein;i<genotype->nproteins-1;i++)
                {
                    genotype->protein_identity[i][0]=genotype->protein_identity[i+1][0];                    
                    genotype->Kd[i]=genotype->Kd[i+1];                    
                }                           
                /*
                 * UPDATE output_protein_id
                 */
                which_output=genotype->protein_identity[which_protein][1];                
                if(which_output!=NON_OUTPUT_PROTEIN)/*if which_protein is an effector protein*/
                {         
                    /*shift output_protein_id and update the id of proteins in output_protein_id*/
                    for(i=which_output;i<genotype->n_output_proteins-1;i++)
                        genotype->output_protein_id[i]=genotype->output_protein_id[i+1]-1; //output_protein_id is ordered ascendingly
                    /*update protein_identity similarly*/
                    for(i=which_protein;i<genotype->nproteins-1;i++)
                        genotype->protein_identity[i][1]=(genotype->protein_identity[i+1][1]==NON_OUTPUT_PROTEIN)?NON_OUTPUT_PROTEIN:genotype->protein_identity[i+1][1]-1;                    
                    genotype->n_output_proteins--;
                }
                else 
                {
                    /*just update the id of protein in output_protein_id*/
                    for(i=0;i<genotype->n_output_proteins;i++)
                        genotype->output_protein_id[i]=(genotype->output_protein_id[i]<which_protein)?genotype->output_protein_id[i]:genotype->output_protein_id[i]-1;
                    /*remove which_protein from protein_identity*/
                    for(i=which_protein;i<genotype->nproteins-1;i++)
                        genotype->protein_identity[i][1]=genotype->protein_identity[i+1][1];
                }  
                /*one less protein*/
                genotype->nproteins--;
                /* in the case, all genes need to recalc binding sites*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)                    
                    genotype->recalc_TFBS[i]=YES; /* recalc BS */ 
            }  
            else /*if the protein has more than one genes*/
            {
                /*
                 * UPDATE protein_pool for protein>=which_protein
                 */
                /* find where is this which_gene*/
                i=0;
                while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++; 
                /*shift protein_pool to overwrite which_gene, and update ids of genes*/
                j=i;
                for(;i<genotype->protein_pool[which_protein][0][0]-1;i++)
                    genotype->protein_pool[which_protein][1][i]=(genotype->protein_pool[which_protein][1][i+1]>which_gene)?genotype->protein_pool[which_protein][1][i+1]-1:genotype->protein_pool[which_protein][1][i+1];// deletion changes the ids of genes!!!
                /*also update the ids for genes before j*/
                for(i=0;i<j;i++)
                    genotype->protein_pool[which_protein][1][i]=(genotype->protein_pool[which_protein][1][i]>which_gene)?genotype->protein_pool[which_protein][1][i]-1:genotype->protein_pool[which_protein][1][i];                            
                /*one less copy encoding which_protein*/
                genotype->protein_pool[which_protein][0][0]--;                
                /*update the ids of genes in protein>which_protein*/                
                for(i=which_protein+1;i<genotype->nproteins;i++)  
                {   
                    for(j=0;j<genotype->protein_pool[i][0][0];j++)
                        genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]>which_gene)?genotype->protein_pool[i][1][j]-1:genotype->protein_pool[i][1][j];//note that deletion changes the ids of genes!!!
                }
                /*
                 * UPDATE which_protein
                 */
                /*shift which_protein to delete which_gene*/
                for(i=which_gene;i<genotype->ngenes-1;i++)
                    genotype->which_protein[i]=genotype->which_protein[i+1];                
            }
            /*
             * UPDATE protein_pool for protein<which_protein
             */           
            for(i=N_SIGNAL_TF;i<which_protein;i++)
            {
                for(j=0;j<genotype->protein_pool[i][0][0];j++)
                    genotype->protein_pool[i][1][j]=(genotype->protein_pool[i][1][j]<which_gene)?genotype->protein_pool[i][1][j]:genotype->protein_pool[i][1][j]-1;
            }            
            break;
        case 'd': /*a gene duplication*/
            /* add it to protein_pool, but do not change nproteins*/    
            genotype->protein_pool[which_protein][1][genotype->protein_pool[which_protein][0][0]]=genotype->ngenes; //append newly duplicated gene to the end
            genotype->protein_pool[which_protein][0][0]++;             
            /*update which_protein*/          
            genotype->which_protein[genotype->ngenes]=which_protein;                         
            break;
        case 'c': /*mutation in tf binding seq, creating a new tf and new tf familiy*/
            /* remove this copy of gene from the original protein_pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            /*shift to delete which_gene*/
            for(;i<genotype->protein_pool[which_protein][0][0]-1;i++) 
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            /*one less gene copy to encoding which_protein*/
            genotype->protein_pool[which_protein][0][0]--;                               
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins; //put the new protein to the end
            genotype->protein_pool[genotype->nproteins][0][0]=1;
            genotype->protein_pool[genotype->nproteins][1][0]=which_gene;
            /* make Kd for the new protein*/
            genotype->Kd[genotype->nproteins]=genotype->Kd[which_protein]; 
            /* make a new tf family and add in the new tf*/
            genotype->TF_family_pool[genotype->nTF_families][0][0]=1;
            genotype->TF_family_pool[genotype->nTF_families][1][0]=genotype->nproteins;
            genotype->which_TF_family[genotype->nproteins]=genotype->nTF_families;
            genotype->nTF_families++;
            /* update activator or repressor numbers, and protein_identity*/
            if(genotype->protein_identity[which_protein][0]==ACTIVATOR) //mutation to binding seq does not change the identity of a tf
                genotype->N_act++;
            else
                genotype->N_rep++;
            genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];  
            /* Does this gene encodes effector? */                      
            if(genotype->protein_identity[which_protein][1]!=NON_OUTPUT_PROTEIN) /*Yes*/
            {
                /*the new protein still encodes effector*/
                genotype->protein_identity[genotype->nproteins][1]=genotype->n_output_proteins;
                /*effector is now encoded by an extra protein*/
                genotype->output_protein_id[genotype->n_output_proteins]=genotype->nproteins;                   
                genotype->n_output_proteins++;
            }
            else
            {
                genotype->protein_identity[genotype->nproteins][1]=NON_OUTPUT_PROTEIN;
            }              
            /* finally, update protein numbers*/
            genotype->nproteins++;            
            break;            
        case 'e': /*mutation in the identity of a TF, creating a new tf and a new tf family*/
            /* remove this copy of gene from the original protein*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0]-1;i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; 
            }
            genotype->protein_pool[which_protein][0][0]--;                       
            /* create a new protein and link it to which_gene*/
            genotype->which_protein[which_gene]=genotype->nproteins; 
            genotype->protein_pool[genotype->nproteins][0][0]=1;
            genotype->protein_pool[genotype->nproteins][1][0]=which_gene;
            /* update Kd*/
            genotype->Kd[genotype->nproteins]=genotype->Kd[which_protein];
            /* make a new tf family and add in the new tf*/
            genotype->TF_family_pool[genotype->nTF_families][0][0]=1;
            genotype->TF_family_pool[genotype->nTF_families][1][0]=genotype->nproteins;
            genotype->which_TF_family[genotype->nproteins]=genotype->nTF_families;
            genotype->nTF_families++;
            /* update protein_identity*/
            if(genotype->protein_identity[which_protein][0]==ACTIVATOR) 
                genotype->N_rep++;  /* an activator turns into a repressor */
            else
                genotype->N_act++;
            genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];
            /* Does this gene encodes effector? */                      
            if(genotype->protein_identity[which_protein][1]!=NON_OUTPUT_PROTEIN) 
            {
                /*the new protein still encodes effector*/
                genotype->protein_identity[genotype->nproteins][1]=genotype->n_output_proteins;
                /*effector is now encoded by an extra protein*/
                genotype->output_protein_id[genotype->n_output_proteins]=genotype->nproteins;                   
                genotype->n_output_proteins++;
            }
            else
            {
                genotype->protein_identity[genotype->nproteins][1]=NON_OUTPUT_PROTEIN;
            }  
            /* finally, update protein numbers*/
            genotype->nproteins++;           
            break;            
        case 'f': /*mutation in tf Kd. We only call update_protein_pool if >1 gene 
                   *copies encode the same TF, in which case a new tf is created,
                   but does not create a new tf family*/
            /* remove this copy of gene from the original protein pool*/
            i=0;
            while(genotype->protein_pool[which_protein][1][i]!=which_gene) i++;
            for(;i<genotype->protein_pool[which_protein][0][0]-1;i++) 
            {
                genotype->protein_pool[which_protein][1][i]= genotype->protein_pool[which_protein][1][i+1]; /* rearrange data array */
            }
            genotype->protein_pool[which_protein][0][0]--; 
            /* create a new protein and link it to this gene*/
            genotype->which_protein[which_gene]=genotype->nproteins; 
            genotype->protein_pool[genotype->nproteins][0][0]=1;
            genotype->protein_pool[genotype->nproteins][1][0]=which_gene;   
            /* add the new protein to the original protein's family*/
            tf_family_id=genotype->which_TF_family[which_protein];
            genotype->TF_family_pool[tf_family_id][1][genotype->TF_family_pool[tf_family_id][0][0]]=genotype->nproteins;
            genotype->TF_family_pool[tf_family_id][0][0]++;
            genotype->which_TF_family[genotype->nproteins]=tf_family_id;
            /* update protein_identity*/
            if(genotype->protein_identity[which_protein][0]==ACTIVATOR) 
                genotype->N_act++;  
            else
                genotype->N_rep++;
            genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];  
            /* Does this gene encodes effector? */                      
            if(genotype->protein_identity[which_protein][1]!=NON_OUTPUT_PROTEIN) 
            {
                /*the new protein still encodes effector*/
                genotype->protein_identity[genotype->nproteins][1]=genotype->n_output_proteins;
                /*effector is now encoded by an extra protein*/
                genotype->output_protein_id[genotype->n_output_proteins]=genotype->nproteins;                   
                genotype->n_output_proteins++;
            }
            else
            {
                genotype->protein_identity[genotype->nproteins][1]=NON_OUTPUT_PROTEIN;
            }
            /* finally, update protein numbers*/
            genotype->nproteins++;
            /* NOTE: this mutation does not change the number of genes*/
            break;
        /*
         *
         * 
         * 
         * 
         * Mutate an effector to regular TF complicates the concept of TF family.
         * Needs to evaluate the necessity of modeling such mutation 
         * 
         * 
         * 
         * 
         * 
         */
        case 't': /*an effector gene mutate into regular TF*/
            if(genotype->protein_pool[which_protein][0][0]>1) /*if there are more than 1 gene encoding the TF*/
            {
                /*the mutation creates a new non-output protein*/
                gene_id=0;
                while(genotype->protein_pool[which_protein][1][gene_id]!=which_gene) gene_id++;
                /*remove this gene from the original protein pool*/
                for(i=gene_id;i<genotype->protein_pool[which_protein][0][0];i++)
                    genotype->protein_pool[which_protein][1][i]=genotype->protein_pool[which_protein][1][i+1];
                genotype->protein_pool[which_protein][0][0]--;
                /*make a new protein*/
                genotype->protein_pool[genotype->nproteins][0][0]=1;
                genotype->protein_pool[genotype->nproteins][1][0]=which_gene;
                genotype->which_protein[which_gene]=genotype->nproteins;
                genotype->Kd[genotype->nproteins]=genotype->Kd[which_protein];
                genotype->protein_identity[genotype->nproteins][0]=genotype->protein_identity[which_protein][0];
                genotype->protein_identity[genotype->nproteins][1]=NON_OUTPUT_PROTEIN;    
                if(genotype->protein_identity[genotype->nproteins][0]==ACTIVATOR)
                    genotype->N_act++;
                else
                    genotype->N_rep++;
                genotype->nproteins++;
            }
            else /*just remove the gene from output_protein_id*/
            {
                which_output=genotype->protein_identity[which_protein][1];
                for(i=0;i<genotype->n_output_proteins;i++)
                    genotype->protein_identity[genotype->output_protein_id[i]][1]=(genotype->protein_identity[genotype->output_protein_id[i]][1]<which_output)?
                                                                                    genotype->protein_identity[genotype->output_protein_id[i]][1]:genotype->protein_identity[genotype->output_protein_id[i]][1]-1;
                for(i=which_output;i<genotype->n_output_proteins-1;i++)
                    genotype->output_protein_id[i]=genotype->output_protein_id[i+1];
                
                genotype->protein_identity[which_protein][1]=NON_OUTPUT_PROTEIN;
                genotype->n_output_proteins--;
            }            
            break;
    }
}

/*To reduce the amount of calculation on the probability of binding distributions, we group gene copies
 *that are created by whole gene duplication. We call such a group a cis-reg cluster because gene copies 
 *in the group should have the same cis-reg sequence. For each cis-reg cluster we only need to calculate
 *the probability of binding distributions once. However, substitutions is cis-reg sequence can create/remove
 *binding sites, therefore we need to check whether a gene copy is still in the original cis-reg cluster 
 *after mutation.We use cisreg_cluster and which_cluster to track the bi-way relation between a gene and 
 *a cis-reg cluster.*/
void update_cisreg_cluster(Genotype *genotype, int which_gene, char mut_type, int new_clusters[NGENES][NGENES], int N_new_clusters, int original_cluster_id)
{
    /*In a cis-reg cluster, gene copies are ordered ascendingly by their ids. There are no empty slots in the list 
     *of gene copies. Empty slots after the list are marked by -1. We do not track the number of gene copies in a
     *cluster. In order to count the number of gene copies correctly, marking the empty slots accurately become important.
     *Therefore, we always set all the slots in a cluster to -1, when deleting or overwrting the cluster. 
     */
    int cisreg_seq_cluster_id, cisreg_seq_cluster_id_copy, i, j, last_cluster; 
    
    switch(mut_type)
    {
        case 'c': //a mutation to the binding sequence of a TF happened. this can create new clusters. 
            /*reset the original cisreg_cluster*/
            i=0;
            while(genotype->cisreg_cluster[original_cluster_id][i]!=NA)
            {
                genotype->cisreg_cluster[original_cluster_id][i]=NA;
                i++;
            }
            /*assign the first new cluster to the original cluster*/
            i=0;
            while(new_clusters[0][i]!=NA)
            {
                genotype->cisreg_cluster[original_cluster_id][i]=new_clusters[0][i];
                genotype->which_cluster[new_clusters[0][i]]=original_cluster_id;
                i++;
            }
            /*find a empty slot in cisreg_cluster for the rest of new clusters*/    
            last_cluster=N_SIGNAL_TF;
            while(genotype->cisreg_cluster[last_cluster][0]!=NA)last_cluster++;//an empty cluster has all of its slots marked by -1
            /*assign the rest of new clusters into empty slots*/       
            for(i=1;i<N_new_clusters;i++)
            {
                /*reset the empty cluster, just in case*/
                j=0;
                while(genotype->cisreg_cluster[last_cluster][j]!=NA)j++;
                /*assign new clusters*/
                j=0;
                while(new_clusters[i][j]!=NA)
                {
                    genotype->cisreg_cluster[last_cluster][j]=new_clusters[i][j];
                    genotype->which_cluster[new_clusters[i][j]]=last_cluster;
                    j++;
                }
                last_cluster++;
            }        
            break;
            
        case 's':
            /*find the cis-reg cluster of which gene*/
            cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
            if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=NA) //if a cluster contains more than 1 gene
            {   
                /*then remove which_gene from the original cluster*/        
                i=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=which_gene) i++;
                /*shift the slots to overwrite which_gene*/
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=NA)
                {
                    genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1];
                    i++;
                }                     
                /* and create a new cluster*/
                last_cluster=N_SIGNAL_TF; // start from the first non-signal-tf gene
                while(genotype->cisreg_cluster[last_cluster][0]!=NA) last_cluster++; 
                genotype->cisreg_cluster[last_cluster][0]=which_gene;
                genotype->which_cluster[which_gene]=last_cluster;
            }
            break;
            
        case 'w': /*is gene deletion*/
            /*find the cis-reg cluster of which gene*/
            cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
            if(genotype->cisreg_cluster[cisreg_seq_cluster_id][1]!=NA) //if a cluster contains more than 1 gene
            {   
                /*then remove which_gene from the original cluster*/        
                i=0;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=which_gene) i++;
                /*shift to overwrite which_gene, and update ids of the remaining genes*/
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=NA)
                {      
                    /*note that gene copies are ordered ascendingly*/
                    genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1]==NA)?NA:(genotype->cisreg_cluster[cisreg_seq_cluster_id][i+1]-1); 
                    i++;
                }               
                /*take care of clusters>cisreg_seq_cluster_id*/
                cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id+1;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=NA)
                {
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=NA)
                    {
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]<which_gene)?genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]:genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]-1;
                        i++;
                    }
                    cisreg_seq_cluster_id_copy++;
                }              
                /*update which_cluster*/
                for(i=which_gene;i<genotype->ngenes;i++)                              
                    genotype->which_cluster[i]=genotype->which_cluster[i+1];               
            }
            else//if a cluster contains only 1 gene
            {   
                /*need to shift cisreg_cluster*/ 
                cisreg_seq_cluster_id_copy=cisreg_seq_cluster_id;
                while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][0]!=NA)
                {
                    /*reset cluster=cisreg_seq_cluster_id_copy*/
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]!=NA)
                    {
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=NA;
                        i++;
                    }
                    /*then copy from cluster=cisreg_seq_cluster_id_copy+1*/                    
                    i=0;
                    while(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]!=NA)
                    {                        
                        genotype->cisreg_cluster[cisreg_seq_cluster_id_copy][i]=(genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]<which_gene)?genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]:genotype->cisreg_cluster[cisreg_seq_cluster_id_copy+1][i]-1; /*note deletion changes gene thread_ID*/
                        i++;
                    }
                    cisreg_seq_cluster_id_copy++;
                }
                /*shift which_cluster and update cluster ids for gene>=which_gene*/                
                for(i=which_gene;i<genotype->ngenes;i++)                                
                    genotype->which_cluster[i]=(genotype->which_cluster[i+1]<cisreg_seq_cluster_id)?genotype->which_cluster[i+1]:genotype->which_cluster[i+1]-1;
                /*update cluster ids for gene<which_gene*/
                for(i=N_SIGNAL_TF;i<which_gene;i++)
                    genotype->which_cluster[i]=(genotype->which_cluster[i]>cisreg_seq_cluster_id)?genotype->which_cluster[i]-1:genotype->which_cluster[i];
            }
            /*update gene ids in clusters<cisreg_seq_cluster_id*/
            for(i=N_SIGNAL_TF;i<cisreg_seq_cluster_id;i++)
            {
                j=0;
                while(genotype->cisreg_cluster[i][j]!=NA)
                {
                    genotype->cisreg_cluster[i][j]=(genotype->cisreg_cluster[i][j]<which_gene)?genotype->cisreg_cluster[i][j]:genotype->cisreg_cluster[i][j]-1;
                    j++;
                }
            }      
            break;
            
        case 'd': /*gene duplication*/     
            cisreg_seq_cluster_id=genotype->which_cluster[which_gene];
            /*find an empty slot in cisreg_seq_cluster_id*/
            i=0;
            while(genotype->cisreg_cluster[cisreg_seq_cluster_id][i]!=NA) i++;
            /*the duplicate is always appended to the end of gene list. Note that 
             *ngenes has not been updated when update_cisreg_cluster is called.*/
            genotype->cisreg_cluster[cisreg_seq_cluster_id][i]=genotype->ngenes;
            /*update which_cluster*/
            genotype->which_cluster[genotype->ngenes]=cisreg_seq_cluster_id;            
            break;  
    }
}

/****************** end of mutation functions *********************************/

