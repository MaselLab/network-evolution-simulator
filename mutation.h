/*  
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2018 Arizona Board of Regents (University of Arizona)
 */

#ifndef MUTATION_H
#define MUTATION_H

#include "RngStream.h"
#include "netsim.h"

#define SIMPLE_SUBSTITUTION 1

/*mutation rate*/
extern float DUPLICATION;   
extern float SILENCING;  

/*mutation effect*/
extern float miu_ACT_TO_INT_RATE;
extern float miu_protein_syn_rate;
extern float miu_Kd;


/*Wrapper of mutation functions*/
void mutate(Genotype *, RngStream, Mutation *);

/*Calculate mutation flux for each type of mutation, and pick one based the flux*/
void draw_mutation(Genotype *, char *, RngStream);

/*single nucleotide mutation in cisreg-sequence*/
void mut_substitution(Genotype *, Mutation *, RngStream);

/*delete a gene and its cisreg_sequence*/
void mut_whole_gene_deletion(Genotype *,Mutation *, RngStream);

/*duplicate a gene and its cisreg_sequence*/
void mut_duplication(Genotype *,Mutation *, RngStream);

/*mutate the concensus binding sequence of a TF*/
void mut_binding_sequence(Genotype *,Mutation *, RngStream);

/*mutate r_mRNA_deg, r_protein_deg, r_transl, r_ACT_to_INT, min_N_activator_to_transc*/
void mut_kinetic_constant(Genotype *, Mutation *, RngStream);

/*mutate locus_length*/
void mut_locus_length(Genotype *, Mutation *, RngStream);

/*mutate Kd of a TF*/
void mut_Kd(Genotype *, Mutation *, RngStream);

/*activator to repressor or the reverse*/
void mut_identity(Genotype *, Mutation *, RngStream);

/* Below are functions used to replay mutations*/
void reproduce_mutate(Genotype *, Mutation *);

void reproduce_substitution(Genotype *, Mutation *);

void reproduce_whole_gene_deletion(Genotype *,Mutation *);

void reproduce_gene_duplication(Genotype *,Mutation *);

void reproduce_mut_binding_sequence(Genotype *,Mutation *);

void reproduce_mut_kinetic_constant(Genotype *, Mutation *);

void reproduce_mut_identity(Genotype *, Mutation *);

void reproduce_mut_Kd(Genotype *, Mutation *);

void reproduce_mut_locus_length(Genotype *, Mutation *);
#endif /* MUTATION_H */

