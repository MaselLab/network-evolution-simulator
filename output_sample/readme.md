This directory lists samples of the output files. We ran the simulation on Haswell V3 28 core processor. Users can run a simulation without modifying the source files, and compare the results with the provided samples. 

Files in the four subdirectories are generated when running the simulation at other modes. See the readme of source codes for explanations of different simulation modes. 


# Details of each file
` 1. evo_summary_81.txt`
This file gives an overview of the evolution. The number “81” in the filename is the random number seed used in a simulation. Each row lists the vital information at each evolutionary step. The columns are:

**Step**: current evolutionary step

**N_tot_mut_tried**: number of mutations tried since the beginning

**N_mut_tried_this_step**: number of mutations tried at the current step 

**N_hit_bound**: accumulative number of mutations that push a kinetic constant to the biological bounds

**accepted_mut**: the mutation that is accepted at the current step

**selection_coeff**: selection coefficient of the accepted mutant

**avg_fitness**: fitness. 

**fitness1**: should be the same as **avg_fitness**. avg_fitness used to be the average of fitness1 and fitness2, which are the fitness of the network under two environments. We only made a selection condition for pulse generator under environment 1 and left environment 2 unused. 

**fitness2**: should be all zero.  

**se_avg_fitness**: standard error of the fitness. 

**se_fitness1**: should be the same as ** se_avg_fitness **.

**se_fitness2**: should be all zero.

**N_genes**: number of genes and variants (the signal is counted as a gene)

**N_effector_genes**: number of effector genes and variants

**N_proteins**: number of proteins (include the signal. Gene variants are treated as the same protein)
	 

`2. networks.txt`
This file summarizes the topology of the TRN evolved at an evolutionary step. We use a table to record the number of TFBs of each TF on each gene. For example:
```
step 50
Gene   a0     r1     r2     A3     which_protein AND_gate_capable
1      3      2      1      2           A3             N
2      1      0      1      0           r1             N
3      3      2      1      2           A3             N
4      3      2      1      1           A3             N
5      1      0      1      0           r1             N
6      1      1      2      1           r2             N
7      1      1      2      1           r2             N
8      1      1      2      1           r2             N
 ```
The table shows the TRN at step 50. Each row is a gene. The header labels TFs by **a**ctivator or **r**epressor (if the TF is also an effector then the letter is in upper case), followed by a number to specify their identity in the proteome. The signal is always an activator and is always the first protein in the proteome, therefore is labeled a0. The column **which_protein** specifies which protein the gene encodes.  In this example, the first row means that gene 1 contains 3 TFBS of the signal, 2 TFBS of repressor 1, 1 TFBSs of repressor 2, and 2 TFBSs of activator 3 which is also the effector. Row 1 also states that gene 1 encodes activator 3, and is capable of being regulated via an AND gate logic between activators (meaning that one activator is able to enhance gene expression).

The program by default does not summarize TRN at every evolutionary step. Use OUTPUT_INTERVAL (see readme of the code) to control the frequency of output.

`3. N_motifs.txt`
Each row of the file lists the numbers of motifs in the TRN evolved at an evolutionary step. We count 17 types of motifs (see **motifs.png** in this folder), which corresponds to the first 17 columns of the files. 

`4. accepted_mutation_81.txt`
Each row is the mutation that is accepted at an evolutionary step. The columns are: 

Col1: mutation type, i.e. nucleotide substitution, gene deletion, gene duplication, consensus binding sequence, gene-specific kinetic constant, identity of TF, affinity of TF (Kd), and length of gene.

Col2: the id of gene that is mutated

Col3: the nucleotide that is substituted

Col4: the new nucleotide

Col5: the code of kinetic constant that is mutated, 0 for *r_Act_to_Int*, 1 for *r_mRNA_deg*, 2 for *r_protein_syn*, and 3 for *r_protien_deg*. 

Col6: new value for mutated quantitative variable. Stored in hexadecimal form.

` 5. sim_setup_81.txt`
This file records the initial condition and selection condition of a simulation.

` 6. init_mutatable_parameters.txt`
Values of gene-specific variables at initialization. Each row is a gene (the first row is always the signal), and columns are values of *r_Act_to_Int*, *r_mRNA_deg*, *r_protein_syn*, *r_protien_deg*, *l*, and log10(Kd(0)).

` 7. end_mutatable_parameters.txt`
Values of gene-specific variables at the end of a simulation.

`8. RngSeeds.txt`
The state of random number generator at the end of an evolutionary step. The state is stored every 20 evolutionary steps by default. Modify OUTPUT_INTERVAL (see readme of the code) to change the frequency of saving.

`9. precise_fitness.txt`
Each row is the fitness of the resident genotype, expressed in hexadecimal form.

` 10. saving_point.txt`
In case the program has to be terminated before completion (e.g. running at windfall mode on a HPC), a “saving point” is made periodically so that the program can be continued. The first number in the file indicates the last evolutionary step before the program is terminated prematurely; the second number is **N_tot_mut_tried** at the last evolutionary step.  When a simulation is continued, the program replays mutation stored in `accepted_mutation_81.txt` up to the evolutionary step specified in `saving_point.txt`, load fitness of the resident genotype from `precise_fitness.txt`, and set the random number generator to the state accordingly. The file is generated every 20 evolutionary steps by default. Modify OUTPUT_INTERVAL (see readme of the code) to change the frequency of saving.

` 11. all_mutations.txt`
Similar to accepted_mutation_81.txt, but additionally list mutations that were not accepted. 

Col1: the evolutionary step at which the mutation is generated

Col2: total number of mutations generated so far

Col3 – Col8 are the same as Col - Col6 of accepted_mutation_81.txt

` 12. fitness_all_mutants.txt`
Lists the low-resolution genotype fitness of every mutant, regardless of whether the mutant is accepted. 

**Under CLEAN_UP_NETWORK mode**
` 13. after_perturbation.txt`
Fitness of a network after removing 2-mismatch TFBSs. The columns are :
Col1: evolutionary step

Col2: code that indicates which 2-mistmatch TFBSs are considered non-adaptive after we compare the original fitness of the network with the fitness when 2-mismatch TFBSs are removed. We represent TFBSs of TF X in the cis-regulatory sequence of gene Y as X -> Y. For example, signal -> effector means all the 2-mismatch TFBSs of the signal that are found in the regulatory sequences of all effector genes. We consider five basic situations and list their code in the table below. We add the codes to represent that more than one type of 2-mismatch TFBSs are non-adaptive. For example, 300 means two types of 2-mismatch TFBSs: effector -> non-effector and signal -> non-effector.

```
0 all TFBSs are adaptive
50 non-effector -> effector
100 effector -> non-effector
200 signal -> non-effector
400 effector -> effector
```
Col3: original fitness of the network

Col4: fitness when signal -> non-effector are removed

Col5: fitness when effector -> non-effector are removed

Col6: fitness when effector -> effector are removed

Col7: fitness when non-effector -> effector are removed

Col8 – Col12: standard error of fitness. Same order as Col3 – Col7 

` 14. networks_clean.txt`
Topology of TRNs after non-adaptive 2-mismatch TFBSs are removed. The format is identical to that of ‘networks.txt’.

`15. N_motifs_clean.txt`
Number of network motifs after non-adaptive 2-mismatch TFBSs are removed. The format is identical to that of ‘N_motifs.txt’.

** Under CLASSIFY_MUTATION**
`16. after_perturbation_all_mutations.txt`
Lists fitness before and after 2-mismatch TFBSs are removed. Each row is a mutation, starting from the first mutation that is trialed. The columns are:
Col1: original fitness of the network
Col2: fitness when signal -> non-effector are removed
Col3: fitness when effector -> non-effector are removed
Col4: fitness when effector -> effector are removed
Col5: fitness when non-effector -> effector are removed
Col6 – Col10: standard error of fitness. Same order as Col3 – Col7
If a mutation is accepted, then the corresponding row is all -1 because the relevant fitness is already recorded in ‘after_perturbation.txt’. 

`17. mutation_classification.txt`
This file record whether a mutation creates or destroys a given network motif. Each row is a mutation, starting from the first mutation that is trialed. The columns are:

Col1: 1 means the mutation is accepted

Col2: which 2-mismatch TFBSs in the mutant are non-adaptive. The code is identical to that in ‘after_perturbation.txt’

Col3: 1 means the mutation creates one or more I1FFLs (the ancestral genotype has no I1FFL). -1 means the mutation removes all I1FFLs. 

Col4: 1 means the mutation creates one or more NFBLs (the ancestral genotype has no NFBL). -1 means the mutation removes all NFBLs.

Col5: 1 means the mutation creates one or more overlapping I1FFLs (the ancestral genotype has no overlapping I1FFL). -1 means the mutation removes all overlapping I1FFLs.

Col6: 1 means the mutation creates one or more I1FFL+NFBL conjugates (the ancestral genotype has no I1FFL+NFBL conjugate). -1 means the mutation removes all I1FFL+NFBL conjugates.

**Under SAMPLE_EFFECTOR_EXPRESSION_LVL**
`18. effector_expression_lvl.txt`
Each line is a time course of the effector’s expression levels that is observed at an evolutionary step. The first number of each row marks the evolutionary step, followed by expression levels at each time point.  The time course starts from the beginning of stage 1 and lasts till the end of growth simulation, which is 360 minutes in total by default. The default sampling interval is 1 minute, therefore a time course contains 360 time points plus one more point that is sampled at time 0. At each time point, we average the expression levels of the effector from multiple replicates (200, by default). Note that the expression level is the total expression levels of all effector genes. Note that github limits file size, so the sample file only list time course in the initial genotype and the first 9999 evolved gentoypes.
 
**Under SAMPLE_GENE_EXRESSION**
`19. gene_i.txt`
Each row of the file is the time course of the expression levels of gene i that is observed in a replicate. See effector_expression_lvl.txt for the number of time points. 

`20. protein_i.txt`
Similar to `gene_i.txt`, but show the total expression levels of all the genes that encode a given protein i.  


` 21. gene_and_protein.txt `
This file lists which gene encodes which protein, and whether a protein functions as an effector.
