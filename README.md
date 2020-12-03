We provide the computational model as source code, which is written in C. The source files must be compiled to produce the simulation program. We mainly used Intel C compiler (icc, version 16.0.4), but the GNU C compiler (gcc) will also work (although the outcome of a simulation will change due to different optimization to numerical calculations). 
We note the simulating the evolution of I1FFLs requires settings that are different from evolving C1FFLs, therefore the source code is also different. Below we explain how to setup the simulations described in our manuscript.

# Installation (Run the default mode)
By default, the program evolves TRNs under selection for a pulse generator. The program runs on 10 CPU cores (Haswell V3 28 core processor), and takes 1-2 days. 

We suggest using a Linux system to facilitate the installation. To run the program in the default mode, follow these steps:

1. Copy all source files (files with suffix .c and .h, and the *makefile*) to one directory. 

2. Under the same directory, create a folder and name it **result**. The folder will be used to hold output files.

3. Change directory into the directory that contains the source file. Compile source files using the command
```
    make simulator CC=icc
```
This command will create several files with suffix .o and an executable program named simulator. “CC=icc” compiles the source files with icc. By default, the compiling is done with -O3 for simulation speed, and -fp-model precise -fp-model source to ensure arithmetic operations are accurate and reproducible.

To compile with gcc, change “CC=icc” to “CC=gcc” and comment out -fp-model precise -fp-model source in *makefile*. We have noticed that when compiled with gcc, the simulation produces result different from when compiled with icc, even for the same random number seed. Enabling safe arithmetic options in gcc may solve the problem, but we haven’t tested it.

4. Execute simulator to start. On Linux, this is done with the following command
```
./simulator
```

# Output 
The simulation will generate several files when it begins to run. The size of some files, e.g. evolutionar_summary_81.txt, will keep increasing. Samples of output files and a description to their content can be found in folder **output_sample**.

# Change the parameters of selection for pulse generator
The relevant parameters are listed in main.c. Line 114 – 117 are the duration and strength of the input signal: 
```c
selection.env1.signal_strength_stage1=100.0;    
selection.env1.signal_strength_stage2=1000.0;
selection.env1.t_stage1=120.0;
selection.env1.t_stage2=2000.0; 
```
Line 125 – 130 are the parameters that are used to calculate the fitness of a pulse:
```c
selection.env1.opt_peak_response=10000.0; //optimal peak expression level of the effector
selection.env1.effector_level_before_stage2=0.1;  // target effector expression level (as a fraction of the peak level) at the end of stage 1. 
selection.env1.fitness_decay_constant=0.693; // this is the sigma squared in Eq. 1
selection.env1.min_reduction_relative_to_peak=0.2;  // minimal reduction in the effector expression level (as a fraction of the peak level) at the end of simulation  
selection.env1.window_size=5; // use a 5-minute window to average the expression level of the effector
selection.env1.max_t_mid=60.0; 
```
Line 140 - 143 specify whether to prevent I1FFL or NFBL from evolving:
```c
selection.effector_is_TF=1; // whether effector is a TF. 1 means NFBLs can evolve
resident.flag_effector_is_TF=1; // make sure also change the resident genotype    
selection.signal_ctrl_repressor=1; // whether the signal can regulate repressor. 1 means I1FFLs can evolve. If set to 0, make sure mutation in TF identity is also disabled 
resident.flag_signal_ctrl_repressor=1;
```
To prevent the evolution of I1FFLs, you also need to disable the mutational conversion between activator and repressor. This is done by setting the mutation rate at line 34 of mutation.c to zero:
```c
static const float MUT_identity= 0.0;
```
Finally, line 51 and 53 specify whether the effector is initialized as an activator (default) or repressor.  
```c
int init_N_output_act=1; //number of effector genes that are activator at initialization
int init_N_output_rep=0;
```

# mode CLEAN_UP_NETWORK

This mode cleans up the network by removing non-adaptive weak TFBSs. The program replays mutation and searches for non-adaptive weak TFBSs in four types of TFBSs that can form I1FFLs and NFBLs at each evolutionary steps (see manuscript for details). The program will remove all instances of a type of TFBSs, and if this does not significantly reduce the fitness of TRN, then the type of TFBSs are considered non-adaptive. After removing all four types of non-adaptive weak TFBSs, the program will output a new TRN and the number of motifs. 
To enable the perturbation mode, following these steps:

1. Modify line 25 and 46 of netsim.h to 
```c
#define PERTURB 1
```
```c
#define CLEAN_UP_NETWORK 1 //identify and exclude non-adaptive 2-mismatch TFBSs
```

2. Specify how much reduction in fitness is considered significant by modify line 51 of netsim.h:
```c
#define CUT_OFF_NONADAPTIVE_TFBS -0.01 
``` 

3. Make sure the random number seed and the parameters of selection for pulse generator are the same as those used to evolve the TRNs.

4. Copy *accepted_mutation_x.txt* file and *evo_summary_x.txt* to result. 

5. Compile the source files and run simulator. 

Because the program needs to measure the fitness of many TRNs, it is recommended to run the program with multiple CPUs (usually take 10 cpus two days). See **readme** in folder **output_sample** for the output files.

# mode CLASSIFY_MUTATION 

In this mode, the program determines check all mutations that were proposed during evolution simulation, and classify them into I1FFL-creating, I1FFL-destroying, NFBL-creating, and NFBL-destroying.
Before running this mode, users must first run “CLEAN_UP_NETWORK” (see above). To classify mutations, following these steps:

1. Set line 25 and 47 of netsim.h to 1:
```c
#define PHENOTYPE 1
```
```c
#define CUT_OFF_NONADAPTIVE_TFBS -0.01
```

2. Make sure the random number seed and the parameters of selection for pulse generator are the same as those used to evolve the TRNs.

3. Copy *accepted_mutation_x.txt*, *evo_summary_x.txt*, *all_mutations.txt*, *fitness_all_mutants.txt*, and *after_perturbation.txt* (generated by running “Remove non-adaptive weak TFBSs”) to result. 

4. Compile the source files and run simulator. 

Because the program needs to measure the fitness of many TRNs, it is recommended to run the program with multiple CPUs (usually take 10 cpus two to three days). See **readme** in folder **output_sample** for the output files.

# mode SAMPLE_GENE_EXPRESSION
This mode outputs a time course of the expression levels of all the genes at a particular evolutionary step. It uses the accepted_mutation_x.txt (here x is the random number seed of the simulation. We provide accepted_mutation_81.txt in folder **output_sample** as an example) of a previous simulation to replay evolution, and reproduce the genotype at a given evolutionary step. To enable this mode, following these steps:

1. Modify line 24 and 40 of netsim.h to
```c
#define PHENOTYPE 1
```
```c
#define SAMPLE_GENE_EXPRESSION 1
```

2. Copy accepted_mutaton_81.txt file (see folder **output_sample**) to result. 

3. Set the target evolutionary step at argument 6 of line 203 of main.c
```c
show_phenotype(&resident, &mut_record, &selection, init_mRNA, init_protein, 50000, RS_parallel); 
```

4. Make sure the random number seed and the parameters of selection for pulse generator are the same as those used to evolve the TRNs.

5. Compile the source code and run the simulator

This mode can be run on one or multiple CPUs and finishes in minutes. See **readme_output.pdf** in folder **output_sample** for the output files.

# mode SAMPLE_EFFECTOR_EXPRESSION
This mode is similar to SAMPLE_GENE_EXPRESSION, but it only outputs the expression time course of the effector during the first N evolutionary steps.  To enable this mode, following these steps:

1. Modify line 24 and 41 of netsim.h to
```c
#define PHENOTYPE 1
```
```c
#define SAMPLE_EFFECTOR_EXPRESSION 1
```

2. Copy accepted_mutaton_81.txt file (see folder **output_sample**) to result. 

3. Set the target evolutionary step at argument 6 of line 203 of main.c
```c
show_phenotype(&resident, &mut_record, &selection, init_mRNA, init_protein, 50000, RS_parallel); 
```

4. Make sure the random number seed and the parameters of selection for pulse generator are the same as those used to evolve the TRNs.

5. Compile the source code and run the simulator

Because the program needs to measure the fitness of many TRNs, it is recommended to run the program with multiple CPUs (usually take 10 cpus several hours). See **readme** in folder **output_sample** for the output files.

# Additional settings
## 1. Change random number seed
Random number seed is set at line 37 of main.c. It mainly controls the initial genotypes.

## 2. Change the number of parallel threads
By default, the program runs on 10 threads. To change, modify line 33 of netsim.h. Note that N_REPLICATES (line 34 of netsim.h) must be divisible by N_THREADS! 

## 3. Change output interval
By default, the program pools results of 50 evolutionary steps before writing to disk. This can be changed by modifying OUTPUT_INTERVAL at line 35 of netsim.h.
