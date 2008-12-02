#!/bin/bash
### Set the job name
#PBS -N network500

### Request email when job begins and ends
#PBS -m bea

### Specify email address to use for notification.
#PBS -M alexlanc@arizona.edu

### Specify the PI group found with va command
#PBS -W group_list=masel

### Set the queue to submit this job., high_priority disabled
#PBS -q high_priority 

### Set the number of cpus that will be used.
#PBS -l select=1:ncpus=1

### Important!!! Include this line for your 1p job.
### Without it, the whole cpu node containing 8 cpus will be allocated.
#PBS -l place=pack:shared 

### Specify up to a maximum of 1600 hours total cpu time for 1-processor job
#PBS -l cput=240:0:0

### Specify up to a maximum of 240 hours walltime for the job
#PBS -l walltime=240:0:0

date
/usr/bin/time /home/u30/alexlanc/src/network-code/netsim-full-500 -r 8 -p 2 -d /home/u30/alexlanc/src/network-code/multiple-pops-500 -c 1.0 -n -s 20 --timesphase 1.0 --timeg2phase 0.0 --random-replication --growthscaling=10.0 --kon 0.2225 --konafter 1e-4
date
