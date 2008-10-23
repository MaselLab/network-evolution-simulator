#!/bin/bash
### Set the job name
#PBS -N networkruntime

### Request email when job begins and ends
#PBS -m bea

### Specify email address to use for notification.
#PBS -M alexlanc@email.arizona.edu

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
outfile=runtime-profile.txt
# print out header
echo "kon rand replicate real user sys" > $outfile
for kon in 0.00001 0.0001 0.001 0.01 0.1
do 

    for ((rand=4;rand<104;rand+=1))
    do
	for replicate in 1
	do
	    echo -n "$kon $rand $replicate " >> $outfile
	    /usr/bin/time -f "%e %U %S" -a -o $outfile /home/u30/alexlanc/src/network-code/netsim-selection -r $rand -p 2 -d selection -c -1.0 --kon $kon
	done
    done
done
date
