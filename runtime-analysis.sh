#!/bin/bash
outfile=$1
# print out header
echo "kon rand replicate real user sys" > $outfile
for kon in 0.0001 0.0005 0.001
do 

    for ((rand=4;rand<104;rand+=1))
    do
	for replicate in 1
	do
	    echo -n "$kon $rand $replicate " >> $outfile
	    /usr/bin/time -f "%E %U %S" -a -o $outfile ./netsim-selection -r $rand -p 2 -d selection -c -1.0 --kon $kon
	done
    done
done