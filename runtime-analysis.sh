#!/bin/bash
outfile=$1
for kon in 0.0001 0.0005
do 
    for rand in 4 5
    do
	echo -n "$kon $rand " >> $outfile
	/usr/bin/time -f "%E %U %S" -a -o $outfile ./netsim-selection -r 4 -p 2 -d selection -c -1.0 --kon $kon
    done
done