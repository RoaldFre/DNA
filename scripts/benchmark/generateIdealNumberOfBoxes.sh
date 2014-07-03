#!/bin/bash

set -e

main=hairpin

N=$1 #number of monomers
P=$2 #time period to simulate (use smaller value with lots of monomers)
Bmax=$3 # Maximum number of boxes
nRuns=$4 # Number of runs

I=$P

for B in `seq 1 1 $Bmax`
do
	echo -n "$B $P"
	seq=`perl -e "print 'A'x$N"`
	for run in `seq 1 $nRuns`
	do
		/usr/bin/time -f "%U" $main -s $seq -I $I -P $P -b $B > /dev/null 2>tmp$N
		echo -n " `cat tmp$N`"
	done
	echo ""
done
