#!/bin/bash

# Run me to dump performance information for a fixed number of monomers in 
# the world, in function of the number of boxes to stdout
#
# Eg:
# $ ./generateIdealNumberOfBoxes.sh > ideal1000

set -e

main=../../src/main

N=1000 #number of monomers
P=0.00675 #time period to simulate (use smaller value with lots of monomers)
I=$P

nRuns=4

for B in `seq 1 1 200`
do
	echo -n "$B $P"
	seq=`perl -e "print 'A'x$N"`
	for run in `seq 1 $nRuns`
	do
		/usr/bin/time -f "%U" $main -s $seq -I $I -P $P -b $B > /dev/null 2>tmp
		echo -n " `cat tmp`"
	done
	echo ""
done
