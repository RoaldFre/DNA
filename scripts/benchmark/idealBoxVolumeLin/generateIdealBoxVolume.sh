#!/bin/bash

set -e

main=hairpin

N=$1 #number of monomers
P=$2 #time period to simulate (use smaller value with lots of monomers)
sweeps=$3 # MC sweeps as initial relaxation
S=$4 #world size
Bmin=$5 # Minimum number of boxes
Bstep=$6
Bmax=$7 # Maximum number of boxes
nRuns=$8 # Number of runs

I=$P

for B in `seq $Bmin $Bstep $Bmax`
do
	echo -n "$B $P $S"
	seq=`perl -e "print 'A'x$N" 2>/dev/null`
	configfile=config$N
	tmpfile=tmp$N
	for run in `seq 1 $nRuns`
	do
		$main -s $seq -S $S -i m -m $sweeps -w $configfile -b $Bmax > /dev/null
		mv $configfile ${configfile}_2
		$main -I 0.2 -P 0.2 -d ${configfile}_2 -w $configfile  -b $Bmax -t 10 -g '5e11' > /dev/null
		/usr/bin/time -f "%U" $main -d $configfile -I $I -P $P -b $B -t 10 -g '5e11' > /dev/null 2>$tmpfile
		echo -n " `cat $tmpfile`"
	done
	echo ""
	rm $tmpfile $configfile ${configfile}_2
done
