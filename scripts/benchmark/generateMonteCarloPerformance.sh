#!/bin/bash

set -e

main=hairpin

N=$1
sweeps=$2
nRuns=$3

tmpfile=monteCarloPerformance_tmp

echo -n "$N $sweeps"
seq=`perl -e "print 'A'x$N"`
for run in `seq 1 $nRuns`
do
	/usr/bin/time -f "%U" $main -s $seq -i m -m $sweeps > /dev/null 2> $tmpfile
	echo -n " `cat  $tmpfile`"
done
echo ""

rm $tmpfile

