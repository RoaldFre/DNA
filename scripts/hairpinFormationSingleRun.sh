#!/bin/sh

set -e

destinationDirRoot=$1
N=$2
temperature=$3
allowUnbounded=$4
time=$5 # this is a cut off: if there is no hairpin by now: give up! XXX SET LARGE ENOUGH!
suffix=$6 # job label when running on cluster

set -e

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=10
wait=0

fps=10
#render="-r -R 3 -f $fps"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/formation_T${temperature}_dt${timestep}_time${time}_allowUnb${allowUnbounded}/N${N}"
destinationFile="$destinationDir/formation_${suffix}"

mkdir -p $destinationDir

seq=`perl -e "print 'C'x$N,'G'x$N"`
$main -t $timestep -E $temperature -T $temperature -I $interval -P $time -W $wait -s $seq -D $outputFile -H $allowUnbounded -X f
mv ${outputFile} ${destinationFile}

rm -rf $outputBaseDir

