#!/bin/sh

set -e

destinationDirRoot=$1
S=$2
L=$3
zippingTemp=$4
unzippingTemp=$5
allowUnbounded=$6
allowBounded=$7
zippedRelaxation=$8
time=$9 # this is a cut off: if there is no hairpin by now: give up! XXX SET LARGE ENOUGH!
suffix=${10} # job label when running on cluster

set -e

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=20
wait=0

fps=10
#render="-r -R 3 -f $fps"

requiredBound=`echo $S - $allowUnbounded | bc`

outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/formation_CG${S}_A${L}_zipT${zippingTemp}_unzipT${unzippingTemp}_allowUnb${allowUnbounded}_allowB${allowBounded}_zippedRel${zippedRelaxation}/dt${timestep}_time${time}/N${N}"
destinationFile="$destinationDir/formation_${suffix}"

mkdir -p $destinationDir

seq=`perl -e "print 'C'x$S,'A'x$L,'G'x$S"`
$main -t $timestep -E $zippingTemp -T $zippingTemp -I $interval -P $time -W $wait -s $seq -D $outputFile -H $requiredBound -M $allowBounded -O $zippingTemp -Q $unzippingTemp -U $zippedRelaxation -X f
mv ${outputFile} ${destinationFile}

rm -rf $outputBaseDir

