#!/bin/sh

set -e

destinationDirRoot=$1
N=$2
zippingTemp=$3
unzippingTemp=$4
allowUnbounded=$5
allowBounded=$6
zippedRelaxation=$7
time=$8 # this is a cut off: if there is no hairpin by now: give up! XXX SET LARGE ENOUGH!
suffix=$9 # job label when running on cluster

set -e

main=hairpin_NO_EXCLUSION # Should be in current working dir, or in $PATH
timestep=15
interval=20
wait=0

fps=10
#render="-r -R 3 -f $fps"

requiredBound=`echo "$N - $allowedUnBounded" | bc`


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/formation_zipT${zippingTemp}_unzipT${unzippingTemp}_allowUnb${allowUnbounded}_allowB${allowBounded}_zippedRel${zippedRelaxation}/dt${timestep}_time${time}/N${N}"
destinationFile="$destinationDir/formation_${suffix}"

mkdir -p $destinationDir

seq=`perl -e "print 'C'x$N,'G'x$N"`
$main -t $timestep -T $zippingTemp -I $interval -P $time -W $wait -s $seq -D $outputFile -H $requiredBound -M $allowBounded -O $zippingTemp -Q $unzippingTemp -U $zippedRelaxation -X f
mv ${outputFile} ${destinationFile}

rm -rf $outputBaseDir

