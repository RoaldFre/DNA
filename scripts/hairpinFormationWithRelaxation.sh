#!/bin/sh

echo "Starting hairpinFormationWithRelaxation.sh wrapper"

set -e

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
fullSequence=$2
N=$3
relaxTemp=$4
zippingTemp=$5
unzippingTemp=$6
allowUnbounded=$7
zippedRelaxation=$8
time=$9 # this is a cut off: if there is no hairpin by now: give up! XXX SET LARGE ENOUGH!
suffix=${10} # job label when running on cluster

allowBounded=0

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=50

echo "Getting wait time"

wait=`$scriptsdir/hairpinRelaxationTime.py $N`

echo "Wait time: $wait"


requiredBound=`echo "$N - $allowUnbounded" | bc`

outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/relaxT${relaxTemp}_zipT${zippingTemp}_unzipT${unzippingTemp}_allowUnb${allowUnbounded}_zippedRel${zippedRelaxation}/dt${timestep}_time${time}/N${N}"
destinationFile="$destinationDir/formation_${suffix}"


echo "Temporary directory is $outputBaseDir"
echo "Making destination directory: $destinationDir"

mkdir -p $destinationDir

echo "Making sequence"

seq=`$scriptsdir/hairpinFromStem.py $fullSequence $N`

echo "Sequence is: $seq"
echo "Starting main program!"

$main -t $timestep -T $relaxTemp -I $interval -P $time -W $wait -s $seq -D $outputFile -H $requiredBound -M $allowBounded -O $zippingTemp -Q $unzippingTemp -U $zippedRelaxation -X f -k -e

echo "Moving output files!"

mv ${outputFile} ${destinationFile}
mv "${outputFile}_temperature" "${destinationFile}_temperature" #for verification / extra info
mv "${outputFile}_endToEnd" "${destinationFile}_endToEnd" #for verification / extra info

rm -rf $outputBaseDir

echo "All done! *High five*!"
