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
allowedUnboundedFraction=$7
zippedRelaxation=$8
minZippingTime=$9
suffix=${10} # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=100

echo "Computing required bounded base pairs"
# Yes, ugly.
requiredBounded=`python -c "print(int(round($N * (1 - $allowedUnboundedFraction))))"`
allowedBounded=0
echo "requiredBounded: $requiredBounded"


echo "Getting wait time"
wait=`$scriptsdir/hairpinRelaxationTime.py $N`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/relaxT${relaxTemp}_zipT${zippingTemp}_unzipT${unzippingTemp}_allowUnb${allowedUnboundedFraction}_zippedRel${zippedRelaxation}_minZipTime${minZippingTime}/dt${timestep}/N${N}"
destinationFile="$destinationDir/formation_${suffix}"
seedFile="$destinationDir/formation_${suffix}_randSeed"

echo "Temporary directory is $outputBaseDir"

echo "Making sequence"

seq=`$scriptsdir/hairpinFromStem.py $fullSequence $N`

echo "Sequence is: $seq"

echo "Making destination directory: $destinationDir"
mkdir -p $destinationDir

echo "Starting main program!"
$main -t $timestep -T $relaxTemp -I $interval -W $wait -s $seq -D $outputFile -H $requiredBounded -M $allowedBounded -O $zippingTemp -Q $unzippingTemp -U $zippedRelaxation -X f -k -e -y $seedFile -Z $minZippingTime


echo "Moving output files!"
mv ${outputFile} ${destinationFile}
mv "${outputFile}_temperature" "${destinationFile}_temperature" #for verification / extra info
mv "${outputFile}_endToEnd" "${destinationFile}_endToEnd" #for verification / extra info

echo "Removing seed file"
rm $seedFile

echo "Deleting temporary directory"
rm -rf $outputBaseDir

echo "All done! *High five*!"
