#!/bin/sh

echo "Starting hairpinFormationWithRelaxationXY.sh wrapper"

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

NXY=4
NtotStem=`python -c "print($N + $NXY)"` # use this for relaxation time etc

echo "Computing required bounded base pairs, N=$N, NXY=$NXY, NtotStem=$NtotStem"
# Yes, ugly.
requiredBounded=`python -c "print(int(round(($NtotStem) * (1 - $allowedUnboundedFraction))))"`
allowedBounded=$NXY
echo "requiredBounded: $requiredBounded"


echo "Getting wait time"
wait=`$scriptsdir/hairpinRelaxationTime.py $NtotStem`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/relaxT${relaxTemp}_zipT${zippingTemp}_unzipT${unzippingTemp}_allowUnb${allowedUnboundedFraction}_zippedRel${zippedRelaxation}_minZipTime${minZippingTime}_NXY${NXY}/dt${timestep}/N${N}"
destinationFile="$destinationDir/formation_${suffix}"
seedFile="$destinationDir/formation_${suffix}_randSeed"

echo "Temporary directory is $outputBaseDir"

echo "Making sequence"

seq=`$scriptsdir/hairpinFromStemXY.py $fullSequence $N $NXY`

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
