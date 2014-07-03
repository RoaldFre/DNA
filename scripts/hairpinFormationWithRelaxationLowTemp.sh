#!/bin/sh

echo "Starting hairpinFormationWithRelaxationLowTemp.sh wrapper"
echo "Running on host $HOSTNAME"

set -ex

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
fullSequence=$2
N=$3
monteCarloSweeps=$4
zippingTemp=$5
unzippingTemp=$6
allowedUnboundedFraction=$7
nucleationFraction=$8
zippedRelaxation=$9
minZippingTime=${10}
suffix=${11} # job label when running on cluster

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=100
relaxFactor=2

echo "Computing required bounded base pairs"
# Yes, ugly.
requiredBounded=`python -c "print(int(round($N * (1 - $allowedUnboundedFraction))))"`
nucleationBounded=`python -c "print(int(round($N * $nucleationFraction)))"`
allowedBounded=0
echo "requiredBounded: $requiredBounded"


echo "Getting wait time"
wait=`$scriptsdir/hairpinRelaxationTimeVariable.py $N $relaxFactor`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/dt${timestep}/zipT${zippingTemp}_unzipT${unzippingTemp}_allowUnb${allowedUnboundedFraction}_nucl${nucleationFraction}/MCsweeps${monteCarloSweeps}_relaxFactor${relaxFactor}/N${N}_minZipTime${minZippingTime}"
destinationFile="$destinationDir/formation_${suffix}"
seedFile="$destinationDir/formation_${suffix}_randSeed"

MCfile="${destinationFile}_MCrelax"
finalRelaxFile="${destinationFile}_finalRelax"
endFile="${destinationFile}_endConfig"


echo "Temporary directory is $outputBaseDir"

echo "Making sequence"

seq=`$scriptsdir/hairpinFromStem.py $fullSequence $N`

echo "Sequence is: $seq"

echo "Making destination directory: $destinationDir"
mkdir -p $destinationDir


echo
echo "Starting initial monte carlo relaxation of $monteCarloSweeps sweeps, no base pairing, writing to $MCfile"
$main -T $zippingTemp -i m -m $monteCarloSweeps -s $seq -w "$MCfile" -n

#just to be sure!
mkdir -p $destinationDir

echo "Starting final langevin relaxation, no base pairing"
$main -t $timestep -T $zippingTemp -I $interval -P $wait -d "$MCfile" -w "$finalRelaxFile" -D "/dev/null" -n

#just to be sure!
mkdir -p $destinationDir

echo "Starting main program!"
$main -t $timestep -T $zippingTemp -I $interval -d "$finalRelaxFile" -w "$endFile" -D $outputFile -H $requiredBounded:$nucleationBounded -M $allowedBounded -O $zippingTemp -Q $unzippingTemp -U $zippedRelaxation -X f -k -e -y $seedFile -Z $minZippingTime


echo "Converting to native Octave format"
$scriptsdir/hairpinFormationToNative.sh $outputFile

#just to be sure!
mkdir -p $destinationDir

echo "Moving output files!"
mv "$outputFile" "$destinationFile"
mv "${outputBaseDir}"/* "${destinationDir}"

echo "Removing seed file"
rm $seedFile

echo "Deleting temporary directory"
rm -rf $outputBaseDir

echo "All done! *High five*!"
