#!/bin/sh

echo "Starting hairpinStateWithRelaxationXY.sh wrapper"

set -e

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
fullSequence=$2
N=$3
relaxTemp=$4
samplingTemp=$5
samplingTime=$6
suffix=$7 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=100

NXY=4
NtotStem=`python -c "print($N + $NXY)"` # use this for relaxation time etc

echo "Getting wait time"
wait=`$scriptsdir/hairpinRelaxationTime.py $NtotStem`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/relaxT${relaxTemp}_sampleT${samplingTemp}_time${samplingTime}_NXY${NXY}/dt${timestep}/N${N}"
destinationFile="$destinationDir/state_${suffix}"

echo "Temporary directory is $outputBaseDir"

echo "Making sequence"

seq=`$scriptsdir/hairpinFromStemXY.py $fullSequence $N $NXY`

echo "Sequence is: $seq"
echo "Starting main program!"

$main -t $timestep -T $relaxTemp -I $interval -W $wait -P $samplingTime -s $seq -D $outputFile -a $samplingTemp -X s

echo "Converting output to native Octave format"
octave -p $scriptsdir --eval "toNative('${outputFile}')"

echo "Making destination directory: $destinationDir"
mkdir -p $destinationDir

echo "Moving output files!"
mv ${outputFile} ${destinationFile}

echo "Deleting temporary directory"
rm -rf $outputBaseDir

echo "All done! *High five*!"
