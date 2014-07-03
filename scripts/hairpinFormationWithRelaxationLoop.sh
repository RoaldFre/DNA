#!/bin/sh

echo "Starting hairpinFormationWithRelaxationLoop.sh wrapper"

set -e

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
L=$2
S=$3
relaxTemp=$4
zippingTemp=$5
unzippingTemp=$6
allowedUnbounded=$7
zippedRelaxation=$8
suffix=$9 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=100

echo "Computing required bounded base pairs"
# Yes, ugly.
requiredBounded=`python -c "print($S - $allowedUnbounded)"`
allowedBounded=0
echo "requiredBounded: $requiredBounded"


echo "Getting wait time"
halfLength=`python -c "print(int(round($S + $L/2)))"`
wait=`$scriptsdir/hairpinRelaxationTime.py $halfLength`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/relaxT${relaxTemp}_zipT${zippingTemp}_unzipT${unzippingTemp}_zippedRel${zippedRelaxation}/dt${timestep}/S${S}/L${L}"
destinationFile="$destinationDir/formation_${suffix}"

echo "Temporary directory is $outputBaseDir"

echo "Making sequence"

seq=`perl -e "print 'C'x$S,'A'x$L,'G'x$S"`

echo "Sequence is: $seq"
echo "Starting main program!"

$main -t $timestep -T $relaxTemp -I $interval -W $wait -s $seq -D $outputFile -H $requiredBounded -M $allowedBounded -O $zippingTemp -Q $unzippingTemp -U $zippedRelaxation -X f -k -e


echo "Making destination directory: $destinationDir"
mkdir -p $destinationDir

echo "Moving output files!"
mv ${outputFile} ${destinationFile}
mv "${outputFile}_temperature" "${destinationFile}_temperature" #for verification / extra info
mv "${outputFile}_endToEnd" "${destinationFile}_endToEnd" #for verification / extra info

echo "Deleting temporary directory"
rm -rf $outputBaseDir

echo "All done! *High five*!"
