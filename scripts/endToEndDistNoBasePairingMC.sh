#!/bin/bash

set -ex

destinationDirRoot=$1
fullSequence=$2
N=$3
temperature=$4
sweeps=$5
suffix=$6 # job label when running on cluster

scriptsdir="$HOME/DNA/scripts"
main=hairpin # Should be in current working dir, or in $PATH
interval=0.1 # 1 == once per MC sweep


echo "Running on $HOSTNAME"

echo "Making temporary directory"
outputBaseDir=`mktemp -d`
echo "Made temporary directory"
outputFile="$outputBaseDir/outputFrom_${suffix}"
#destinationDir="$destinationDirRoot/$fullSequence/T${temperature}/MCsweeps${sweeps}/N${N}"
destinationDir="$destinationDirRoot/MCendToEndRelax/T${temperature}/MCsweeps${sweeps}/N${N}"
destinationFile="$destinationDir/endToEnd_${suffix}"

echo "Making sequence"
#seq=`$scriptsdir/hairpinFromStem.py $fullSequence $N`
seq=`perl -e "print 'A'x$N"`
echo "Made sequence"

echo "Starting main program"
echo ""
time $main -T $temperature -I $interval -s $seq -D $outputFile -e -n -i m -m $sweeps
echo ""
echo "End of main program"

echo "Making destination directory"
mkdir -p $destinationDir
echo "Made destination directory"

echo "Converting to native format"
$scriptsdir/toNative.sh "${outputFile}_endToEnd"

echo "Moving output files"
mv "${outputFile}_endToEnd" "${destinationFile}"

echo "Removing temporary directory"
rm -rf $outputBaseDir

echo "All done! Have an awesome day!"

