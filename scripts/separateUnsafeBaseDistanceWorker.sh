#!/bin/bash

dir=$1
unsafeDir=$2
base=$3

scriptsdir="$HOME/DNA/scripts"
checkCmd="$scriptsdir/hasSafeBaseDistance.sh"

#echo "************ WORKER ***********:"
#echo "dir       $dir"
#echo "unsafeDir $unsafeDir"
#echo "base      $base"

if ! [ -e "${base}_finalRelax" ]
then
	echo "$base has no finalRelax file! Not moving!"
	exit
fi

if "$checkCmd" "${base}_finalRelax"
then
	echo "OK        $base"
	mv "$base"* $dir
else
	echo "**WRONG** $base"
	mv "$base"* $unsafeDir
fi

