#!/bin/bash

# Check whether the bases of the given world config file are not too close 
# (not within the repulsive part of the base pairing interaction -- more 
# specifically: they don't have a positive base pairing energy)

worldConfigFile=$1

tmpdir=`mktemp -d hasSafeBaseDist.XXXXXX --tmpdir=/tmp`

if ! hairpin -d "$worldConfigFile" -t '1e-20' -I '1e-23' -P '1e-26' -p -D "$tmpdir/data" > /dev/null
then
	echo "The file '$worldConfigFile' had an error! Returning failure!"
	exit 2
fi


scriptsdir="$HOME/DNA/scripts"
status=`octave -q -p $scriptsdir --eval "hasSafeBaseDistance('$tmpdir/data_basePairingEnergy')" | grep -v warning`

rm -rf $tmpdir

#echo $status
exit $status
