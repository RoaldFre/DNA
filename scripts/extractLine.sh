#!/bin/bash

# Extract a line from the cached resultfiles
# Currently hardcoded to make the first column of T=10 a 2-row multirow, 
# and the first column of all other T's (eg T=30) a blank space (to work 
# with the multirow of T=10)

name=$1
T=$2
scalingType=$3 #withTime or withSize
sqdev=$4 # 1 or 0
Ns=$5
meanStartN=$6

loglog=1
SA=0
resultsdir=$HOME/DNA/scripts/results

extractField() {
	field="$1"
	field=`printf '%q' "$field"` # escape special characters
	#field=`printf '%q' "$field"` # escape special characters
	singleExponent=$2

	file="$resultsdir/${name}T${T}_result${Ns}meanStartN${meanStartN}_${scalingType}_singexp${singleExponent}_2times0_sqdev${sqdev}_loglog${loglog}_SA${SA}"

	if ! sed -n "s#^${field}:##p" $file
	then
		#echo "--"
		#echo "*"
		echo "*"
		echo -e "ERROR: couldn't extract field '$field' from file $file" >&2
	fi
}


if [ "$scalingType" == withTime ]
then
	observable=`extractField 'observable' 0`

	betaSing=`extractField '\beta' 0`
	#deltaSing=`extractField 'theoreticalDelta' 0`
	qualSing=`extractField 'qual' 0`

	betaFull=`extractField '\beta' 1`
	#alphaFull=`extractField '\alpha' 1`
	alphaInvFull=`extractField '1/\alpha' 1`
	qualFull=`extractField 'qual' 1`

	#pmstring: string with original \pm plusminus
	pmString="$betaSing & $qualSing & $betaFull & $alphaInvFull & $qualFull"
elif [ "$scalingType" == withSize ]
then
	observable=`extractField 'observable' 0`

	betaSing=`extractField '\beta' 0`
	deltaSing=`extractField 'theoreticalDelta' 0`
	qualSing=`extractField 'qual' 0`

	betaFull=`extractField '\beta' 1`
	deltaFull=`extractField '\delta' 1`
	qualFull=`extractField 'qual' 1`

	#pmstring: string with original \pm plusminus
	pmString="$betaSing & \$$deltaSing\$ & $qualSing & $betaFull & $deltaFull & $qualFull"
else
	echo "UNKNOWN SCALING TYPE!"
	exit 1
fi

pmString="$T\$^\\circ\$C & $pmString"
if [ $T -eq 10 ]
then
	pmString="\\multirow{2}{*}{\$$observable\$} & $pmString"
else
	pmString=" & $pmString"
fi

echo "$pmString" | sed -e 's/ *\\pm */\&/g'
