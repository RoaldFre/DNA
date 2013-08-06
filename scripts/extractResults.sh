#!/bin/bash

meanStartN=20
Ns=5.6.7.8.10.12.14.17.20.24.28.34.40.48
loglog=1
SA=0


resultsdir=$HOME/DNA/scripts/results

makeSuffix() {
	singleExponent=$1
	sqDev=$2
	echo "${scalingType}_singexp${singleExponent}_2times0_sqdev${sqDev}_loglog${loglog}_SA${SA}"
}

extractField() {
	field="$1"
	field=`printf '%q' "$field"` # escape special characters
	#field=`printf '%q' "$field"` # escape special characters
	singleExponent=$2
	sqDev=$3

	file="$resultsdir/${name}T${T}_result${Ns}meanStartN${meanStartN}_`makeSuffix $singleExponent $sqDev`"

	if ! sed -n "s#^${field}:##p" $file
	then
		#echo "--"
		#echo "*"
		echo "*"
		echo -e "ERROR: couldn't extract field '$field' from file $file" >&2
	fi
}

extractLine() {
	# Extract a line from the cached resultfiles
	# Currently hardcoded to make the first column of T=10 a 2-row 
	# multirow, and the first column of all other T's (eg T=30) a blank 
	# space (to work with the multirow of T=10)

	name=$1
	T=$2
	sqDev=$3 # 1 or 0

	if [ "$scalingType" == withTime ]
	then
		observable=`extractField 'observable' 0 $sqDev`

		betaSing=`extractField '\beta' 0 $sqDev`
		#deltaSing=`extractField 'theoreticalDelta' 0 $sqDev`
		qualSing=`extractField 'qual' 0 $sqDev`

		betaFull=`extractField '\beta' 1 $sqDev`
		#alphaFull=`extractField '\alpha' 1 $sqDev`
		alphaInvFull=`extractField '1/\alpha' 1 $sqDev`
		qualFull=`extractField 'qual' 1 $sqDev`

		#pmstring: string with original \pm plusminus
		pmString="$betaSing & $qualSing & $betaFull & $alphaInvFull & $qualFull"
	elif [ "$scalingType" == withSize ]
	then
		observable=`extractField 'observable' 0 $sqDev`

		betaSing=`extractField '\beta' 0 $sqDev`
		deltaSing=`extractField 'theoreticalDelta' 0 $sqDev`
		qualSing=`extractField 'qual' 0 $sqDev`

		betaFull=`extractField '\beta' 1 $sqDev`
		deltaFull=`extractField '\delta' 1 $sqDev`
		qualFull=`extractField 'qual' 1 $sqDev`

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
}



echo 
echo 
scalingType=withTime # This is a global, also used in functions above!
echo    '\begin{tabular}{cc | r@{$\pm$}l  r@{$\pm$}l | r@{$\pm$}l  r@{$\pm$}l  r@{$\pm$}l}'

echo -n '\multirow{2}{*}{observable} &'
echo -n '\multirow{2}{*}{$T$} &'
echo -n '\multicolumn{4}{|c}{fixed $\alpha = 1/\beta$} &'
echo    '\multicolumn{6}{|c}{fitted $\alpha$} \\'

echo -n ' & &'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '\multicolumn{2}{c|}{$Q$} &'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '\multicolumn{2}{c}{$1/\alpha$} &'
echo    '\multicolumn{2}{c}{$Q$} \\\hline'

echo "`extractLine state 10 0`\\\\"
echo "`extractLine state 30 0`\\\\\\hline"
echo "`extractLine energy 10 1`\\\\"
echo "`extractLine energy 30 1`\\\\\\hline"
echo "`extractLine gyrad 10 1`\\\\"
echo "`extractLine gyrad 30 1`\\\\\\hline"
echo "`extractLine endToEnd 10 1`\\\\"
echo "`extractLine endToEnd 30 1`\\\\"

echo    '\end{tabular}'



echo 
echo 
echo 


scalingType=withSize # This is a global, also used in functions above!
echo    '\begin{tabular}{cc | r@{$\pm$}l c r@{$\pm$}l | r@{$\pm$}l  r@{$\pm$}l  r@{$\pm$}l}'

echo -n '\multirow{2}{*}{observable} &'
echo -n '\multirow{2}{*}{$T$} &'
echo -n '\multicolumn{5}{|c}{fixed $\delta$} &'
echo    '\multicolumn{6}{|c}{fitted $\delta$} \\'

echo -n ' & &'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '$\delta$ &'
echo -n '\multicolumn{2}{c|}{$Q$} &'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '\multicolumn{2}{c}{$\delta$} &'
echo    '\multicolumn{2}{c}{$Q$} \\\hline\hline'

echo "`extractLine state 10 0`\\\\"
echo "`extractLine state 30 0`\\\\\\hline"
echo "`extractLine energy 10 1`\\\\"
echo "`extractLine energy 30 1`\\\\\\hline"
echo "`extractLine gyrad 10 0`\\\\"
echo "`extractLine gyrad 30 0`\\\\\\hline"
echo "`extractLine gyrad 10 1`\\\\"
echo "`extractLine gyrad 30 1`\\\\\\hline"
echo "`extractLine endToEnd 10 1`\\\\"
echo "`extractLine endToEnd 30 1`\\\\"

echo    '\end{tabular}'
