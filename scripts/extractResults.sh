#!/bin/bash

# We need the 'multirow' and 'pdflscape' package loaded for our LaTeX code to be able to compile!

meanStartN=20
#Ns=5.6.7.8.10.12.14.17.20.24.28.34.40.48
#Ns=20.24.28.34.40.48
#Ns=20.24.28.34.40
Ns=20.22.24.26.28.31.34.37.40.44.48
loglog=1
SA=0
#dropLog=3
dropLog=5
clustSize=2
#clustSize=3
#suff=NoDih1.7
suff=''

resultsdir=$HOME/DNA/scripts/results
dummyFile="$resultsdir/dummy"
resultPlotsFile=/tmp/resultPlots
collapsePlotsFile=/tmp/collapsePlots
plotInputDir=images/plots

IFS='.' read -a NsArr <<< "$Ns"
numNs=${#NsArr[@]}

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

	if ! [ -e $file ]
	then
		file="$dummyFile"
	fi

	if ! sed -n "s#^${field}:##p" $file
	then
		#echo "--"
		#echo "*"
		echo "*"
		echo -e "ERROR: couldn't extract field '$field' from file $file" >&2
	fi
}

makeResultPlotRef() {
	# This also dumps the input command to the resultPlotsFile
	singleExponent=$1

	plotFile="${name}_T${T}_clustSize${clustSize}_dropLog${dropLog}_results_`makeSuffix $singleExponent $sqDev`"
	echo "\\input{$plotInputDir/${plotFile}}" >> $resultPlotsFile
	echo "\\ref{$plotFile}"
}

writeCollapsePlotFigs() {
	singleExponent=$1

	for i in `seq 0 $(( $numNs - $clustSize ))`
	do
		clustNs="${NsArr[@]:$i:$clustSize}" # sub array in bash
		clustNs=`echo $clustNs | tr ' ' '.'`

		plotFile="${name}_T${T}_clustSize${clustSize}_dropLog${dropLog}_collap${clustNs}_`makeSuffix $singleExponent $sqDev`"
		echo "\\input{$plotInputDir/${plotFile}}" >> $collapsePlotsFile
	done
	echo "\\FloatBarrier" >> $collapsePlotsFile
}

extractLine() {
	# Extract a line from the cached resultfiles
	# Currently hardcoded to make the first column of T=10 a 2-row 
	# multirow, and the first column of all other T's (eg T=30) a blank 
	# space (to work with the multirow of T=10)
	#
	# This also adds relevant \input{}'s to the $resultPlotsFile and 
	# references those plots.
	# Lastly, \input{}'s for the individual clustered collapses are 
	# written to $collapsePlotsFile.

	name=$1
	T=$2
	sqDev=$3 # 1 or 0

	if [ "$scalingType" == withTime ]
	then
		observable=`extractField 'observable' 1 $sqDev`

		betaSing=`extractField '\beta' 1 $sqDev`
		#deltaSing=`extractField 'theoreticalDelta' 1 $sqDev`
		qualSing=`extractField 'qual' 1 $sqDev`
		plotRefSing=`makeResultPlotRef 1`

		betaFull=`extractField '\beta' 0 $sqDev`
		#alphaFull=`extractField '\alpha' 0 $sqDev`
		alphaInvFull=`extractField '1/\alpha' 0 $sqDev`
		qualFull=`extractField 'qual' 0 $sqDev`
		plotRefFull=`makeResultPlotRef 0`

		#pmstring: string with original \pm plusminus
		pmString="$betaSing & $qualSing & $plotRefSing & $betaFull & $alphaInvFull & $qualFull & $plotRefFull"
	elif [ "$scalingType" == withSize ]
	then
		observable=`extractField 'observable' 1 $sqDev`

		betaSing=`extractField '\beta' 1 $sqDev`
		deltaSing=`extractField 'theoreticalDelta' 1 $sqDev`
		qualSing=`extractField 'qual' 1 $sqDev`
		plotRefSing=`makeResultPlotRef 1`

		betaFull=`extractField '\beta' 0 $sqDev`
		deltaFull=`extractField '\delta' 0 $sqDev`
		qualFull=`extractField 'qual' 0 $sqDev`
		plotRefFull=`makeResultPlotRef 0`

		#pmstring: string with original \pm plusminus
		pmString="$betaSing & \$$deltaSing\$ & $qualSing & $plotRefSing & $betaFull & $deltaFull & $qualFull & $plotRefFull"
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

	writeCollapsePlotFigs 1
	writeCollapsePlotFigs 0
	echo "\\FloatBarrier" >> $resultPlotsFile
}


echo -n "" > $resultPlotsFile
echo -n "" > $collapsePlotsFile


echo 
echo 

scalingType=withTime # This is a global, also used in functions above!
echo '\begin{landscape}'
echo '\begin{table}[t]'
echo '\caption{Combined result of clustered collapse with clusters of size' "$clustSize," 'scaling with time $t$.}'
echo '\label{tblClusteredResultsWithTime}'
echo '\begin{center}'
echo '\begin{tabular}{cc || r@{$\pm$}l  r@{$\pm$}l c | r@{$\pm$}l  r@{$\pm$}l  r@{$\pm$}l c}'

echo -n '\multirow{2}{*}{observable} &'
echo -n '\multirow{2}{*}{$T$} &'
echo -n '\multicolumn{5}{|c}{fixed $\alpha = 1/\beta$} &'
echo    '\multicolumn{7}{|c}{fitted $\alpha$} \\'

echo -n ' & &'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '\multicolumn{2}{c}{$Q$} &'
echo -n 'Fig.&'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '\multicolumn{2}{c}{$1/\alpha$} &'
echo -n '\multicolumn{2}{c}{$Q$} &'
echo    'Fig.\\\hline\hline'

echo "`extractLine state 10 0`\\\\"
echo "`extractLine state 30 0`\\\\\\hline"
echo "`extractLine energy 10 1`\\\\"
echo "`extractLine energy 30 1`\\\\\\hline"
echo "`extractLine gyrad 10 1`\\\\"
echo "`extractLine gyrad 30 1`\\\\\\hline"
echo "`extractLine endToEnd 10 1`\\\\"
echo "`extractLine endToEnd 30 1`\\\\"

#echo "`extractLine state$suff 10 0`\\\\"
##echo "`extractLine state$suff 30 0`\\\\\\hline"
#echo "`extractLine energy$suff 10 1`\\\\"
##echo "`extractLine energy$suff 30 1`\\\\\\hline"
#echo "`extractLine gyrad$suff 10 0`\\\\"
##echo "`extractLine gyrad$suff 30 0`\\\\\\hline"
#echo "`extractLine gyrad$suff 10 1`\\\\"
##echo "`extractLine gyrad$suff 30 1`\\\\\\hline"
#echo "`extractLine endToEnd$suff 10 1`\\\\"
##echo "`extractLine endToEnd$suff 30 1`\\\\"

echo '\end{tabular}'
echo '\end{center}'
echo '\end{table}'
echo '\end{landscape}'



echo 
echo 
echo 


scalingType=withSize # This is a global, also used in functions above!
echo '\begin{landscape}'
echo '\begin{table}[b]'
echo '\caption{Combined result of clustered collapse with clusters of size' "$clustSize," 'scaling with size $N$.}'
echo '\label{tblClusteredResultsWithSize}'
echo '\begin{center}'
echo '\begin{tabular}{cc || r@{$\pm$}l c r@{$\pm$}l c | r@{$\pm$}l  r@{$\pm$}l  r@{$\pm$}l c}'

echo -n '\multirow{2}{*}{observable} &'
echo -n '\multirow{2}{*}{$T$} &'
echo -n '\multicolumn{6}{|c}{fixed $\delta$} &'
echo    '\multicolumn{7}{|c}{fitted $\delta$} \\'

echo -n ' & &'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '$\delta$ &'
echo -n '\multicolumn{2}{c}{$Q$} &'
echo -n 'Fig.&'
echo -n '\multicolumn{2}{c}{$\beta$} &'
echo -n '\multicolumn{2}{c}{$\delta$} &'
echo -n '\multicolumn{2}{c}{$Q$} &'
echo    'Fig.\\\hline\hline'

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

#echo "`extractLine state$suff 10 0`\\\\"
##echo "`extractLine state$suff 30 0`\\\\\\hline"
#echo "`extractLine energy$suff 10 1`\\\\"
##echo "`extractLine energy$suff 30 1`\\\\\\hline"
#echo "`extractLine gyrad$suff 10 0`\\\\"
##echo "`extractLine gyrad$suff 30 0`\\\\\\hline"
#echo "`extractLine gyrad$suff 10 1`\\\\"
##echo "`extractLine gyrad$suff 30 1`\\\\\\hline"
#echo "`extractLine endToEnd$suff 10 1`\\\\"
##echo "`extractLine endToEnd$suff 30 1`\\\\"

echo '\end{tabular}'
echo '\end{center}'
echo '\end{table}'
echo '\end{landscape}'

echo 
cat $resultPlotsFile
echo 
cat $collapsePlotsFile
echo 
