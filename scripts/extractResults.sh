#!/bin/bash

meanStartN=20
Ns=5.6.7.8.10.12.14.17.20.24.28.34.40.48

extractLine="$HOME/DNA/scripts/extractLine.sh"



echo 
echo 
scalingType=withTime
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

echo "`$extractLine state 10 $scalingType 0 $Ns $meanStartN`\\\\"
echo "`$extractLine state 30 $scalingType 0 $Ns $meanStartN`\\\\\\hline"
echo "`$extractLine energy 10 $scalingType 1 $Ns $meanStartN`\\\\"
echo "`$extractLine energy 30 $scalingType 1 $Ns $meanStartN`\\\\\\hline"
echo "`$extractLine gyrad 10 $scalingType 1 $Ns $meanStartN`\\\\"
echo "`$extractLine gyrad 30 $scalingType 1 $Ns $meanStartN`\\\\\\hline"
echo "`$extractLine endToEnd 10 $scalingType 1 $Ns $meanStartN`\\\\"
echo "`$extractLine endToEnd 30 $scalingType 1 $Ns $meanStartN`\\\\"

echo    '\end{tabular}'



echo 
echo 
echo 


scalingType=withSize
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

echo "`$extractLine state 10 $scalingType 0 $Ns $meanStartN`\\\\"
echo "`$extractLine state 30 $scalingType 0 $Ns $meanStartN`\\\\\\hline"
echo "`$extractLine energy 10 $scalingType 1 $Ns $meanStartN`\\\\"
echo "`$extractLine energy 30 $scalingType 1 $Ns $meanStartN`\\\\\\hline"
echo "`$extractLine gyrad 10 $scalingType 0 $Ns $meanStartN`\\\\"
echo "`$extractLine gyrad 30 $scalingType 0 $Ns $meanStartN`\\\\\\hline"
echo "`$extractLine gyrad 10 $scalingType 1 $Ns $meanStartN`\\\\"
echo "`$extractLine gyrad 30 $scalingType 1 $Ns $meanStartN`\\\\\\hline"
echo "`$extractLine endToEnd 10 $scalingType 1 $Ns $meanStartN`\\\\"
echo "`$extractLine endToEnd 30 $scalingType 1 $Ns $meanStartN`\\\\"

echo    '\end{tabular}'
