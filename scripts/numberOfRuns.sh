#!/bin/bash

dir="/home/other/roald/clusterdata_properVelocityInit_bpRepulsiveCORRECT/"
#dir="/home/other/roald/clusterdata_properVelocityInit"

cd $dir
for subdir in `find -type d | grep "NXY4$" | sort`
do
	echo $subdir
	for pattern in "*.itf11" "*iks" "*_basePairingEnergy" "*_endToEnd" "*_gyrationRadius" "*_forceVelFric"
	#for pattern in "*_forceVelFric" "*_forceVelFricP"
	do
		echo -e "`find $subdir/$pattern | wc -l`\t$pattern"
	done
done
