#!/bin/bash

file=$1
scriptsdir="$HOME/DNA/scripts"

if [ `file $file -b --mime-type` != 'text/plain' ]
then
	echo "The file '$file' is not a plain text file! Bailing out!"
	exit 1
fi

# Ok, this is really really REALLY ugly, but it looks like sometimes octave 
# gets stuck. So we spawn it off in the background, wait for a sufficiently 
# long time and if it still hasn't completed yet -> kill it and try again.
while true
do
	#octave -q -p $scriptsdir --eval "toNative('$file')" &
	octave-3.6.4 -q -p $scriptsdir --eval "toNative('$file')" &
	PID=$!

	sleep 1s

	for sleeptime in 1s 2s 4s 8s 15s 30s 1m 2m 3m 4m 5m 10m
	do
		if ps -p $PID
		then
			#still running! Sleep some more.
			sleep 1m
		else
			exit 0
		fi
	done

	
	if ps -p $PID
	then
		#still running! KILL, KILL, KILL!
		kill $PID
		sleep 5s
		kill -9 $PID
		sleep 5s
		#next iteration in infinite loop
	else
		exit 0 # It stopped nicely after waiting for the last time!
	fi
done

