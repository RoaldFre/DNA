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
	octave -q -p $scriptsdir --eval "toNative('$file')" &
	PID=$!
	sleep 1m
	
	if ps -p $PID
	then
		#still running!
		kill $PID
		sleep 10s
		kill -9 $PID
		sleep 10s
		#next iteration in infinite loop
	else
		exit 0 #We can stop
	fi
done

