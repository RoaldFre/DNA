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
octave -q -p $scriptsdir --eval "toNative('$file')" &
PID=$!

sleep 0.5s

prevTicks=0

while true # keep polling CPU usage
do


	# DIT GEEF TOCH GEEN VOLLEDIGE CHILD TIME ?????
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat
	cut -d ' ' -f 14,15,16,17 /proc/$PID/stat

	cpuPerMil=`ps --cumulative S h o 'cp' p $PID`
	# I can't seem to get the total CPU time (including children [gzip]) in here :X
	if [ -z "$cpuPerMil" ]
	then
		echo "Couldn't find octave (cpuPerMil='$cpuPerMil'). Exiting cleanly! :-)"
		# We couldn't find octave, so it stopped cleanly!
		exit 0
	fi

	echo "Still running!"

	# Still running! Check if we are actually using CPU 
	# cycles, or if we are stuck
	echo "cpuPerMil: $cpuPerMil"
	if [ $cpuPerMil -lt 10 ] # less than 1% cpu average usage
	then
		echo "Seems stuck! :-( KILLING!!!!!"
		# Seems stuck: KILL, KILL, KILL!
		kill $PID
		echo "Sent SIGTERM, waiting..."
		sleep 5s
		kill -9 $PID
		echo "Sent SIGKILL, waiting..."
		sleep 5s

		echo "Recursive call to try again!"
		# Call ourselves again, then exit cleanly
		#$scriptsdir/toNative2.sh "$file"
		./$0 "$file"
		exit 0
	fi

	echo "Still running: sleeping for a bit"
	#Sleep some more.
	sleep 0.5s
done

