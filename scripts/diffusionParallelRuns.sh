destinationDir=$1
gamma=$2
N=$3
nRuns=$4

for run in `seq 1 $nRuns`
do
	./diffusionSingleRun.sh $destinationDir $gamma $N $run &
done

wait

