gamma=$1
N=$2
nRuns=$3

for run in `seq 1 $nRuns`
do
	./diffusionSingleRun.sh $gamma $N $run &
done

wait

