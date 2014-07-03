#/bin/bash

dir="$1"

numpar=10

scriptsdir="$HOME/DNA/scripts"
uncheckedDir="$dir/uncheckedBaseDistance"
unsafeDir="$dir/unsafeBaseDistance"

mkdir -p $uncheckedDir
mv $dir/* $uncheckedDir
#mv $uncheckedDir/$uncheckedDir/* $uncheckedDir # If we already tried this and aborted midway the move
mkdir -p $unsafeDir

echo $uncheckedDir/state*.itf11 | xargs -n 1 -P $numpar $scriptsdir/separateUnsafeBaseDistanceWorker.sh $dir $unsafeDir 


