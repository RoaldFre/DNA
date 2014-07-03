% xs ys zs: every column is one run
% function [time, dists, N, sequence] = parseEndToEndVects(filesglob)
function [time, xs, ys, zs, dists, N, sequence] = parseEndToEndVects(filesglob)


files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

time = load(files{1})(:,1);
[_, sequence] = system(['sed -n "s/# Strand 1: \(\w\)/\1/p" ',files{1}]);
N = numel(sequence);
dists = zeros(numel(time), nRuns);
xs = zeros(numel(time), nRuns);
ys = zeros(numel(time), nRuns);
zs = zeros(numel(time), nRuns);

for run = 1:nRuns
	data = load(files{run});
	thistime = data(:,1);
	[_, thisSequence] = system(['sed -n "0,/# Strand1/s/# Strand 1: \(\w\)/\1/p" ',files{run}]);
	if (thisSequence != sequence)
		error "Trying to load incompatible datasets with different sequence!"
	elseif (!isequal(thistime, time))
		if numel(thistime) != numel(time)
			error "Trying to load incompatible datasets with different number of samples!"
		end
		% times are not exactly the same, but it's 'QK' if they are 
		% off by only a small fraction of the sample interval
		sampleInterval = time(3) - time(2);
		if (sum(thistime - time > sampleInterval / 100) > 1)
			% Too much error
			error "Trying to load incompatible datasets with different times!"
		end
	end
	dists(:, run) = data(:,2);
	xs(:, run) = data(:,3);
	ys(:, run) = data(:,4);
	zs(:, run) = data(:,5);
end

