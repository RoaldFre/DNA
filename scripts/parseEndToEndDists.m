% THIS IS OBSOLETE in favor combineNativeRuns() working on native data. 
% Nonetheless, this function can still be useful as it returns some extra 
% data (N and sequence) and does some extra consistency checks.
%
% dists: every column is one run
% function [time, dists, N, sequence] = parseEndToEndDists(filesglob, decimateFactor)
function [time, dists, N, sequence] = parseEndToEndDists(filesglob, decimateFactor)

if (nargin < 1)
	error "Not enough required arguments given!"
elseif (nargin < 2)
	decimateFactor = 1;
end


files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

fullTime = load(files{1})(:,1);
time = decimateWrapper(fullTime, decimateFactor);
[_, sequence] = system(['sed -n "s/# Strand 1: \(\w\)/\1/p" ',files{1}]);
N = numel(sequence) - 1; % trailing \n
dists = zeros(numel(time), nRuns);

moreWasOn = page_screen_output;
more off;
for run = 1:nRuns
	printf("\rLoading endToEnd data from file %d of %d", run, nRuns);
	data = load(files{run});
	thistime = data(:,1);
	[_, thisSequence] = system(['sed -n "0,/# Strand1/s/# Strand 1: \(\w\)/\1/p" ',files{run}]);
	if (thisSequence != sequence)
		error "Trying to load incompatible datasets with different sequence!"
	elseif (!isequal(thistime, fullTime))
		if numel(thistime) != numel(fullTime)
			error "Trying to load incompatible datasets with different number of samples!"
		end
		% times are not exactly the same, but it's 'QK' if they are 
		% off by only a small fraction of the sample interval
		sampleInterval = time(3) - time(2);
		%[thistime; fullTime; thistime-fullTime; thistime - fullTime > sampleInterval]
		if (sum(thistime - fullTime > sampleInterval / 100) > 1)
			% Too much error
			error "Trying to load incompatible datasets with different times!"
		end
	end
	dists(:, run) = decimateWrapper(data(:,2), decimateFactor);
end

printf("\n");
if moreWasOn
	more on;
end

