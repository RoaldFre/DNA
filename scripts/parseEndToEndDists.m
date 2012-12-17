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
N = numel(sequence);
dists = zeros(numel(time), nRuns);

for run = 1:nRuns
	data = load(files{run});
	thistime = data(:,1);
	[_, thisSequence] = system(['sed -n "0,/# Strand1/s/# Strand 1: \(\w\)/\1/p" ',files{run}]);
	if (!isequal(thistime, fullTime) || thisSequence != sequence)
		error "Trying to load incompatible datasets with different times or sequence!"
	end
	dists(:, run) = decimateWrapper(data(:,2), decimateFactor);
end

