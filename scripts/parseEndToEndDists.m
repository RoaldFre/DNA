% dists: every column is one run
function [time, dists] = parseEndToEndDists(filesglob, decimateFactor)

addpath('octave-forge');

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
time = decimate(fullTime, decimateFactor);
dists = zeros(numel(time), nRuns);

for run = 1:nRuns
	data = load(files{run});
	thistime = data(:,1);
	if (!isequal(thistime, fullTime))
		error "Trying to load incompatible datasets with different times!"
	end
	dists(:, run) = decimate(data(:,2), decimateFactor);
end

