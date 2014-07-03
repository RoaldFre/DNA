% Parse all runs to an array. All runs must have the same starting time and 
% duration!
function [bound, time, state] = parseHairpinStateToArray(filesglob, alsoSaveFullState, dropLogFactor, DLdataStartIndex);

addpath generic

if nargin < 2; alsoSaveFullState = false; end
if nargin < 3; dropLogFactor = -1; end
if nargin < 4; DLdataStartIndex = 1; end

if ischar(filesglob)
	files = glob(filesglob);
elseif iscell(filesglob)
	files = cellfun(@glob,  filesglob, 'UniformOutput', false);
	files = flattenCell(files);
else
	error "I don't understand the filesglob format!"
end

if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

% read first file for the time and the number of samples
load(files{1});
fullTime = data(:, 1);
time = [fullTime(1:DLdataStartIndex-1); dropDataLogspace(fullTime(DLdataStartIndex:end), dropLogFactor)];
nSamples = numel(time);
nBonds = size(data)(2) - 2;

bound = zeros(nRuns, nSamples);
if alsoSaveFullState
	state = false(nRuns, nBonds, nSamples);
else
	state = [];
end

moreWasOn = page_screen_output;
more off
loadedRuns = false(nRuns,1);
for run = 1:nRuns
	try
		printf("\rReading file %d of %d", run, nRuns)
		load(files{run});

		thisTime = data(:,1);
		if not(equalsEpsilon(fullTime, thisTime, 1e-5))
			printf("\nError reading file %s\n", files{run})
			error "Reading datasets with different number of samples or time samples!"
		end

		thisBound = data(:, 2)';
		thisBound = [thisBound(1:DLdataStartIndex-1), dropDataLogspace(thisBound(DLdataStartIndex:end), dropLogFactor)];
	catch
		% Notify user of the error with that file, but carry on reading the rest
		printf("\nError reading file %s\n", files{run})
		continue;
	end
	
	loadedRuns(run) = true;
	bound(run,:) = thisBound;

	if alsoSaveFullState
		thisState = logical(data(:,3:end))';
		thisState = [thisState(1:DLdataStartIndex-1,:); dropDataLogspace(thisState(DLdataStartIndex:end,:), dropLogFactor)];
		state(run,:,:) = thisState;
	end
end
printf("\n");
if moreWasOn
	more on;
end

bound = bound(loadedRuns, :, :);
if alsoSaveFullState
	state = state(loadedRuns, :, :);
end
