% function [time, samplerData, bound, correctlyBound] = parseBasePairingSamplerToArray(filesglob);
%
% Parse all runs to an array. All runs must have the same starting time and 
% duration!
%
% The data files should be made with the basePairingSampler.
% The 'data' array holds the data of the sampler (being energy, distances, 
% or boolean 'bound' values, depending on which sampler mode was used).
function [time, samplerData, bound, correctlyBound] = parseBasePairingSamplerToArray(filesglob, dropLogFactor);

addpath generic

if nargin < 2
	dropLogFactor = -1;
end

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

% read first file for the time and the number of samples
load(files{1});
nBonds = size(data)(2) - 3;
fullTime = data(:,nBonds + 1);
time = dropDataLogspace(fullTime, dropLogFactor);
nSamples = numel(time);

bound = zeros(nRuns, nSamples);
correctlyBound = zeros(nRuns, nSamples);
%samplerData = zeros(nRuns, nBonds, nSamples);
samplerData = zeros(nRuns, nBonds, nSamples, "single");

moreWasOn = page_screen_output;
more off
loadedRuns = false(nRuns,1);
for run = 1:nRuns
	try
		printf("\rReading file %d of %d", run, nRuns)
		load(files{run});
		thisTime = data(:, nBonds + 1);
		if not(equalsEpsilon(fullTime, thisTime, 1e-5))
			printf("\nError reading file %s\n", files{run})
			error "Reading datasets with different number of samples or time samples!"
		end
	catch
		printf("\nError reading file %s\n", files{run})
		continue
	end

	loadedRuns(run) = true;
	bound(run,:) = dropDataLogspace(data(:, nBonds + 2)', dropLogFactor);
	correctlyBound(run,:) = dropDataLogspace(data(:, nBonds + 3)', dropLogFactor);
	samplerData(run,:,:) = dropDataLogspace(data(:,1:nBonds)', dropLogFactor);
end
printf("\n");
if moreWasOn
	more on;
end

bound = bound(loadedRuns,:);
correctlyBound = correctlyBound(loadedRuns,:);
samplerData = samplerData(loadedRuns,:,:);
