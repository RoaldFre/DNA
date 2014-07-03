% This is simila rto parseBasePairingSamplerToArray, except that it already 
% performs a sum over the bases with indices given in basesToSum.
%
% If basesToSum is [] or omitted, then all bases are summed over (WARNING: 
% these can include XY pairs!)
%
% The data files should be made with the basePairingSampler and converted 
% to native Octave format with toNative.sh.
function [time, summedSamplerData, bound, correctlyBound] = parseBasePairingSamplerSumToArray(filesglob, dropLogFactor, DLdataStartIndex, basesToSum, minBoundThreshold);


if nargin < 4 || isempty(basesToSum)
	basesToSum = 1:nBonds;
end
if nargin < 5
	minBoundThreshold = 0;
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
time = [fullTime(1:DLdataStartIndex-1); dropDataLogspace(fullTime(DLdataStartIndex:end), dropLogFactor)];
nSamples = numel(time);



bound = zeros(nRuns, nSamples);
correctlyBound = zeros(nRuns, nSamples);
summedSamplerData = zeros(nRuns, nSamples);

moreWasOn = page_screen_output;
more off
nGoodRuns = 0; % number of accepted files (can be less than nRuns because of minBoundThreshold)
for run = 1:nRuns
	try
		printf("\rReading file %d of %d", run, nRuns)
		load(files{run});
		thisTime = data(:, nBonds + 1);
		if not(equalsEpsilon(fullTime, thisTime, 1e-5))
			printf("\nError reading file %s\n", files{run})
			error "Reading datasets with different number of samples or time samples!"
		end
		thisBound = data(:, nBonds + 2)';
		if sum(thisBound < minBoundThreshold) > 0
			continue % This run doesn't count.
		end
	catch
		printf("\nError reading file %s\n", files{run})
		continue
	end

	nGoodRuns++;
	thisCB = data(:, nBonds + 3)';;
	bound(nGoodRuns,:) = [thisBound(1:DLdataStartIndex-1), dropDataLogspace(thisBound(DLdataStartIndex:end), dropLogFactor)];
	correctlyBound(nGoodRuns,:) = [thisCB(1:DLdataStartIndex-1), dropDataLogspace(thisCB(DLdataStartIndex:end), dropLogFactor)];
	thisData = data(:,basesToSum)';
	thisData = [thisData(:,1:DLdataStartIndex-1), dropDataLogspace(thisData(:,DLdataStartIndex:end), dropLogFactor)];
	summedSamplerData(nGoodRuns,:) = sum(thisData, 1);
end
printf("\n");
if moreWasOn
	more on;
end

runsVSgoodRuns = [nRuns, nGoodRuns]
bound = bound(1:nGoodRuns, :);
correctlyBound = correctlyBound(1:nGoodRuns, :);
summedSamplerData = summedSamplerData(1:nGoodRuns, :);

