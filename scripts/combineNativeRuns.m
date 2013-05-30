% Read all (2D) data from the native files in the filesglob, and combine it 
% into a 3D array, where the first index is the run index. The actual data 
% of the native file must be in a variable named 'data'.
%
% The first column of each individual data set is regarded as a 'time' 
% column. It should be idential for all data sets, and it is returned 
% separately (and only once).
% If there is only one remaining column in the data file, the returned combined matrix is 2D instead of 3D.
%
% function [time, combinedData] = combineNativeRuns(filesglob, singlePrecision);
function [time, combinedData] = combineNativeRuns(filesglob, singlePrecision);

if nargin < 2
	singlePrecision = false;
end

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

% read first file for the time and the number of samples
load(files{1});
time = data(:, 1);
nSamples = numel(time);
nColumns = size(data,2) - 1;

if singlePrecision
	precisionString = "single";
else
	precisionString = "double";
end
combinedData = zeros(nRuns, nColumns, nSamples, precisionString);

moreWasOn = page_screen_output;
more off
for run = 1:nRuns
	printf("\rReading file %d of %d", run, nRuns)
	load(files{run});
	thisTime = data(:, 1);
	if not(isequal(time, thisTime))
		error "Reading datasets with different number of samples or time samples!"
	end
	combinedData(run,:,:) = data(:, 2:end);
end
printf("\n");
if moreWasOn
	more on;
end

combinedData = squeeze(combinedData);
