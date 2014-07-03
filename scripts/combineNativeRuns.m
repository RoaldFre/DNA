% Read all (2D) data from the native files in the filesglob, and combine it 
% into a 3D array, where the first index is the run index. The actual data 
% of the native file must be in a variable named 'data'.
%
% The first column of each individual data set is regarded as a 'time' 
% column. It should be idential for all data sets, and it is returned 
% separately (and only once).
% If there is only one remaining column in the data file, the returned combined matrix is 2D instead of 3D.
% The optional 'columns' argument should contain a list of column indices 
% that should be extracted instead of extracting the full data. Index 1 
% means the first data column (so excluding the time column in the 
% indexing).
%
% function [time, combinedData] = combineNativeRuns(filesglob, singlePrecision, columns);
function [time, combinedData, files] = combineNativeRuns(filesglob, singlePrecision, dropLogFactor, DLdataStartIndex, columns);

addpath generic

if nargin < 2; singlePrecision = false; end
if nargin < 3; dropLogFactor = -1; end
if nargin < 4; DLdataStartIndex = 1; end
allColumns = nargin < 5;

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

% read first file for the time and the number of samples
load(files{1});
if exist("data") != 1
	error(["Couldn't find a 'data' variable in ",files{1}]);
end
fullTime = data(:, 1);
time = [fullTime(1:DLdataStartIndex-1); dropDataLogspace(fullTime(DLdataStartIndex:end), dropLogFactor)];
nSamples = numel(time);

if allColumns
	nColumns = size(data,2) - 1;
	columns = 1:nColumns;
else
	nColumns = numel(columns);
end

if singlePrecision
	precisionString = "single";
else
	precisionString = "double";
end
combinedData = zeros(nRuns, nSamples, nColumns, precisionString);

moreWasOn = page_screen_output;
more off
loadedRuns = false(nRuns,1);
for run = 1:nRuns
	try
		printf("\rReading file %d of %d", run, nRuns)
		load(files{run});
		if exist("data") != 1
			error(["Couldn't find a 'data' variable in ",files{run}]);
		end
		thisTime = data(:, 1);
		if not(equalsEpsilon(fullTime, thisTime, 1e-5))
			printf("\nError reading file %s\n", files{run})
			error "Reading datasets with different number of samples or time samples!"
		end
		thisData = data(:, 1 + columns);
		thisData = [thisData(1:DLdataStartIndex-1,:); dropDataLogspace(thisData(DLdataStartIndex:end,:)', dropLogFactor)'];
	catch
		% Notify user of the error with that file, but carry on reading the rest
		printf("\nError reading file %s\n", files{run})
		continue;
	end

	loadedRuns(run) = true;
	combinedData(run,:,:) = thisData;
end
printf("\n");
if moreWasOn
	more on;
end

combinedData = squeeze(combinedData(loadedRuns, :, :));
files = files{loadedRuns};
