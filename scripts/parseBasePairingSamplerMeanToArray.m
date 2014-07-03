% Parse the mean and error of all runs to an array. All runs must have the 
% same starting time and duration!
%
% The sum of the data over the bases given in basesToSum is returned in the 
% *SummedData output variables.
% If basesToSum is [] or omitted, then all bases are summed over (WARNING: 
% these can include XY pairs!)
%
% The *SummedTranslatedSq variants are translated so that they start at 0 
% and are squared. The error takes into account the covariance between the 
% sums at time t=0 and at later times.
%
% The data files should be made with the basePairingSampler and converted 
% to native Octave format with toNative.sh.
%
% This does two passes over the files, in order to limit memory usage.
function [time, meanSamplerData, errSamplerData, meanOfSummedData, errOfSummedData, meanOfSummedTranslatedSq, errOfSummedTranslatedSq] = parseBasePairingSamplerMeanToArray(filesglob, basesToSum, minBoundThreshold);


files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end
nRuns = numel(files);
% read first file for the time and the number of samples
load(files{1});
nBonds = size(data)(2) - 3;
time = data(:,nBonds + 1);
nSamples = numel(time);



if nargin < 2 || isempty(basesToSum)
	basesToSum = 1:nBonds;
end
if nargin < 3
	minBoundThreshold = 0;
end



moreWasOn = page_screen_output;
more off

% Make a first pass to collect the mean
nGoodRuns = 0; % number of accepted files (can be less than nRuns because of minBoundThreshold)
accumulatedData = zeros(nSamples, nBonds);
accumulatedSummedData = zeros(nSamples, 1); % summed over all bases (and accumulation == sum over runs)
for run = 1:nRuns
	try
		printf("\rPass 1: Reading file %d of %d", run, nRuns)
		load(files{run});
		thisTime = data(:, nBonds + 1);
		if not(equalsEpsilon(fullTime, thisTime, 1e-5))
			printf("\nError reading file %s\n", files{run})
			error "Reading datasets with different number of samples or time samples!"
		end
		bound = data(:, nBonds + 2)';
		%correctlyBound = data(:, nBonds + 3)';
		if sum(bound < minBoundThreshold) > 0
			continue % This run doesn't count.
		end
	catch
		printf("\nError reading file %s\n", files{run})
		% just go to next iteration
		continue
	end

	nGoodRuns++;
	accumulatedData += data(:,1:nBonds);
	accumulatedSummedData += sum(data(:, basesToSum), 2);
end
printf("\nRead %d files, of which %d were accepted\n", nRuns, nGoodRuns);

meanSamplerData = accumulatedData / nGoodRuns;
meanOfSummedData = accumulatedSummedData / nGoodRuns;
% reset the accumulators to get the variance
accumulatedData = zeros(nSamples, nBonds);
accumulatedSummedData = zeros(nSamples, 1);
accumulatedCovarData = zeros(nSamples, 1); % covariance with t=0 and later times

% Make a second pass to collect the error
for run = nRuns:-1:1 % run in reverse to make maximum use of disk cache from previous run
	printf("\rPass 2: Reading file %d of %d    ", run, nRuns)
	load(files{run});
	bound = data(:, nBonds + 2)';
	if sum(bound < minBoundThreshold) > 0
		continue % This run doesn't count.
	end

	accumulatedData += (data(:,1:nBonds) - meanSamplerData).^2;
	accumulatedSummedData += (sum(data(:, basesToSum), 2) - meanOfSummedData).^2;
	accumulatedCovarData += (sum(data(1, basesToSum), 2) - meanOfSummedData(1)) ...
                              * (sum(data(:, basesToSum), 2) - meanOfSummedData   );

end
printf("\n");


errSamplerData = sqrt(accumulatedData) / (nGoodRuns - 1);
errOfSummedData = sqrt(accumulatedSummedData) / (nGoodRuns - 1);
varOfSummedData = accumulatedSummedData / (nGoodRuns - 1);
covarWithFirst = accumulatedCovarData / (nGoodRuns - 1);

meanOfSummedTranslatedSq = (meanOfSummedData - meanOfSummedData(1)).^2;
% For error propagation:
%   |d((x - y)^2)| = 2 |x - y| |d(x - y)|
%   for the error on (x - y):
%   var(x - y) = var(x) + var(y) - 2*covar(x, y)
% here x = meanOfSummedData(t)
%      y = meanOfSummedData(0)
errOfSummedTranslatedSq = 2 * abs(meanOfSummedData - meanOfSummedData(1)) .* sqrt(max(0, varOfSummedData + varOfSummedData(1) - 2*covarWithFirst)) / sqrt(nGoodRuns - 1);


if moreWasOn
	more on;
end

