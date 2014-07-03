%Input: times = cell array of times from each run
%       numBounds = cell array of numBound from each run
%OR: input can be 1D array for time and 2D array for numBounds instead of cells if all runs have the same length. The 2D array should have each run in one column
%OR: input can be 1D arrays like for getZippingTimesFromBound(), in which case this 
%function behaves identically
function [zippingTimes, nucleationTimes, zipFromNuclTimes, zipIndices, nucleationIndices] = getZippingTimesFromBounds(times, numBounds, nuclBoundBPs, zippedBoundBPs)

% If we get data that isn't in cell format:
if not(iscell(times))
	if numel(times) == numel(numBounds)
		% Only one run
		[zippingTimes, nucleationTimes, zipFromNuclTimes, zipIndices, nucleationIndices] = getZippingTimesFromBound(times, numBounds, nuclBoundBPs, zippedBoundBPs);
		return
	else
		% Multiple runs. Convert to cell format so we can use code below
		Nruns = size(numBounds)(2);
		cellTimes = cell(Nruns, 1); % Yes, this is wasteful...
		cellBounds = cell(Nruns, 1);
		for i = 1:Nruns
			cellTimes{i} = times;
			cellBounds{i} = numBounds(:, i);
		end
		times = cellTimes;
		numBounds = cellBounds;
	end
end


N = numel(times);
if (numel(numBounds) != N)
	error "times and numBounds don't have the same number of elements"
end

zippingTimes = zeros(N, 1);
nucleationTimes = zeros(N, 1);
zipFromNuclTimes = zeros(N, 1);
zipIndices = zeros(N, 1);
nucleationIndices = zeros(N, 1);

for i = 1:N
	[zippingTime, nucleationTime, zipFromNuclTime, zipIndex, nucleationIndex] = getZippingTimesFromBound(times{i}, numBounds{i}, nuclBoundBPs, zippedBoundBPs);
	zippingTimes(i) = zippingTime;
	nucleationTimes(i) = nucleationTime;
	zipFromNuclTimes(i) = zipFromNuclTime;
	zipIndices(i) = zipIndex;
	nucleationIndices(i) = nucleationIndex;
end
