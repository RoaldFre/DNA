% Returns the longest possibe data set where all runs still have data
function [time, bound, meanState] = getIntersection(dataCell, N)

if nargin < 2
	N = numel(dataCell{1,3}(1,:));
end

numRuns = size(dataCell)(1);

numSamples = Inf;
for r = 1:numRuns
	numSamples = min(numSamples, numel(dataCell{r,1}));
end

time = dataCell{1,1}(1 : numSamples) - dataCell{1,1}(1);
bound = zeros(numRuns, numSamples);
meanState = zeros(numSamples, N);
for r = 1:numRuns
	bound(r, :) = dataCell{r,2}(1 : numSamples);
	meanState += dataCell{r,3}(1 : numSamples, 1 : N);
end
meanState /= numRuns;

