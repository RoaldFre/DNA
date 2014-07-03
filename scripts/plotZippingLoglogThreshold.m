% Plots a linear regression in the loglog plot corresponding to the points 
% where the numBound in the dataset crossed the given zippedBoundBPs 
% threshold.
%
% times: cell array of time data for the different lengths
% bounds: cell array of numbound data for the different lenghts. Can be the 
%         mean bound (1D array), or the numbound data for all different 
%         runs (2D array). In the latter case: each row represents one 
%         run and the runs are averaged here TODO: 
%         linear average or log average??
% color: color of the plot. If empty: don't plot.
function [zippingTimes, zippingBounds, zipE, zipC] = plotZippingLoglogThreshold(times, bounds, nuclBoundBPs, zippedBoundBPs, color,  regressionBeginIndex)

if nargin < 5
	color = 'b';
end
if nargin < 6
	regressionBeginIndex = 1;
end

N = numel(times);
if (numel(bounds) != N)
	error "times and bounds don't have the same number of elements"
end

zippingBounds = zeros(N, 1);
zippingTimes = zeros(N, 1);
zippingIndices = zeros(N, 1);
meanNumBound = cell(N, 1);

hold on;
%more off;
for i = 1:N
	nRuns = size(bounds{i})(1);

	if nRuns > 1
		% TODO: linear or log averaging?
		meanNumBound{i} = mean(bounds{i}());
		%meanNumBound{i} = exp(mean(log(bounds{i})));
		%meanNumBound{i} = exp(mean(log(bounds{i} + 1))) - 1;
	else
		meanNumBound{i} = bounds{i};
	end

	if not(isempty(color))
		%plot(time{i}, meanNumBound{i}, color);
		loglog(times{i}, meanNumBound{i}, color);
	end

	[zippingTime, nucleationTime, zipFromNuclTime, zipIndex, nucleationIndex] = getZippingTimesFromBounds(times{i}, meanNumBound{i}, nuclBoundBPs(i), zippedBoundBPs(i));
	boundAtThreshold = meanNumBound{i}(zipIndex);
	zippingTimes(i) = zippingTime;
	zippingIndices(i) = zipIndex;
	zippingBounds(i) = boundAtThreshold;
end

%zippingPoints

[zipCte, zipExponent] = loglogRegression(zippingTimes(regressionBeginIndex:end), zippingBounds(regressionBeginIndex:end), 1, 1);

if not(isempty(color))
	loglog(zippingTimes, zippingBounds, '+k');
	loglog(zippingTimes, zipCte * zippingTimes.^zipExponent, '-k');

%	figure
%	hold on;
%	for i = 1:N
%		boundsAtZippingTime = bounds{i}(:, zippingIndices(i))
%		[f, b] = normalizedHist(boundsAtZippingTime, 30);
%		plot(b, f/max(f) + i);
%	end
%	hold off
%	sleep(1e-9);
end

zipC = zipCte;
zipE = 1/zipExponent;

hold off;

