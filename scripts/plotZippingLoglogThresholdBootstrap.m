function [exponent, expErr, zippingTimes, zippingBounds] = plotZippingLoglogThresholdBootstrap(Ns, times, bounds, nuclBoundBPs, zippedBoundBPs, color, regressionBeginIndex)

addpath generic

if nargin < 7
	regressionBeginIndex = 1;
end

bootstrapSamples = 30;

% Plot the mean curves
[originalZippingTimes, originalZippingBounds] = plotZippingLoglogThreshold(times, bounds, nuclBoundBPs, zippedBoundBPs, color, regressionBeginIndex)

% Do the boostrapping
N = numel(times);
resampledMeanBound = cell(size(bounds));
zippingTimes = zeros(N, bootstrapSamples);
zippingBounds = zeros(N, bootstrapSamples);
exponents = zeros(1, bootstrapSamples);
constants = zeros(1, bootstrapSamples);
for i = 1:bootstrapSamples
	printf("resampling %d of %d\n", i, bootstrapSamples);
	for r = 1:N
		resampledMeanBound{r} = bootstrapSampleMean(bounds{r});
	end
	[theseZippingTimes, theseZippingBounds, zipE, zipC] = plotZippingLoglogThreshold(times, resampledMeanBound, nuclBoundBPs, zippedBoundBPs, [], regressionBeginIndex);
	zippingTimes(:, i) = theseZippingTimes;
	zippingBounds(:, i) = theseZippingBounds;

	[cte, exponent] = loglogRegression(Ns(regressionBeginIndex:end), theseZippingTimes(regressionBeginIndex:end), 1, 1)
	exponents(i) = exponent;
	constants(i) = cte;
end

zippingTimeErrs = std(zippingTimes')';
zippingTimeMean = mean(zippingTimes')';


[Ns', originalZippingTimes, zippingTimeMean, zippingTimeErrs]

[cte, exponent] = loglogRegression(Ns(regressionBeginIndex:end), originalZippingTimes(regressionBeginIndex:end), 1, 1);

cteErr = std(constants);
expErr = std(exponents);

% standard deviation of standard deviation
% http://stats.stackexchange.com/questions/631/standard-deviation-of-standard-deviation
n = bootstrapSamples;
stdOfStd = expErr * gamma((n-1)/2)/gamma(n/2) * sqrt((n-1)/2 - (gamma(n/2)/gamma((n-1)/2))^2);
[exponent, expErr, stdOfStd]


figure
hold on
loglogerror(Ns, originalZippingTimes, zippingTimeErrs);
loglog(Ns, cte * Ns.^exponent, 'k');
hold off


% Histograms of bootstrapped times and exponents to check for normality:
% Pick bootstrapSamples sufficiently large (eg 500)!
%figure
%hold on;
%for i = 1:N
%	[f, b] = normalizedHist(zippingTimes(i, :), 20);
%	plot(b, f/max(f) + i);
%end
%title("bootstrapped zippingTimes")
%hold off
%sleep(1e-9);
%
%figure
%normalizedHist(exponents, 30)
%title("exponents")
