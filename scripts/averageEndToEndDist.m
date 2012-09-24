function [tau, tauStddev, offset, offsetStddev, amplitude, amplitudeStddev, time, meanEndToEndDist, errEndToEndDist] = averageEndToEndDist(filesglob, decimateFactor)

if (nargin < 1)
	error "Not enough required arguments given!"
elseif (nargin < 2)
	decimateFactor = 1;
end

[time, dists] = parseEndToEndDists(filesglob, decimateFactor);

nRuns = numel(dists(1,:))
meanEndToEndDist = mean(dists')';
errEndToEndDist = std(dists')' / sqrt(nRuns - 1);

initialLength = meanEndToEndDist(1); % This should have no error, as all strands have the same initial length!

[tau, offset, tauStddev, offsetStddev] = endToEndRegressionFixedInitialLength(time, meanEndToEndDist, initialLength, 1e-7, 1e-8, errEndToEndDist);
tauInfo = [tau, tauStddev]
offsetInfo = [offset, offsetStddev]

fit = offset + (initialLength - offset) * exp(-time/tau);

color = "b";
clf; hold on;
h1 = errorbar(time, meanEndToEndDist, errEndToEndDist);
set(h1, "color", color);
set(h1, "marker", ".");
set(h1, "linestyle", "none");
set(h1, "linewidth", 1);

plot(time, fit, "k", "linewidth", 4);
pause(1e-9);
hold off;
