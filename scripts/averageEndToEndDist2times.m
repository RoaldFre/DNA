% function [tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, amplFrac, amplFracStddev, time, meanEndToEndDist, errEndToEndDist] = averageEndToEndDist2times(filesglob, decimateFactor)
function [tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, amplFrac, amplFracStddev, time, meanEndToEndDist, errEndToEndDist] = averageEndToEndDist2times(filesglob, decimateFactor)

if (nargin < 1)
	error "Not enough required arguments given!"
elseif (nargin < 2)
	decimateFactor = 1;
end

[time, dists, N] = parseEndToEndDists(filesglob, decimateFactor);

% XXX
% This is a hack to make sure that all times are positive (for stretched 
% exponential fit). Due to the decimation, it is possible that we end up 
% with a tiny negative first time value.
time = time - time(1);


nRuns = numel(dists(1,:))
meanEndToEndDist = mean(dists')';
errEndToEndDist = std(dists')' / sqrt(nRuns - 1);

color = "b";
clf; hold on;
h1 = errorbar(time, meanEndToEndDist, errEndToEndDist);
set(h1, "color", color);
set(h1, "marker", ".");
set(h1, "linestyle", "none");
set(h1, "linewidth", 1);


initialLength = meanEndToEndDist(1); % This should have no error, as all strands have the same initial length!

tauGuess = 9.7e-8 * (N / 84.0) ^ 2.088;
tau1Guess = tauGuess;
tau2Guess = tauGuess / 5;

[tau1, tau2, offset, amplFrac, tau1Stddev, tau2Stddev, offsetStddev, amplFracStddev] = endToEndRegressionFixedInitialLength2times(time, meanEndToEndDist, initialLength, tau1Guess, tau2Guess, 1.3e-8, 0.8, errEndToEndDist);
tau1Info = [tau1, tau1Stddev]
tau2Info = [tau2, tau2Stddev]
offsetInfo = [offset, offsetStddev]
amplFracInfo = [amplFrac, amplFracStddev]

fit = offset +   amplFrac  *(initialLength - offset) * exp(-time/tau1) ...
             + (1-amplFrac)*(initialLength - offset) * exp(-time/tau2);

plot(time, fit, "g", "linewidth", 4);
pause(1e-9);
hold off;
