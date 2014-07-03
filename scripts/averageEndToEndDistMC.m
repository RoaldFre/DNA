function [tau, tauStddev, offset, offsetStddev, amplitude, amplitudeStddev, sweeps, meanEndToEndDist, errEndToEndDist] = averageEndToEndDistMC(filesglob, decimateFactor)

addpath("generic");

timePerSweep = 1e-12; % a single sweeep is a pico second

if (nargin < 1)
	error "Not enough required arguments given!"
elseif (nargin < 2)
	decimateFactor = 1;
end

[time, dists] = parseEndToEndDists(filesglob, decimateFactor);
sweeps = time / timePerSweep;

% XXX
% This is a hack to make sure that all times are positive (for stretched 
% exponential fit). Due to the decimation, it is possible that we end up 
% with a tiny negative first time value.
sweeps = sweeps - sweeps(1);


nRuns = numel(dists(1,:))
meanEndToEndDist = mean(dists')';
errEndToEndDist = std(dists')' / sqrt(nRuns - 1);

color = "b";
clf; hold on;
h1 = errorbar(sweeps, meanEndToEndDist, errEndToEndDist);
set(h1, "color", color);
set(h1, "marker", ".");
set(h1, "linestyle", "none");
set(h1, "linewidth", 1);


initialLength = meanEndToEndDist(1); % This should have no error, as all strands have the same initial length!

[tau, offset, amplitude, tauStddev, offsetStddev, amplitudeStddev] = exponentialRelaxation(sweeps, meanEndToEndDist, 10, initialLength / 2, initialLength / 2, errEndToEndDist);
tauInfo = [tau, tauStddev]
offsetInfo = [offset, offsetStddev]
amplitudeInfo = [amplitude, amplitudeStddev]

#fit = offset + (initialLength - offset) * exp(-(sweeps/tau).^beta);
fit = offset + amplitude * exp(-sweeps/tau);

plot(sweeps, fit, "k", "linewidth", 4);
pause(1e-9);
hold off;
