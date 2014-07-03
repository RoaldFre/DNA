function [tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, ampl1, ampl1Stddev, ampl2, ampl2Stddev, sweeps, meanEndToEndDist, errEndToEndDist] = averageEndToEndDist2timesMC(filesglob, decimateFactor)
% FIXED INITIAL LENGTH:
%function [tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, amplFrac, amplFracStddev, sweeps, meanEndToEndDist, errEndToEndDist] = averageEndToEndDist2timesMC(filesglob, decimateFactor)

addpath("generic");

timePerSweep = 1e-12; % a single sweeep is a pico second

if (nargin < 1)
	error "Not enough required arguments given!"
elseif (nargin < 2)
	decimateFactor = 1;
end

[time, dists, N] = parseEndToEndDists(filesglob, decimateFactor);
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

tauGuess = 50;
tau1Guess = tauGuess;
tau2Guess = tauGuess / 5;

% FIXED INITIAL LENGTH
%[tau1, tau2, offset, amplFrac, tau1Stddev, tau2Stddev, offsetStddev, amplFracStddev] = exponentialRelaxationFixedInitial2times(sweeps, meanEndToEndDist, initialLength, tau1Guess, tau2Guess, 1.3e-8, 0.8, errEndToEndDist);
%tau1Info = [tau1, tau1Stddev]
%tau2Info = [tau2, tau2Stddev]
%offsetInfo = [offset, offsetStddev]
%amplFracInfo = [amplFrac, amplFracStddev]
%
%fit = offset +   amplFrac  *(initialLength - offset) * exp(-sweeps/tau1) ...
%             + (1-amplFrac)*(initialLength - offset) * exp(-sweeps/tau2);


[tau1, tau2, offset, ampl1, ampl2, tau1Stddev, tau2Stddev, offsetStddev, ampl1Stddev, ampl2Stddev] = exponentialRelaxation2times(sweeps, meanEndToEndDist, tau1Guess, tau2Guess, initialLength / 3, initialLength / 3, initialLength / 3, errEndToEndDist);
tau1Info = [tau1, tau1Stddev]
tau2Info = [tau2, tau2Stddev]
offsetInfo = [offset, offsetStddev]
ampl1info = [ampl1, ampl1Stddev]
ampl2info = [ampl2, ampl2Stddev]

fit = offset + ampl1 * exp(-sweeps/tau1) ...
             + ampl2 * exp(-sweeps/tau2);



plot(sweeps, fit, "g", "linewidth", 4);
pause(1e-9);
hold off;
