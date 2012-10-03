% Fit the stretched exponential function
%   ys = offset + (initialLength - offset) * exp[-(xs/tau)^beta]
% for a fixed value of initialLength.
% function [tau, offset, beta, tauStddev, offsetStddev, betaStddev] = endToEndRegressionFixedInitialLength(xs, ys, initialLength, guessTau, guessOffset, guessBeta, yerr)
function [tau, offset, beta, tauStddev, offsetStddev, betaStddev] = endToEndRegressionFixedInitialLength(xs, ys, initialLength, guessTau, guessOffset, guessBeta, yerr)

addpath('octave-forge');

if (nargin < 6)
       error("not enough required arguments!");
end
if (nargin < 7)
       yerr = ones(length(ys), 1);
end

% To improve conditioning of the problem: express time in 10^-6s and 
% lengths in 10^-8 m
timeUnit = 1e-6;
lengthUnit = 1e-8;
xs = xs / timeUnit;
ys = ys / lengthUnit;
guessTau = guessTau / timeUnit;
guessOffset = guessOffset / lengthUnit;
initialLength = initialLength / lengthUnit;


weights = yerr.^(-1);

acceptedError = 1e-10;
maxIterations = 1000;


[fr, p, kvg, iter, corp, covp, covr, stdresid, Z, r2] = leasqr(
               xs, ys, [guessTau, guessOffset, guessBeta],
               @(x,p)(p(2) + (initialLength - p(2))*exp(-(x/p(1)).^p(3))),
               acceptedError,
               maxIterations,
               weights);

tau = p(1);
offset = p(2);
beta = p(3);
disp("This should be +/- diagonal:")
disp(covp);
errs = sqrt(diag(covp));
tauStddev = errs(1);
offsetStddev = errs(2);
betaStddev = errs(3);

% Back to proper units
tau = tau * timeUnit;
tauStddev = tauStddev * timeUnit
offset = offset * lengthUnit;
offsetStddev = offsetStddev * lengthUnit;
