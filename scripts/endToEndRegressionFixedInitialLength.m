% Fit the function ys = offset + (initialLength - offset) * exp(-xs/tau)
% for a fixed value of initialLength.
% function [tau, offset, tauStddev, offsetStddev] = endToEndRegressionFixedInitialLength(xs, ys, initialLength, guessTau, guessOffset, yerr)
function [tau, offset, tauStddev, offsetStddev] = endToEndRegressionFixedInitialLength(xs, ys, initialLength, guessTau, guessOffset, yerr)

addpath('octave-forge');

if (nargin < 4)
       error("not enough required arguments!");
end
if (nargin < 5)
       yerr = ones(length(ys), 1);
end


weights = yerr.^(-1);

acceptedError = 1e-10;
maxIterations = 1000;

[fr, p, kvg, iter, corp, covp, covr, stdresid, Z, r2] = leasqr(
               xs, ys, [guessTau, guessOffset],
               @(x,p)(p(2) + (initialLength - p(2))*exp(-x/p(1))),
               acceptedError,
               maxIterations,
               weights);

tau = p(1);
offset = p(2);
disp("This should be +/- diagonal:")
disp(covp);
errs = sqrt(diag(covp));
tauStddev = errs(1);
offsetStddev = errs(2);

