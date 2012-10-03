% Fit the sum of two exponential relaxation functions
%   ys = offset +    amplFrac    * (initialLength - offset) * exp(-xs/tau1)
%               + (1 - amplFrac) * (initialLength - offset) * exp(-xs/tau2)
% for a fixed value of initialLength, and where tau1 >= tau2.
% function [tau, offset, beta, tauStddev, offsetStddev, betaStddev] = endToEndRegressionFixedInitialLength(xs, ys, initialLength, guessTau, guessOffset, guessBeta, yerr)
function [tau1, tau2, offset, amplFrac, tau1Stddev, tau2Stddev, offsetStddev, amplFracStddev] = endToEndRegressionFixedInitialLength(xs, ys, initialLength, guessTau1, guessTau2, guessOffset, guessAmplFrac, yerr)

addpath('octave-forge');

if (nargin < 7)
       error("not enough required arguments!");
end
if (nargin < 8)
       yerr = ones(length(ys), 1);
end

% To improve conditioning of the problem: express time in 10^-6s and 
% lengths in 10^-8 m
timeUnit = 1e-6;
lengthUnit = 1e-8;
xs = xs / timeUnit;
ys = ys / lengthUnit;
guessTau1 = guessTau1 / timeUnit;
guessTau2 = guessTau2 / timeUnit;
guessOffset = guessOffset / lengthUnit;
initialLength = initialLength / lengthUnit;


weights = yerr.^(-1);

acceptedError = 1e-10;
maxIterations = 1000;


[fr, p, kvg, iter, corp, covp, covr, stdresid, Z, r2] = leasqr(
               xs, ys, [guessTau1, guessTau2, guessOffset, guessAmplFrac],
               @(x,p)(p(3) +   p(4)  *(initialLength - p(3))*exp(-x/p(1)) ...
                           + (1-p(4))*(initialLength - p(3))*exp(-x/p(2))),
               acceptedError,
               maxIterations,
               weights);

tau1 = p(1);
tau2 = p(2);
offset = p(3);
amplFrac = p(4);
disp("This should be +/- diagonal:")
disp(covp);
errs = sqrt(diag(covp));
tau1Stddev = errs(1);
tau2Stddev = errs(2);
offsetStddev = errs(3);
amplFracStddev = errs(4);

% Back to proper units
tau1 = tau1 * timeUnit;
tau2 = tau2 * timeUnit;
tau1Stddev = tau1Stddev * timeUnit
tau2Stddev = tau2Stddev * timeUnit
offset = offset * lengthUnit;
offsetStddev = offsetStddev * lengthUnit;

if tau1 < tau2
	temp = tau1;
	tau1 = tau2;
	tau2 = temp;
	temp = tau1Stddev;
	tau1Stddev = tau2Stddev;
	tau2Stddev = temp;
end
