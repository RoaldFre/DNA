% Fit the function ys = offset + amplitude * exp(-xs/tau)
% function [tau, offset, amplitude, tauStddev, offsetStddev, amplitudeStddev] = endToEndRegression(xs, ys, guessTau, guessOffset, guessAmplitude, yerr)
function [tau, offset, amplitude, tauStddev, offsetStddev, amplitudeStddev] = endToEndRegression(xs, ys, guessTau, guessOffset, guessAmplitude, yerr)

addpath('octave-forge');

if (nargin < 4)
       error("not enough required arguments!");
end
if (nargin < 5)
       yerr = ones(length(ys), 1);
end


%Attempt at a better conditioning of the problem:
%xfact = mean(xs);
%yfact = mean(xs);
%xs = xs / xfact;
%ys = ys / yfact;
%guessTau = guessTau / xfact;
%guessOffset = guessOffset / yfact;
%guessAmplitude = guessAmplitude / yfact;


weights = yerr.^(-1);

acceptedError = 1e-10;
maxIterations = 1000;

[fr, p, kvg, iter, corp, covp, covr, stdresid, Z, r2] = leasqr(
               xs, ys, [guessTau, guessOffset, guessAmplitude],
               @(x,p)(p(2) + p(3)*exp(-x/p(1))),
               acceptedError,
               maxIterations,
               weights);

tau = p(1);
offset = p(2);
amplitude = p(3);
disp("This should be +/- diagonal:")
disp(covp);
errs = sqrt(diag(covp));
tauStddev = errs(1);
offsetStddev = errs(2);
amplitudeStddev = errs(3);


%tau = tau * xfact;
%offset = offset * yfact;
%amplitude = amplitude * yfact;
%tauStddev = tauStddev * xfact;
%offsetStddev = offsetStddev * yfact;
%amplitudeStddev = amplitudeStddev * yfact;
