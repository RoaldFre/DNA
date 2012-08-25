% Fit the function ys = 1 / expcte * xs.^exponent
%function [A, B, C, Astddev, Bstddev, Cstddev] = meltingTempRegression(xs, ys, guessA, guessB, guessC, yerr)
function [A, B, Astddev, Bstddev] = meltingTempRegression(xs, ys, guessA, guessB, yerr)

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
               xs, ys, [guessA, guessB],
               @(x,p)(0.5 - 0.5*tanh(p(1).*(x - p(2)))),
               acceptedError,
               maxIterations,
               weights);

A = p(1);
B = p(2);
%C = p(3);
disp("This should be +/- diagonal:")
disp(covp);
errs = sqrt(diag(covp));
Astddev = errs(1);
Bstddev = errs(2);
%Cstddev = errs(3);
