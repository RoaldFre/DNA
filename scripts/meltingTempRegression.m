% Fit the function ys = (0.5 - 0.5*tanh(A.*(xs - B)))
% function [A, B, Astddev, Bstddev] = meltingTempRegression(xs, ys, guessA, guessB, yerr)
function [A, B, Astddev, Bstddev] = meltingTempRegression(xs, ys, guessA, guessB, yerr)

addpath('octave-forge');

if (nargin < 5)
       error("not enough required arguments!");
end

[fr, p, pErr] = leasqrError(
               xs, ys, yerr, [guessA, guessB],
               @(x,p)(0.5 - 0.5*tanh(p(1).*(x - p(2)))));

A = p(1);
B = p(2);
Astddev = pErr(1);
Bstddev = pErr(2);
