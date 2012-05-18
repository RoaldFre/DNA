% Fit the function ys = cte * xs.^exponent
function [cte, exponent, cteStddev, exponentStddev] = loglogregression(xs, ys, guessCte, guessExponent, yerr)

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
		xs, ys, [guessCte, guessExponent],
		@(x,p)(p(1) * x.^p(2)),
		acceptedError,
		maxIterations,
		weights);

cte = p(1);
exponent = p(2);
disp("This should be +/- diagonal:")
disp(covp);
errs = sqrt(diag(covp));
cteStddev = errs(1);
exponentStddev = errs(2);
