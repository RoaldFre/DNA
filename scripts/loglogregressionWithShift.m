% Fit the function ys = cte * xs.^exponent
function [cte, exponent, shift, cteStddev, exponentStddev, shiftStddev] = loglogregressionWithShift(xs, ys, guessCte, guessExponent, guessShift, yerr)

if (nargin < 5)
	error("not enough required arguments!");
end
if (nargin < 6)
	yerr = ones(length(ys), 1);
end

weights = yerr.^(-1);

acceptedError = 1e-10;
maxIterations = 1000;


[fr, p, kvg, iter, corp, covp, covr, stdresid, Z, r2] = leasqr(
		xs, ys, [guessCte, guessExponent, guessShift],
		@(x,p)(p(1) * max(x - p(3), 1e-100).^p(2)), %hack to avoid imaginary result
		acceptedError,
		maxIterations,
		weights);

cte = p(1);
exponent = p(2);
shift = p(3);
disp("This should be +/- diagonal:")
disp(covp);
errs = sqrt(diag(covp));
cteStddev = errs(1);
exponentStddev = errs(2);
shiftStddev = errs(3);
end


