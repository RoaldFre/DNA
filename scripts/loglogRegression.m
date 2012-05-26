%xs must be a column vector
function [cte, exponent, cteStddev, exponentStddev] = loglogRegression(xs, ys, guessCte, guessExponent, yerr)

if (nargin < 4)
	error("Not enough required arguments!");
end
if (nargin == 4)
	yerr = ones(size(ys));
end


logys = log(ys);
logxs = log(xs);
logysErr = yerr ./ ys;

%figure
%errorbar(logxs/log(10), logys/log(10), logysErr);

weights = logysErr.^(-1);

F = [ones(numel(logxs), 1), logxs];

[p, yVar, r, pVar] = LinearRegression(F, logys, weights);
cteExponent = p(1);
cteExponentStddev = sqrt(pVar(1));
exponent = p(2);
exponentStddev = sqrt(pVar(2));

cte = exp(cteExponent);
cteStddev = cte * cteExponentStddev;

