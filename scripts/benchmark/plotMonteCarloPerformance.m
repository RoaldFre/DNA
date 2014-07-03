data = load('monteCarloPerformance');

addpath ../generic/

Ns = data(:,1);
sweeps = data(:,2);
times = data(:,3:end);

nNs = size(times, 1);
nRuns = size(times, 2);

meanTimes = mean(times')';
errTimes = std(times')' / sqrt(nRuns - 1);

loglogerror(Ns, meanTimes ./ sweeps, errTimes ./ sweeps);

allTimesPerSweep = zeros(1, nRuns * nNs);
allNs = zeros(1, nRuns * nNs);
for n = 1:nNs
	allNs(1 + (n-1)*nRuns : n*nRuns) = Ns(n);
	allTimesPerSweep(1 + (n-1)*nRuns : n*nRuns) = times(n, :) / sweeps(n);
end

[cte, exponent, cteStddev, exponentStddev] = loglogRegressionBootstrap(allNs, allTimesPerSweep, 1, 1);

[exponent, exponentStddev]
[cte, cteStddev]

hold on;
loglog(Ns, cte*Ns.^exponent);
hold off;

