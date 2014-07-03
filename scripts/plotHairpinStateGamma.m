% function plotHairpinStateGamma(N, gamma)
function plotHairpinStateGamma(N, gamma)

[bound, time] = parseHairpinStateToArray(['/home/other/roald/clusterdata/hairpinState_gamma/relaxT90C_sampleT10C_time10000_sweeps100000_NXY4/dt15_gamma',gamma,'/N',num2str(N),'/*.itf11']);


time = time - time(1);

bound = bound - 4; % correct for XY pairs

% remove runs where XY hasn't zipped yet
goodRuns = find(bound(:,1) >= 0);
nRunsAll = numel(bound(:,1))
nRuns = numel(goodRuns)
bound = bound(goodRuns, :);


meanBound = mean(bound);
errBound = std(bound) / sqrt(nRuns - 1);

errFrac = 1 + errBound ./ meanBound;

%clf;
hold on;
loglog(time, meanBound, 'b');
loglog(time, meanBound .* errFrac, 'r');
loglog(time, meanBound ./ errFrac, 'r');

% last part
%is=find(t1<time & time<t2);
%xs=time(is);
%ys=meanBound(is);
%[cte, exponent, cteStddev, exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6)
%loglog(xs, cte*xs.^exponent, 'k', 'linewidth', 3)

