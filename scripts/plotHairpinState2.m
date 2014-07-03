function plotHairpinState2(N, loadFromCache, timeShiftFactor)

if nargin < 2
	loadFromCache = false;
end
if nargin < 3
	timeShiftFactor = 0;
end

%timeShiftFactor = 0;

init_t1 = -1; % guard

if N == 7
	totTime = 300;
	timeShift = 0e-9;
	t1=4.0e-9;
	t2=1.8e-8;
elseif N == 10
	totTime = 80;
	timeShift = 2e-9;
	t1=4.0e-9;
	t2=1.8e-8;
elseif N == 14
	totTime = 300;
	timeShift = 0e-9;
	t1=3.0e-9;
	t2=1.8e-8;
elseif N == 20
	totTime = 200;
	timeShift = 1e-8;
	init_t1 = 1.4e-8;
	init_t2 = 2.7e-8;
	t1 = 3.5e-8;
	t2 = 6.9e-8;
elseif N == 28
	totTime = 500;
	timeShift = 0e-8;
	t1=5.5e-8;
	t2=1.3e-7;
elseif N == 40
	totTime = 500;
	timeShift = 1e-8;
	init_t1 = 1.5e-8;
	init_t2 = 1.1e-7;
	t1 = 1.2e-7;
	t2 = 2.4e-7;
elseif N == 57
	totTime = 800;
	timeShift = 1e-8;
	init_t1 = 0.2e-7;
	init_t2 = 1.4e-7;
	t1 = 2.5e-7;
	t2 = 4.4e-7;
else
	error blah
end



if not(loadFromCache)
	glob = ['/home/other/roald/clusterdata/hairpinState2/relaxT90C_sampleT10C_time',num2str(totTime),'_sweeps10000_NXY4/dt15/N',num2str(N),'/*itf11'];
	[bound, time] = parseHairpinStateToArray(glob);
	save('-binary', '-z', ["./data/hairpinState2N",num2str(N),"_onlyBound"], 'bound', 'time')
else
	load(["./data/hairpinState2N",num2str(N),"_onlyBound"]);
end


time = time - time(1) + timeShift * timeShiftFactor;

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
is=find(t1<time & time<t2);
xs=time(is);
ys=meanBound(is);
[cte, exponent, cteStddev, exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6)
loglog(xs, cte*xs.^exponent, 'k', 'linewidth', 3)


% initial part
if init_t1 > 0
	is=find(init_t1<time & time<init_t2);
	xs=time(is);
	ys=meanBound(is);
	[init_cte, init_exponent, init_cteStddev, init_exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6)
	loglog(xs, init_cte*xs.^init_exponent, 'g', 'linewidth', 3)
end


hold off;


%figure
%clf; hold on;
%plot(time, meanBound, 'b');
%plot(time, meanBound + errBound, 'r');
%plot(time, meanBound - errBound, 'r');
%hold off;

if init_t1 > 0
	[1/init_exponent, init_exponentStddev/init_exponent^2]
end
[1/exponent, exponentStddev/exponent^2]


saveData = [time - time(1), meanBound', errBound'];
save('-ascii', ['./data/meanBound2N',num2str(N)], 'saveData');
