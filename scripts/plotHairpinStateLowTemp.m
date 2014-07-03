%function plotHairpinStateLowTemp(N, T, sweeps, loadFromCache)
function [time, meanBound, errBound, bound] = plotHairpinStateLowTemp(N, T, totTime, sweeps, wait, loadFromCache, plot)

addpath("~/DNA/scripts/generic");

if nargin < 6
	loadFromCache = false;
end
if nargin < 7
	plot = true;
end

if ischar(wait)
	waitStr = wait;
else
	waitStr = num2str(wait);
end

init_t1 = -1; % guard
end_t1 = -1; % guard

if T == 10 || T == 30
	if N == 7
		%totTime = 100;
		t1=4.0e-9;
		t2=1.8e-8;
	elseif N == 10
		%totTime = 100;
		init_t1 = 2e-10;
		init_t2 = 3e-9;
		t1=4.0e-9;
		t2=1.8e-8;
	elseif N == 14
		%totTime = 150;
		init_t1 = 2.5e-10;
		init_t2 = 3e-9;
		t1=4.0e-9;
		t2=2.2e-8;
	elseif N == 20
		%totTime = 200;
		init_t1 = 3.6e-10;
		init_t2 = 1.1e-8;
		t1 = 2.0e-8;
		t2 = 7.0e-8;
	elseif N == 28
		%totTime = 300;
		init_t1 = 3.6e-10;
		init_t2 = 2.0e-8;
		t1 = 3.3e-8;
		t2 = 1.3e-7;
	elseif N == 40
		%totTime = 500;
		init_t1 = 4.0e-10;
		init_t2 = 2.4e-8;
		t1 = 3.6e-8;
		t2 = 1.45e-7;
		end_t1 = 1.6e-7;
		end_t2 = 2.0e-7;
	elseif N == 57
		%totTime = 800;
		init_t1 = 4.0e-10;
		init_t2 = 3.3e-8;
		t1 = 5.0e-8;
		t2 = 3.4e-7;
		end_t1 = 3.6e-7;
		end_t2 = 5.0e-7;
	else
		error blah
	end
end


%sweeps=100000;
%sweeps=10000;
dir = ["/home/other/roald/clusterdata_properVelocityInit/hairpinState_lowTemp/T",num2str(T),"C_time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_wait",waitStr,"_NXY4/dt15/N",num2str(N)]
%dir = ["/home/other/roald/clusterdata/hairpinState_lowTemp_gamma/T",num2str(T),"C_time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_NXY4/dt15/g5e11/N",num2str(N)]
%dir = ["/home/other/roald/clusterdata/hairpinState_lowTemp/T",num2str(T),"C_time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_relaxFactor100_NXY4/dt15/N",num2str(N)]
if not(loadFromCache)
	[bound, time] = parseHairpinStateToArray([dir, '/*.itf11']);
	%save('-binary', '-z', [dir,"/hairpinStateN",num2str(N),"_onlyBound"], 'bound', 'time')
else
	error yada yada
	load([dir,"/hairpinStateN",num2str(N),"_onlyBound"]);
end


time = time - time(1);

bound = bound - 4; % correct for XY pairs

% remove runs where XY hasn't zipped yet
goodRuns = find(bound(:,1) >= 0);
nRunsAll = numel(bound(:,1))
nRuns = numel(goodRuns)
bound = bound(goodRuns, :);




% XXX
initialIndex = 3;
time = time(initialIndex : end);
bound = bound(:, initialIndex:end);




meanBound = mean(bound);
errBound = std(bound) / sqrt(nRuns - 1);

errFrac = 1 + errBound ./ meanBound;


if not(plot)
	return
end


%clf;
hold on;
loglog(time, meanBound, 'b');
loglog(time, meanBound .* errFrac, 'r');
loglog(time, meanBound ./ errFrac, 'r');

% last part
is=find(t1<time & time<t2);
xs=time(is);
ys=meanBound(is);
yerrs=errBound(is)';
[cte, exponent, cteStddev, exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
loglog(xs, cte*xs.^exponent, 'k', 'linewidth', 3)


% initial part
if init_t1 > 0
	is=find(init_t1<time & time<init_t2);
	xs=time(is);
	ys=meanBound(is);
	yerrs=errBound(is)';
	[init_cte, init_exponent, init_cteStddev, init_exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
	loglog(xs, init_cte*xs.^init_exponent, 'g', 'linewidth', 3)
end

% final part
if end_t1 > 0
	is=find(end_t1<time & time<end_t2);
	xs=time(is);
	ys=meanBound(is);
	yerrs=errBound(is)';
	[end_cte, end_exponent, end_cteStddev, end_exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
	loglog(xs, end_cte*xs.^end_exponent, 'm', 'linewidth', 3)
end

xlabel("time (s)")
ylabel("average number of bound base pairs");

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
if end_t1 > 0
	[1/end_exponent, end_exponentStddev/end_exponent^2]
end










time = time + time(2) - time(1); %don't start at t=0 for loglog regression

%[exponent, x1, xOffset, yOffset, cutOff, cutOffWidth, exponentErr, x1err, xOffsetErr, yOffsetErr, cutOffErr, cutOffWidthErr, f] = plOffsetCutoffLoglogFit(time, meanBound, errBound, 0.8, 1e-7);
%[1/exponent, x1, xOffset, yOffset, cutOff, cutOffWidth; exponentErr, x1err, xOffsetErr/exponent^2, yOffsetErr, cutOffErr, cutOffWidthErr]

%[exp1, x1, exp2, x2, xOffset, yOffset, cutOff, cutOffWidth, exp1err, x1err, exp2err, x2err, xOffsetErr, yOffsetErr, cutOffErr, cutOffWidthErr, f] = pl2OffsetCutoffFit(time, meanBound, errBound, 0.8, 1e-9, 0.6, 1e-9);
%[1/exp1, x1, 1/exp2, x2, xOffset, yOffset, cutOff, cutOffWidth; exp1err, x1err, exp2err, x2err, xOffsetErr, yOffsetErr, cutOffErr, cutOffWidthErr]

%[exp1, exp2, x1, xCrossOver, cutOff, cutOffWidth, exp1err, exp2err, x1err, xCrossOverErr, cutOffErr, cutOffWidthErr, f] = pl2CutoffFit(time', meanBound, errBound, 1, 0.5, 2e-9, 5e-8, 1e-7, 1e-8);
%[exp1, exp2, x1, xCrossOver, cutOff, cutOffWidth; exp1err, exp2err, x1err, xCrossOverErr, cutOffErr, cutOffWidthErr]

%[exp1, exp2, x1, xCrossOver, cutOff, cutOffWidth, xOffset, yOffset, exp1err, exp2err, x1err, xCrossOverErr, cutOffErr, cutOffWidthErr, xOffsetErr, yOffsetErr, f] = pl2OffsetCutoffFit(time', meanBound, errBound, 1, 2.5, 3e-9, 1e-8, 1e-7, 1e-8);
[exp1, exp2, x1, xCrossOver, cutOff, cutOffWidth, xOffset, yOffset, exp1err, exp2err, x1err, xCrossOverErr, cutOffErr, cutOffWidthErr, xOffsetErr, yOffsetErr, f] = pl2OffsetCutoffLoglogFit(time', meanBound, errBound, 1, 1/2.5, -8.5, -8);
[1/exp1, 1/exp2, x1, xCrossOver, cutOff, cutOffWidth, xOffset, yOffset; exp1err, exp2err, x1err, xCrossOverErr, cutOffErr, cutOffWidthErr, xOffsetErr, yOffsetErr]


numSamples = numel(f)
isReal = isreal(f)
numNegative = sum(f <= 0)

hold on;
loglog(time, f, 'k', "linewidth", 2);
loglog(time + xOffset, meanBound - yOffset, 'b');
loglog(time + xOffset, meanBound .* errFrac - yOffset, 'r');
loglog(time + xOffset, meanBound ./ errFrac - yOffset, 'r');
loglog(time + xOffset, f - yOffset, 'g', 'linewidth', 2);
hold off;




%saveData = [time - time(1), meanBound', errBound'];
%save('-ascii', ['./data/meanBoundN',num2str(N)], 'saveData');
