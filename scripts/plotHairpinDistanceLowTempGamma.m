% function plotHairpinStateLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, loadFromCache)
function plotHairpinStateLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, loadFromCache)

addpath("generic");

doLinRegs = true;

if nargin < 8
	loadFromCache = false;
end

if ischar(wait)
	waitStr = wait;
else
	waitStr = num2str(wait);
end

baseLineGamma = 5e12;
gamma = str2num(gammaStr);

init_t1 = -1; % guard
init_t2 = -1; % guard
t1 = -1; %guard
t2 = -1; %guard
end_t1 = -1; % guard
end_t2 = -1; % guard

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
	end
end


%basedir = '/home/other/roald/clusterdata/hairpinState_lowTemp_gamma';
%basedir = '/home/other/roald/clusterdata_properVelocityInit_wrongResizing/hairpinState_lowTemp_gamma';
basedir = '/home/other/roald/clusterdata_properVelocityInit/hairpinState_lowTemp_gamma';
dir = [basedir,"/T",num2str(T),"C/g",gammaStr,"/N",num2str(N),"/dt",num2str(dt),"/time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_wait",waitStr,"_NXY4"]

filesglob = [dir, '/*.itf11'];
if not(loadFromCache)
	[bound, time] = parseHairpinStateToArray(filesglob);
	bound = bound - 4; % correct for XY pairs

	%below: problem with '*' in wait!
	%save('-binary', '-z', [dir,"/hairpinStateN",num2str(N),"_onlyBound"], 'bound', 'time')
	save('-binary', '-z', ['./data/hairpinStateN',num2str(N)], 'bound', 'time');
else
	%error blargh blargh blargh
	%load([dir,"/hairpinStateN",num2str(N),"_onlyBound"]);
	load(['./data/hairpinStateN',num2str(N)]);
end


time = time - time(1);

% Rescale time to 'base line'
%time = time * baseLineGamma/gamma;
%time = time * 9; # for g = 5e11
%time = time * 18; # for g = 1.6e11
init_t1 /= baseLineGamma/gamma;
init_t2 /= baseLineGamma/gamma;
t1 /= baseLineGamma/gamma;
t2 /= baseLineGamma/gamma;
end_t1 /= baseLineGamma/gamma;
end_t2 /= baseLineGamma/gamma;


% remove runs where XY hasn't zipped yet
%goodRuns = find(bound(:,1) >= 0);
%nRunsAll = numel(bound(:,1))
%nRuns = numel(goodRuns)
%bound = bound(goodRuns, :);

[runsWithUnboundXY, _] = find(bound < 0);
sort(runsWithUnboundXY);
runsWithUnboundXY = unique(runsWithUnboundXY);
% remove runs where XY unzipped at least once at any moment in time
goodRuns = setdiff(1:numel(bound(:,1)), runsWithUnboundXY);

%files = glob(filesglob);
%for i = 1:numel(runsWithUnboundXY)
%	numFuckedUp = sum(bound(runsWithUnboundXY(i),:) < 0);
%	printf("%d %s\n", numFuckedUp, files{runsWithUnboundXY(i)});
%end
%
%clf; hold on
%for i = runsWithUnboundXY
%	plot(time, bound(i, :));
%end

nRuns = numel(goodRuns)
bound = bound(goodRuns, :);





meanBound = mean(bound);
errBound = std(bound) / sqrt(nRuns - 1);

errFrac = 1 + errBound ./ meanBound;

%clf;
hold on;
%plot(time, meanBound, 'b');
%plot(time, meanBound + errBound, 'r');
%plot(time, meanBound - errBound, 'r');
%return
loglog(time, meanBound, 'b');
loglog(time, meanBound .* errFrac, 'r');
loglog(time, meanBound ./ errFrac, 'r');


if doLinRegs

% initial part
if init_t1 > 0
	is=find(init_t1<time & time<init_t2);
	xs=time(is);
	ys=meanBound(is);
	yerrs=errBound(is)';
	[init_cte, init_exponent, init_cteStddev, init_exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
	loglog(xs, init_cte*xs.^init_exponent, 'g', 'linewidth', 3)
end

if t1 > 0
	% middle part
	is=find(t1<time & time<t2);
	xs=time(is);
	ys=meanBound(is);
	yerrs=errBound(is)';
	[cte, exponent, cteStddev, exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
	loglog(xs, cte*xs.^exponent, 'k', 'linewidth', 3)
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


%figure
%clf; hold on;
%plot(time, meanBound, 'b');
%plot(time, meanBound + errBound, 'r');
%plot(time, meanBound - errBound, 'r');
%hold off;

if init_t1 > 0
	[1/init_exponent, init_exponentStddev/init_exponent^2]
end
if t1 > 0
	[1/exponent, exponentStddev/exponent^2]
end
if end_t1 > 0
	[1/end_exponent, end_exponentStddev/end_exponent^2]
end

end


hold off;

% Go back to original time for saving
%time = time / baseLineGamma/gamma;
saveData = [time - time(1), meanBound', errBound'];
save('-ascii', ['./data/meanBoundN',num2str(N)], 'saveData');

