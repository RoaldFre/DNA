function plotHairpinUnzippingState(N, relaxTemp, unzipTemp, totTime, gammaStr, loadFromCache)

addpath("generic");

if nargin < 6
	loadFromCache = false;
end

baseLineGamma = 5e12;
gamma = str2num(gammaStr);

init_t1 = -1; % guard
t1 = -1; %guard
end_t1 = -1; % guard

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


# TODO TODO TODO TODO TODO MERGE DIRECTORY /home/other/roald/clusterdata/hairpinState_lowTemp_gamma/relaxT10C_unzipT90C_time250_relaxFact1.0 !!!

relaxFact = "1.0"; #TODO TODO difference between higher relax??
dir = ["/home/other/roald/clusterdata/hairpinUnzipState/relaxT",num2str(relaxTemp),"C_unzipT",num2str(unzipTemp),"C_time",num2str(totTime),"_relaxFact",relaxFact,"/dt15/g",gammaStr,"/N",num2str(N)]
if not(loadFromCache)
	[bound, time] = parseHairpinStateToArray([dir, '/*.itf11']);
	save('-binary', '-z', [dir,"/hairpinStateN",num2str(N),"_onlyBound"], 'bound', 'time')
else
	load([dir,"/hairpinStateN",num2str(N),"_onlyBound"]);
end


time = time - time(1);

% Rescale time to 'base line'
time = time * baseLineGamma/gamma;


nRuns = numel(bound(:,1))
meanBound = mean(bound);
% Convert to number of *unzipped* BPs! TODO IS THIS THE BEST WAY??
meanBound = meanBound(1) - meanBound;

errBound = std(bound) / sqrt(nRuns - 1);

errFrac = 1 + errBound ./ meanBound;

%clf;
hold on;
loglog(time, meanBound, 'b');
loglog(time, meanBound .* errFrac, 'r');
loglog(time, meanBound ./ errFrac, 'r');

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
if t1 > 0
	[1/exponent, exponentStddev/exponent^2]
end
if end_t1 > 0
	[1/end_exponent, end_exponentStddev/end_exponent^2]
end


% Go back to original time for saving
time = time / baseLineGamma/gamma;
saveData = [time - time(1), meanBound', errBound'];
save('-ascii', ['./data/meanBoundN',num2str(N)], 'saveData');
