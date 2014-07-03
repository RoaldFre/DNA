% FOR MIN ZIP, LOW TEMP
function plotZippingBoundNuclThreshold(N, minTime, nuclBoundBPs)

addpath "generic"

dir="/home/other/roald/clusterdata/hairpinFormationLowTemp_native/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/dt15/zipT10C_unzipT90C_allowUnb0.2_nucl0.2/MCsweeps100000_relaxFactor2/"



glob = [dir,'N',num2str(N),'_minZipTime',num2str(minTime),'/*itf11'];

[zippingt, unzippingt, zipFromNuclt, nucleationt, preZipping, zipping, unzipping] = loadZippingTimeWithNucleationDataCustomThreshold(glob, nuclBoundBPs, N*0.8);

numRuns = numel(zippingt);
[time, zippingBound] = getIntersection(zipping);


% XXX Correcting for the already "nucleation bounded" base pairs to get proper power law without bias
zippingBound = zippingBound - nuclBoundBPs;
% TODO this makes a pretty significant difference on the exponents. What to do? Factor it into the uncertainty?? TODO TODO
% TODO: not adding this seems to get a better fit at the relaxation portion (although the initial couple of data points 'curl upwards' -- adding the correction fixes the initial data points, but the entire relaxation phase gets a slight curve in it!)



init_t1 = -1; % guard
end_t1 = -1; % guard
if N == 7
	t1=4.0e-9;
	t2=1.8e-8;
elseif N == 10
	%init_t1 = 2e-10;
	%init_t2 = 3e-9;
	%t1=4.0e-9;
	%t2=1.8e-8;
	t1 = 2e-10;
	t2 = 5.5e-9;
elseif N == 14
	init_t1 = 2.5e-10;
	init_t2 = 3e-9;
	t1=4.0e-9;
	t2=1.8e-8;
elseif N == 20
	init_t1 = 1.3e-10;
	init_t2 = 5.1e-9;
	t1 = 1.0e-8;
	t2 = 4.4e-8;
elseif N == 28
	init_t1 = 5.3e-10;
	init_t2 = 1.8e-8;
	t1 = 3.3e-8;
	t2 = 9.8e-8;
elseif N == 40
	init_t1 = 6.0e-10;
	init_t2 = 2.5e-8;
	t1 = 7.5e-8;
	t2 = 2.1e-7;
	%end_t1 = 1.6e-7;
	%end_t2 = 2.0e-7;
elseif N == 57
	init_t1 = 1.0e-10;
	init_t2 = 3.3e-8;
	t1 = 5.0e-8;
	t2 = 3.4e-7;
	end_t1 = 3.6e-7;
	end_t2 = 5.0e-7;
else
	error blah
end



meanBound = mean(zippingBound);
errBound = std(zippingBound) / sqrt(numRuns - 1);

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



return



[unzipTime, unzippingBound] = getIntersection(unzipping);
unzippingBound = N*0.83 - unzippingBound;
meanUnzipBound = mean(unzippingBound);
errUnzipBound = std(unzippingBound) / sqrt(numRuns - 1);
errUnzipFrac = 1 + errUnzipBound ./ meanUnzipBound;

clf
hold on;
loglog(unzipTime, meanUnzipBound, 'b');
loglog(unzipTime, meanUnzipBound .* errUnzipFrac, 'r');
loglog(unzipTime, meanUnzipBound ./ errUnzipFrac, 'r');
xlabel("time (s)")
ylabel("average number of bound base pairs");

[unzip_cte, unzip_exponent, unzip_cteStddev, unzip_exponentStddev] = loglogRegression(unzipTime(2:end), meanUnzipBound(2:end), 1e5, 0.6, errUnzipBound(2:end))
loglog(unzipTime, unzip_cte*unzipTime.^unzip_exponent, 'g', 'linewidth', 3)
hold off;

[1/unzip_exponent, unzip_exponentStddev/unzip_exponent^2]





return


figure
imagesc(meanZippingState');
figure
hold on
plot(meanZippingState(1,:), 'r')
plot(meanZippingState(round(end*0.25),:), 'g')
plot(meanZippingState(round(end*0.5),:), 'b')
hold off
