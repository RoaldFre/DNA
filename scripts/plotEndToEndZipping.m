function [zippingTime, zippingEndToEnd] = plotEndToEndZipping(N, T, totTime, sweeps, wait, gammaStr, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff)

addpath("generic");

if nargin < 10; loadFromCache = false; end
if nargin < 11; color = 'b'; end
if not(exist('basedir'))
	%basedir = '/home/other/roald/clusterdata/hairpinState_lowTemp_gamma';
	%basedir = '/home/other/roald/clusterdata_properVelocityInit_wrongResizing/hairpinState_lowTemp_gamma';
	%NOTE DANGER DANGER '*' in basedir!! ----------------------v !!
	basedir = '/home/other/roald/clusterdata_properVelocityInit*/hairpinState_lowTemp_gamma';
end
if not(exist('suff')) suff = ''; end


doPlot = true


if ischar(wait)
	waitStr = wait;
else
	waitStr = num2str(wait);
end

if loadFromCache
	load(['./data/hairpinEndToEnd',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)]);
else
	filesglob = [basedir,'/T',num2str(T),'C/g',gammaStr,'/N',num2str(N),'/dt',num2str(dt),'/time',num2str(totTime),'_sweeps',num2str(sweeps, 10),'_wait',waitStr,'_NXY4/*_endToEnd']

	[time, distsAndVectors] = combineNativeRuns(filesglob, false, dropLogFactor, DLdataStartIndex);
	dists = distsAndVectors(:,:,1);
	vectors = distsAndVectors(:,:,2:4);
	save('-binary', '-z', ['./data/hairpinEndToEnd',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'dists', 'vectors');
end

meanDists = mean(dists);
errDists = std(dists) / sqrt(size(dists, 1) - 1);
smoothed = smooth(meanDists, 15);

if doPlot
	if false
		hold on;
		%loglog(time, (mean(meanDists(1:1)) - meanDists).^2)
		%loglog(time, (mean(smoothed(1:1)) - smoothed).^2, 'g', 'linewidth', 3)
		%loglog(time, meanDists)
		%semilogx(time, meanDists)
		plot(time, meanDists)
		%loglog(time, smoothed, 'g', 'linewidth', 3)
		hold off;
	else
		[distSqDev, distSqDevErr, sqDevTime] = squaredMeanDeviation(dists, time, false);
		distSqDevErrFrac = 1 + distSqDevErr ./ distSqDev;
		hold on;
		loglog(sqDevTime, distSqDev);
		loglog(sqDevTime, distSqDev .* distSqDevErrFrac, 'r');
		loglog(sqDevTime, distSqDev ./ distSqDevErrFrac, 'r');
		%plot(sqDevTime, distSqDev);
		%plot(sqDevTime, distSqDev .* distSqDevErrFrac, 'r');
		%plot(sqDevTime, distSqDev ./ distSqDevErrFrac, 'r');
		%plot(sqDevTime, distSqDev + distSqDevErr, 'k');
		%plot(sqDevTime, distSqDev - distSqDevErr, 'k');
		hold off;
	end
end

