function plotHairpinGyradLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff)

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


if ischar(wait)
	waitStr = wait;
else
	waitStr = num2str(wait);
end



dir = [basedir,"/T",num2str(T),"C/g",gammaStr,"/N",num2str(N),"/dt",num2str(dt),"/time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_wait",waitStr,"_NXY4"]

filesglob = [dir, '/*_gyrationRadius'];
if not(loadFromCache)
	[time, gyrad] = combineNativeRuns(filesglob, false, dropLogFactor, DLdataStartIndex);
	% correct for XY base pairs
	save('-binary', '-z', ['./data/hairpinGyrad',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'gyrad');
else
	load(['./data/hairpinGyrad',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)]);
end


nRuns = size(gyrad, 1);

if false
	meanGyrad = squeeze(mean(gyrad, 1));
	meanGyradErr = squeeze(std(gyrad, 0, 1)) / sqrt(nRuns - 1);
	meanGyradErrFrac = 1 + meanGyradErr ./ meanGyrad;

	hold on;
	loglog(time, meanGyrad);
	loglog(time, meanGyrad .* meanGyradErrFrac, 'r');
	loglog(time, meanGyrad ./ meanGyradErrFrac, 'r');
	hold off;

else

	[gyradSqDev, gyradSqDevErr, sqDevTime] = squaredMeanDeviation(gyrad, time);
	gyradSqDevErrFrac = 1 + gyradSqDevErr ./ gyradSqDev;
	hold on;
	%loglog(sqDevTime, gyradSqDev);
	%loglog(sqDevTime, gyradSqDev .* gyradSqDevErrFrac, 'r');
	%loglog(sqDevTime, gyradSqDev ./ gyradSqDevErrFrac, 'r');
	plot(sqDevTime, gyradSqDev);
	plot(sqDevTime, gyradSqDev .* gyradSqDevErrFrac, 'r');
	plot(sqDevTime, gyradSqDev ./ gyradSqDevErrFrac, 'r');
	plot(sqDevTime, gyradSqDev + gyradSqDevErr, 'k');
	plot(sqDevTime, gyradSqDev - gyradSqDevErr, 'k');
	hold off;
end

sleep(1e-9);
