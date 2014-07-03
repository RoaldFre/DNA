% totalEnergy = *negative* of total base pairing energy! (to make it positive for loglog plots)
% function plotHairpinEnergyLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, loadFromCache)
function plotHairpinEnergyLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff)

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

doLinRegs = false;
onlySum = true;
minBoundThreshold = 4;
basesToSum = 1:N;

gamma = str2num(gammaStr);

dir = [basedir,"/T",num2str(T),"C/g",gammaStr,"/N",num2str(N),"/dt",num2str(dt),"/time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_wait",waitStr,"_NXY4"]

filesglob = [dir, '/*_basePairingEnergy'];
if not(loadFromCache)
	if onlySum
		[time, totalEnergy, bound, correctlyBound] = parseBasePairingSamplerSumToArray(filesglob, dropLogFactor, DLdataStartIndex, basesToSum, minBoundThreshold);
		totalEnergy = -totalEnergy; % XXX make it positive for loglog plots!
		% correct for XY base pairs
		bound -= minBoundThreshold;
		correctlyBound -= minBoundThreshold;
		save('-binary', '-z', ['./data/hairpinTotalEnergy',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'totalEnergy', 'bound', 'correctlyBound');
	else
		[time, energy, bound, correctlyBound] = parseBasePairingSamplerToArray(filesglob, dropLogFactor);
		% correct for XY base pairs
		bound -= minBoundThreshold;
		correctlyBound -= minBoundThreshold;
		save('-binary', '-z', ['./data/hairpinEnergy',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'energy', 'bound', 'correctlyBound');
	end
else
	if onlySum
		load(['./data/hairpinTotalEnergy',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)]);
	else
		load(['./data/hairpinEnergy',suff,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)]);
	end
end


if not(onlySum)
	energyThreshold = -1e-20;
	bound = squeeze(sum(energy <= energyThreshold, 2));
	% correct for XY base pairs
	bound = bound - minBoundThreshold;
	nRuns = size(energy, 1);

	if minBoundThreshold > 0
		% remove runs where XY unzipped at least once at any moment in time
		[runsWithUnboundXY, _] = find(bound < 0);
		%sort(runsWithUnboundXY)
		runsWithUnboundXY = unique(runsWithUnboundXY);
		goodRuns = setdiff(1:numel(bound(:,1)), runsWithUnboundXY);

		nRuns = numel(goodRuns)
		bound = bound(goodRuns, :);
		energy = energy(goodRuns, :, :);
	end
	
	% Mean energy per base in function of time
	% meanEnergy(base, time)
	meanEnergy = squeeze(mean(energy(:, :, :), 1)); % TODO investigate this!

	% Total energy of all bases in time for each run
	% totalEnergy(run, time)
	totalEnergy = squeeze(sum(energy(:, basesToSum, :), 2));
	totalEnergy = -totalEnergy; % XXX make it positive for loglog plots!
end

% The data at our disposal at this point:
%   time(timeIndex)
%   totalEnergy(runIndex, timeIndex)


%[totEnergyMean, totEnergyErr] = meanOfSamples(totalEnergy);
%totEnergyErrFrac = 1 + totEnergyErr ./ totEnergyMean;
%hold on;
%loglog(time, totEnergyMean, 'linewidth', 1);
%loglog(time, totEnergyMean .* totEnergyErrFrac, 'r');
%loglog(time, totEnergyMean ./ totEnergyErrFrac, 'r');
%hold off;
%title("energy")
%
%return
%
%figure

[totEnergySqDiff, totEnergySqDiffErr, sqDiffTime] = squaredMeanDeviation(totalEnergy(:,2:end), time(2:end));

totEnergySqDiffErrFrac = 1 + totEnergySqDiffErr ./ totEnergySqDiff;

hold on;
loglog(sqDiffTime, totEnergySqDiff, 'linewidth', 1);
loglog(sqDiffTime, totEnergySqDiff .* totEnergySqDiffErrFrac, 'r');
loglog(sqDiffTime, totEnergySqDiff ./ totEnergySqDiffErrFrac, 'r');
hold off;
title("energy sqDiff")

