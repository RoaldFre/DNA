function plotHairpinFvfPLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, dropLogFactor, DLdataStartIndex, loadFromCache)

addpath("generic");


if nargin < 10
	loadFromCache = false;
end

if ischar(wait)
	waitStr = wait;
else
	waitStr = num2str(wait);
end



%basedir = '/home/other/roald/clusterdata/hairpinState_lowTemp_gamma';
%basedir = '/home/other/roald/clusterdata_properVelocityInit_wrongResizing/hairpinState_lowTemp_gamma';
%NOTE DANGER DANGER '*' in basedir!! ----------------------v !!
basedir = '/home/other/roald/clusterdata_properVelocityInit*/hairpinState_lowTemp_gamma';
dir = [basedir,"/T",num2str(T),"C/g",gammaStr,"/N",num2str(N),"/dt",num2str(dt),"/time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_wait",waitStr,"_NXY4"]

filesglob = [dir, '/*.itf11_forceVelFricP'];
if not(loadFromCache)
	columns = 1 + 3*(0:(2*N+11)); %TODO hardcoded for length 2N + 2x4XY + 4loop monomers!!

	%% Split to force, velocity, friction because it is too big otherwise
	%[time, force, files] = combineNativeRuns(filesglob, true, dropLogFactor, DLdataStartIndex, columns);
	%save('-binary', '-z', ['./data/hairpinForceT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'force', 'files');
	%clear force
	%[time, vel, files] = combineNativeRuns(filesglob, true, dropLogFactor, DLdataStartIndex, columns + 1);
	%save('-binary', '-z', ['./data/hairpinVelT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'vel', 'files');
	%clear vel
	%[time, fric, files] = combineNativeRuns(filesglob, true, dropLogFactor, DLdataStartIndex, columns + 2);
	%save('-binary', '-z', ['./data/hairpinFricT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'fric', 'files');
	%clear fric

	[time, data, files] = combineNativeRuns(filesglob, true, dropLogFactor, DLdataStartIndex);
	force = data(:,:,columns);
	vel = data(:,:,columns + 1);
	fric = data(:,:,columns + 2);
	save('-binary', '-z', ['./data/hairpinPForceT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'force', 'files');
	save('-binary', '-z', ['./data/hairpinPVelT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'vel', 'files');
	save('-binary', '-z', ['./data/hairpinPFricT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'fric', 'files');

	fric2 = fric .* (force.^2) ./ (vel.^2); % I FUCKED UP IN THE SAMPLER! >_<
	clear vel fric force
	save('-binary', '-z', ['./data/hairpinPFric2T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'fric2', 'files');
end




%load(['./data/hairpinPForceT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)])
%load(['./data/hairpinPVelT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)])
%load(['./data/hairpinPFricT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)])
%fric2 = fric .* (force.^2) ./ (vel.^2); % I FUCKED UP IN THE SAMPLER! >_<
%clear vel fric force
%save('-binary', '-z', ['./data/hairpinPFric2T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'fric2', 'files');




legendName = 'phosphate';
dataStartIndex = 1;

%plotHairpinTension('./data/hairpinPVel', 'vel', legendName, T, N, dropLogFactor, dataStartIndex, 1)
%plotHairpinTension('./data/hairpinPVel', 'vel', legendName, T, N, dropLogFactor, dataStartIndex, 10)
plotHairpinTension('./data/hairpinPVel', 'vel', legendName, T, N, dropLogFactor, dataStartIndex, 100)
%plotHairpinTension('./data/hairpinPForce', 'force', legendName, T, N, dropLogFactor, dataStartIndex, 1)
%plotHairpinTension('./data/hairpinPForce', 'force', legendName, T, N, dropLogFactor, dataStartIndex, 10)
plotHairpinTension('./data/hairpinPForce', 'force', legendName, T, N, dropLogFactor, dataStartIndex, 100)
%plotHairpinTension('./data/hairpinPFric', 'fric', legendName, T, N, dropLogFactor, dataStartIndex, 1)
%plotHairpinTension('./data/hairpinPFric', 'fric', legendName, T, N, dropLogFactor, dataStartIndex, 10)
plotHairpinTension('./data/hairpinPFric', 'fric', legendName, T, N, dropLogFactor, dataStartIndex, 100)

plotHairpinTension('./data/hairpinPFric2', 'fric2', legendName, T, N, dropLogFactor, dataStartIndex, 100)
