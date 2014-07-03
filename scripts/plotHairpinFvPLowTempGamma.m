function plotHairpinFvPLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, dropLogFactor, DLdataStartIndex, loadFromCache)

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

filesglob = [dir, '/*.itf11_forceVelP'];
if not(loadFromCache)
	columns = 1 + 4*(0:(2*N+10)); %TODO hardcoded for length 2N + 2x4XY + 4loop monomers!! (-1 because sampler cannot sample one monomer at the beginnin (or end, I forgot) of the strand)

	[time, data, files] = combineNativeRuns(filesglob, true, dropLogFactor, DLdataStartIndex);
	forcePara = data(:,:,columns);
	forcePerp = data(:,:,columns + 1);
	velPara = data(:,:,columns + 2);
	velPerp = data(:,:,columns + 3);
	save('-binary', '-z', ['./data/hairpinFvPForceT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'forcePara', 'forcePerp', 'files');
	save('-binary', '-z', ['./data/hairpinFvPVelT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'velPara', 'velPerp', 'files');
	clear data;
end


legendName = 'phosphate';
dataStartIndex = 2;

plotHairpinTension('./data/hairpinFvPForce', 'forcePara', legendName, T, N, dropLogFactor, dataStartIndex, 50)
plotHairpinTension('./data/hairpinFvPForce', 'forcePerp', legendName, T, N, dropLogFactor, dataStartIndex, 50)
plotHairpinTension('./data/hairpinFvPVel', 'velPara', legendName, T, N, dropLogFactor, dataStartIndex, 50)
plotHairpinTension('./data/hairpinFvPVel', 'velPerp', legendName, T, N, dropLogFactor, dataStartIndex, 50)

