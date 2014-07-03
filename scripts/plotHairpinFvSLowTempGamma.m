function plotHairpinFvSLowTempGamma(N, T, totTime, sweeps, wait, gammaStr, dt, dropLogFactor, DLdataStartIndex, loadFromCache)

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

filesglob = [dir, '/*.itf11_forceVelS'];
if not(loadFromCache)
	columns = 1 + 6*(0:(2*N+10)); %TODO hardcoded for length 2N + 2x4XY + 4loop monomers!! (-1 because sampler cannot sample one monomer at the beginnin (or end, I forgot) of the strand)

	[time, data, files] = combineNativeRuns(filesglob, true, dropLogFactor, DLdataStartIndex);
	forcePara = data(:,:,columns);
	forceBase = data(:,:,columns + 1);
	forcePerp = data(:,:,columns + 2);
	velPara = data(:,:,columns + 3);
	velBase = data(:,:,columns + 4);
	velPerp = data(:,:,columns + 5);
	save('-binary', '-z', ['./data/hairpinFvSForceT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'forcePara', 'forceBase', 'forcePerp', 'files');
	save('-binary', '-z', ['./data/hairpinFvSVelT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)], 'time', 'velPara', 'velBase', 'velPerp', 'files');
end

legendName = '';
dataStartIndex = 2;

blockAvg = 50;
%blockAvg = 1;
particleType = 'sugars';
observable = 'force on';
yfactor = 1e12;
plotHairpinTension('./data/hairpinFvSForce', 'forcePara', legendName, T, N, dropLogFactor, dataStartIndex, blockAvg, particleType, '$\expect{F_\mathrm{t}(t)}$ (pN)', 'tangential $\ut$', observable, yfactor)
plotHairpinTension('./data/hairpinFvSForce', 'forceBase', legendName, T, N, dropLogFactor, dataStartIndex, blockAvg, particleType, '$\expect{F_\mathrm{b}(t)}$ (pN)', 'base $\ub$', observable, yfactor)
plotHairpinTension('./data/hairpinFvSForce', 'forcePerp', legendName, T, N, dropLogFactor, dataStartIndex, blockAvg, particleType, '$\expect{F_\mathrm{p}(t)}$ (pN)', 'perpendicular $\up$', observable, yfactor)
observable = 'velocity of';
yfactor = 1;
plotHairpinTension('./data/hairpinFvSVel', 'velPara', legendName, T, N, dropLogFactor, dataStartIndex, blockAvg, particleType, '$\expect{v_\mathrm{t}(t)}$ (m/s)', 'tangential $\ut$', observable, yfactor)
plotHairpinTension('./data/hairpinFvSVel', 'velBase', legendName, T, N, dropLogFactor, dataStartIndex, blockAvg, particleType, '$\expect{v_\mathrm{b}(t)}$ (m/s)', 'base $\ub$', observable, yfactor)
plotHairpinTension('./data/hairpinFvSVel', 'velPerp', legendName, T, N, dropLogFactor, dataStartIndex, blockAvg, particleType, '$\expect{v_\mathrm{p}(t)}$ (m/s)', 'perpendicular $\up$', observable, yfactor)

