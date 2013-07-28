addpath generic

generate = true;

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';


T = 30;
Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57, 67, 80];
Ns = [28, 34, 40, 48, 57, 67];
Ns = [28, 34, 40, 48];
%clusterSize = 2;
%Ns = [20, 24, 28, 34, 40, 48, 57, 67, 80];
clusterSize = numel(Ns);
%T = 10;
%Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48];
%Ns = [20, 24, 28, 34, 40, 48, 57, 67];
%clusterSize = numel(Ns);
dropLogFactor = 2;
dataStartIndex = 2; % ignore first sample at 't = 0' which is actually at 't = dt'

opt.bootstrapSamples = 3;
opt.eps = 1e-5;
opt.fixOffsetToZero = true;
opt.singleExponent = true;
opt.twoTimescales = false;
opt.guessAlpha = 1/1.6;
opt.guessAlpha2 = opt.guessAlpha * 0.8;
opt.guessBeta = 1.6;
opt.guessEta = 1;
opt.delta = 1;
opt.squaredDeviation = false;
opt.loglog = true;
opt.simulatedAnnealing = false;



opt.singleExponent = false;
opt.rescaleWithSize = true;





filenamePrefix = ['./data/hairpinStateT',num2str(T)];
variableName = 'bound';
graphFilePrefix = ['stateT',num2str(T),'FSS'];
resultFilePrefix = ['data/',graphFilePrefix];

plotopt.T = T;
plotopt.graphFilePrefix = graphFilePrefix;
plotopt.measurement = 'the number of bound base pairs';
plotopt.quantity = '\nb';
plotopt.delta = 1;

if generate
	[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = ...
	hairpinFSSclustered(Ns, filenamePrefix, variableName, resultFilePrefix, clusterSize, dropLogFactor, opt, plotopt, dataStartIndex)
else
	resultFile = hairpinFSSclusteredFilename(resultFilePrefix, clusterSize, dropLogFactor, opt)
	load(resultFile)
end

meanStartN = 20;
clf;
%finiteSizeResultplots(resultFilePrefix, clusterSize, dropLogFactor, opt, meanStartN)

makeFiniteSizeResultplot(clustNs, clustPs,    clustPErrs,    false, opt, plotopt, meanStartN)
makeFiniteSizeResultplot(clustNs, clustQuals, clustQualErrs, true,  opt, plotopt, meanStartN)






% collapse everything

% collapse starting from meanStartN
 
