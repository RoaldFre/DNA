addpath generic

generate = true

T = 30;
Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57, 67, 80];
Ns = [20, 24, 28, 34, 40, 48, 57, 67, 80];
%T = 10;
%Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48];
clusterSize = 2;
dropLogFactor = 2;
dataStartIndex = 2; % ignore first sample at 't = 0' which is actually at 't = dt'

scalingFunction = @finiteSizeRescaleWithSize;
opt.singleExponent = false;

%scalingFunction = @finiteSizeRescaleWithTime;
%opt.singleExponent = true;

opt.bootstrapSamples = 20;
opt.eps = 1e-5;
opt.fixOffsetToZero = true;
opt.guessAlpha = 1;
opt.guessBeta = 1;
opt.delta = 0.587597;
opt.squaredDeviation = false;
opt.twoTimescales = false;
opt.loglog = true;
opt.simulatedAnnealing = false;

filenamePrefix = ['./data/hairpinGyradT',num2str(T)];
variableName = 'gyrad';
resultFilePrefix = 'data/gyradFSS';

if generate
[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = hairpinFSSclustered(Ns, filenamePrefix, variableName, resultFilePrefix, clusterSize, dropLogFactor, scalingFunction, opt, dataStartIndex)
end

meanStartN = 20;
clf;
finiteSizeResultplots(resultFilePrefix, clusterSize, dropLogFactor, opt, meanStartN)
