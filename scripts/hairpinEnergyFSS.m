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

scalingFunction = @finiteSizeRescaleWithTime;

opt.bootstrapSamples = 20;
opt.eps = 1e-5;
opt.fixOffsetToZero = true;
opt.singleExponent = true;
opt.guessAlpha = 1/1.6;
opt.guessAlpha2 = opt.guessAlpha * 0.8;
opt.guessBeta = 1.6;
opt.delta = 1;
%opt.squaredDeviation = false;
opt.squaredDeviation = true;
opt.twoTimescales = false;
opt.loglog = true;
opt.simulatedAnnealing = false;

filenamePrefix = ['./data/hairpinTotalEnergyT',num2str(T)];
variableName = 'totalEnergy';
resultFilePrefix = 'data/totalEnergyFSS';

if generate
[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = hairpinFSSclustered(Ns, filenamePrefix, variableName, resultFilePrefix, clusterSize, dropLogFactor, scalingFunction, opt, dataStartIndex)
end

meanStartN = 20;
clf;
finiteSizeResultplots(resultFilePrefix, clusterSize, dropLogFactor, opt, meanStartN)
