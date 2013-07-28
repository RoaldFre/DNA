addpath generic

generate = true

T = 30;
Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57, 67, 80];
Ns = [20, 24, 28, 34, 40, 48, 57, 67, 80];
clusterSize = numel(Ns);

dropLogFactor = 1;
dataStartIndex = 2; % ignore first sample at 't = 0' which is actually at 't = dt'

opt.bootstrapSamples = -1;
opt.eps = 1e-5;
opt.fixOffsetToZero = true;
opt.singleExponent = false;
opt.twoTimescales = true;
opt.guessAlpha = 1/1.6;
opt.guessAlpha2 = opt.guessAlpha * 0.8;
opt.guessBeta = 1.6;
opt.guessEta = 1;
opt.delta = 1;
opt.squaredDeviation = false;
opt.loglog = true;
%opt.loglog = false;
opt.simulatedAnnealing = true;
opt.SAnt = 100; % iterations between temperature reductions (e.g. 20)
opt.SAns = 25; % iterations between bounds adjustments (e.g. 5)
opt.SArt = 0.95; % temperature reduction factor
opt.T = T;
if opt.twoTimescales
	scalingFunction = @finiteSizeRescaleWithTwoTimes;
else
	scalingFunction = @finiteSizeRescaleWithTime;
end

filenamePrefix = ['./data/hairpinStateT',num2str(T)];
variableName = 'bound';
resultFilePrefix = 'data/stateFSS';

if generate
[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = hairpinFSSclustered(Ns, filenamePrefix, variableName, resultFilePrefix, clusterSize, dropLogFactor, scalingFunction, opt, dataStartIndex)
end

meanStartN = 20;
clf;
finiteSizeResultplots(resultFilePrefix, clusterSize, dropLogFactor, opt, meanStartN)
