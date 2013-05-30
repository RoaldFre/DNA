addpath generic

%T = 10;
%Ns = [5, 7, 10, 14, 20, 28, 40, 57, 80];
T = 30;
Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57];
clusterSize = 3;
dropLogFactor = 2;
bootstrapSamples = 20;
eps = 1e-5;
dataStartIndex = 2; % ignore first sample at 't = 0' which is actually at 't = dt'
fixOffsetToZero = true;
singleExponent = false;
guessAlpha = 1;
guessBeta = 1;
nu = 2*0.591; %because SQUARED deviation!
scalingFunction = @finiteSizeRescaleWithSize;
squaredDeviation = false;

filenamePrefix = ['./data/hairpinGyradT',num2str(T)];
variableName = 'gyrads';

[clustNs, alphas, alphaErrs, betas, betaErrs] = hairpinFSSclustered(Ns, filenamePrefix, variableName, clusterSize, dropLogFactor, singleExponent, fixOffsetToZero, nu, scalingFunction, squaredDeviation, guessBeta, guessAlpha, eps, bootstrapSamples, dataStartIndex)
