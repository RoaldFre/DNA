function plotBasePairing(filename, plotState, plotNumBasePairs, numBPsDecimateFactor)

addpath("octave-forge");

if nargin == 1
	plotState = true;
	plotNumBasePairs = false;
	numBPsDecimateFactor = 1;
elseif nargin == 2
	plotNumBasePairs = false;
	numBPsDecimateFactor = 1;
elseif nargin == 3
	numBPsDecimateFactor = 1;
end

data = load(filename);
dataSize = size(data);
nSamples = dataSize(1);
nMonomers = dataSize(2) - 3; 
	%last two columns are just total number of (correctly) bonded 
	%pairs, the one before that holds the time

config         = data(:, 1:nMonomers); %monomer configuration
time           = data(:, nMonomers + 1);
allBound       = data(:, nMonomers + 2);
correctlyBound = data(:, nMonomers + 3);

if plotNumBasePairs
	q = numBPsDecimateFactor;
	plot(decimate(time, q), decimate(allBound, q), 'k', ...
	     decimate(time, q), decimate(correctlyBound,q), 'b');
end

if plotState
	imagesc(time, 1:nMonomers, data(1:nSamples, 1:nMonomers)');
	colormap(bone(2));
end
