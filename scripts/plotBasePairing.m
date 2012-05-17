function plotBasePairing(filename)
data = load(filename);
dataSize = size(data);
nSamples = dataSize(1);
nMonomers = dataSize(2) - 3; 
	%last two columns are just total number of (correctly) bonded 
	%pairs, the one before that holds the time

time = data(:, nMonomers + 1);
config = data(:, 1:nMonomers); %monomer configuration



imagesc(time, 1:nMonomers, data(1:nSamples, 1:nMonomers)');
colormap(bone(2));
