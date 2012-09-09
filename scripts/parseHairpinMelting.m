function [averageBoundBasePairs, temperatures, allStates, numMonomers, sampleInterval, timestep, temperature, relaxationTime] = parseHairpinMelting(filesglob, alsoLoadFullData);

if (nargin < 1)
	error "Not enough arguments!"
end
if (nargin == 1)
	alsoLoadFullData = false;
	allStates = 0;
end

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

load(files{1}); % for constants, temperatures, etc
clear hairpinState;


for run = 1:nRuns
	averageBoundBasePairs(:,run) = load(files{run}, "averageBoundBasePairs").averageBoundBasePairs;
	if (alsoLoadFullData)
		allStates(:,:,:,run) = load(files{run}, "hairpinState").hairpinState;
	end
end

