% N = number of matching base pairs
function [avgFraction, errFraction, temperatures] = zippingTemperature(filesglob, N)

[averageBoundBasePairs, temperatures, allStates, numMonomers, sampleInterval, timestep, temperature, relaxationTime] = parseHairpinMelting(filesglob);

numRuns = numel(averageBoundBasePairs(1,:));

avgBound = mean(averageBoundBasePairs');
errBound = std(averageBoundBasePairs') / sqrt(max(1, numRuns - 1));

avgFraction = avgBound / N;
errFraction = errBound / N;
