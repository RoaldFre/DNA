function [timesTillZipping, timesTillUnzipping, timesTillZippingFromNucleation] = parseHairpinFormation(filesglob);

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

for run = 1:nRuns
	timesTillZipping(run) = load(files{run}, "timeTillZipping").timeTillZipping;
	timesTillUnzipping(run) = load(files{run}, "timeTillUnzipping").timeTillUnzipping;

	% time starting from nucleation:
	[zippingTime, zippingState, zippingBound] = parseHairpinFormationState(files{run});
	N = numel(zippingState(1,:));
	nucleationThreshold = round(N * 0.10);
	nucleationIndex = find(zippingBound >= nucleationThreshold, 1, 'first');
	timesTillZippingFromNucleation(run) = zippingTime(end) - zippingTime(nucleationIndex);
end

