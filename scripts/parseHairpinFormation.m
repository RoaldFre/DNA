% if numZippable < 0, then the *fromNucleation data is not calculated
function [timesTillZipping, timesTillUnzipping, timesTillZippingFromNucleation, timesTillNucleation] = parseHairpinFormation(filesglob, numZippable);

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

for run = 1:nRuns
	timesTillZipping(run) = load(files{run}, "timeTillZipping").timeTillZipping;
	timesTillUnzipping(run) = load(files{run}, "timeTillUnzipping").timeTillUnzipping;

	if numZippable >= 0
		[zippingTime, zippingState, zippingBound] = parseHairpinFormationState(files{run});
		nucleationThreshold = round(numZippable * 0.40);
		nucleationIndex = find(zippingBound <= nucleationThreshold, 1, 'last');
		if isempty(nucleationIndex)
			%figure
			%plot(zippingBound)
			disp(["WARNING! Couldn't find nucleation in zipping! Threshold: ",num2str(nucleationThreshold),". Probably already had some zipping before quenching (so in the initial relaxation phase)! Assuming nucleation happened right at the time of the quench..."]);
			nucleationIndex = 1;
		end
		timesTillZippingFromNucleation(run) = zippingTime(end) - zippingTime(nucleationIndex);
		timesTillNucleation(run) = zippingTime(nucleationIndex) - zippingTime(1);
	end
end

