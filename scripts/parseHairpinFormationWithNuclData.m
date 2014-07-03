% if numZippable < 0, then the *fromNucleation data is not calculated
function [timesTillZipping, timesTillUnzipping, timesTillZippingFromNucleation, timesTillNucleation] = parseHairpinFormationWithNuclData(filesglob);

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

moreWasOn = page_screen_output;
more off;
for run = 1:nRuns
	printf("\rloading file %d of %d", run, nRuns);
	load(files{run});
	timesTillZipping(run) = timeTillZipping;
	timesTillUnzipping(run) = timeTillUnzipping;
	timesTillZippingFromNucleation(run) = timeTillZippingFromNucl;
	timesTillNucleation(run) = timeTillZipping - timeTillZippingFromNucl;
end
printf("\n");
if moreWasOn
	more on;
end
