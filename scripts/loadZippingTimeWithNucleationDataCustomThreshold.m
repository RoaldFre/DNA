function [zippingTimes, unzippingTimes, zipFromNuclTimes, nucleationTimes, preZipping, zipping, unzipping] = loadZippingTimeWithNucleationDataCustomThreshold(filesglob, nuclBoundBPs, zippedBoundBPs)

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

preZipping = cell(nRuns, 3); %time, bound, state
zipping    = cell(nRuns, 3); %time, bound, state
unzipping  = cell(nRuns, 3); %time, bound, state

moreWasOn = page_screen_output;
more off;
for run = 1:nRuns
	printf("\rloading file %d of %d", run, nRuns);
	load(files{run});

	if zippingTime(end) == 	relaxTime(1)
		relaxTime  = relaxTime(2:end);
		relaxBound = relaxBound(2:end);
		relaxState = relaxState(2:end, :);
	end

	% Combine zipping and relax into one zipping phase
	zippingPhaseTime  = [zippingTime;  relaxTime];
	zippingPhaseBound = [zippingBound; relaxBound];
	zippingPhaseState = [zippingState; relaxState];

	[zippingTime, nucleationTime, zipFromNuclTime, zipIndex, nucleationIndex] = getZippingTimesFromBound(zippingPhaseTime, zippingPhaseBound, nuclBoundBPs, zippedBoundBPs);

	zippingTimes(run) = zippingTime;
	nucleationTimes(run) = nucleationTime;
	zipFromNuclTimes(run) = zipFromNuclTime;
	unzippingTimes(run) = timeTillUnzipping;

	preZipping{run, 1} = zippingPhaseTime(1 : nucleationIndex-1);
	preZipping{run, 2} = zippingPhaseBound(1 : nucleationIndex-1);
	preZipping{run, 3} = zippingPhaseState(1 : nucleationIndex-1, :);

	zipping{run, 1} = zippingPhaseTime(nucleationIndex : end);
	zipping{run, 2} = zippingPhaseBound(nucleationIndex : end);
	zipping{run, 3} = zippingPhaseState(nucleationIndex : end, :);

	unzipping{run, 1} = unzippingTime;
	unzipping{run, 2} = unzippingBound;
	unzipping{run, 3} = unzippingState;
end
printf("\n");
if moreWasOn
	more on;
end
