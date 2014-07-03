function [zippingTimes, unzippingTimes, zipFromNuclTimes, nucleationTimes, zipping, relax, unzipping] = loadZippingTimeWithNucleationData(filesglob)

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

zipping   = cell(nRuns, 3); %time, bound, state
relax     = cell(nRuns, 3); %time, bound, state
unzipping = cell(nRuns, 3); %time, bound, state

moreWasOn = page_screen_output;
more off;
for run = 1:nRuns
	printf("\rloading file %d of %d", run, nRuns);
	load(files{run});

	zippingTimes(run) = timeTillZipping;
	unzippingTimes(run) = timeTillUnzipping;
	zipFromNuclTimes(run) = timeTillZippingFromNucl;
	nucleationTimes(run) = timeTillZipping - timeTillZippingFromNucl;

	zipping{run, 1} = zippingTime;
	zipping{run, 2} = zippingBound;
	zipping{run, 3} = zippingState;

	relax{run, 1} = relaxTime;
	relax{run, 2} = relaxBound;
	relax{run, 3} = relaxState;

	unzipping{run, 1} = unzippingTime;
	unzipping{run, 2} = unzippingBound;
	unzipping{run, 3} = unzippingState;
end
printf("\n");
if moreWasOn
	more on;
end
