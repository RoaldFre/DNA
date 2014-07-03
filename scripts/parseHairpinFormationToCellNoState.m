% Parse all runs to a cell, don't store the full state
function data = parseHairpinFormationToCellNoState(filesglob);

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

data = cell(nRuns, 9);

for run = 1:nRuns
	timesTillZipping(run) = load(files{run}, "timeTillZipping").timeTillZipping;
	timesTillUnzipping(run) = load(files{run}, "timeTillUnzipping").timeTillUnzipping;

	[zippingTime,   zippingState,   zippingBound ...
         relaxTime,     relaxState,     relaxBound, ...
         unzippingTime, unzippingState, unzippingBound] ...
			 = parseHairpinFormationState(files{run});
	data{run, 1} = zippingTime;
	data{run, 2} = [];
	data{run, 3} = zippingBound;
	data{run, 4} = relaxTime;
	data{run, 5} = [];
	data{run, 6} = relaxBound;
	data{run, 7} = unzippingTime;
	data{run, 8} = [];
	data{run, 9} = unzippingBound;
end

