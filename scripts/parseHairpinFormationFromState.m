function [timesTillZipping, timesTillUnzipping] = parseHairpinFormationFromState(filesglob);

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

%zipfig = figure;
%unzipfig = figure;
hold on;

for run = 1:nRuns
	timesTillZipping(run) = load(files{run}, "timeTillZipping").timeTillZipping;
	timesTillUnzipping(run) = load(files{run}, "timeTillUnzipping").timeTillUnzipping;

	% time starting from nucleation:
	[zippingTime,   zippingState,   zippingBound, ...
         relaxTime,     relaxState,     relaxBound, ...
         unzippingTime, unzippingState, unzippingBound] ...
	  			= parseHairpinFormationState(files{run});
	N = numel(zippingState(1,:));

	plot(zippingTime - zippingTime(1), zippingBound);
	%plot(unzippingTime - unzippingTime(1), unzippingBound);
end

hold off;
