function [timesTillZipping, timesTillUnzipping] = parseHairpinFormation(filesglob);

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

for run = 1:nRuns
	timesTillZipping(run) = load(files{run}, "timeTillZipping").timeTillZipping;
	timesTillUnzipping(run) = load(files{run}, "timeTillUnzipping").timeTillUnzipping;
end

