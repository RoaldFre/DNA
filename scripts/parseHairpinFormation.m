function [timesTillZipping, timesTillUnzipping, temperature, temperatures, data, timestep, allowedUnboundBPs, allowedBoundBps] = parseHairpinFormation(filesglob, alsoLoadFullData);

if (nargin < 1)
	error("Not enough arguments!");
end
if (nargin == 1)
	alsoLoadFullData = false;
	allData = 0;
end

files = glob(filesglob);

nRuns = numel(files);

load(files{1}); % for constants, temperatures, etc
clear data;


for run = 1:nRuns
	timesTillZipping(run) = load(files{run}, "timeTillZipping").timeTillZipping;
	timesTillUnzipping(run) = load(files{run}, "timeTillUnzipping").timeTillUnzipping;
	if (alsoLoadFullData)
		allData(:,:,:,run) = load(files{run}, "data").data;
	end
end

