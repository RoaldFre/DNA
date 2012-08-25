function data = parse2DRuns(filesglob)

#files = glob([directory, "/*"]);
files = glob(filesglob);

nRuns = numel(files);

for run = 1:nRuns
	data(:,:,run) = load(files{run});
end

