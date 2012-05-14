function data = parse2DRuns(directory)

files = glob([directory, "/*"]);

nRuns = numel(files);

for run = 1:nRuns
	data(:,:,run) = load(files{run});
end

