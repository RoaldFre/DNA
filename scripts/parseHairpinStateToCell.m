% Parse all runs to a cell
function data = parseHairpinStateToCell(filesglob, alsoSaveFullState);

if nargin < 2
	alsoSaveFullState = true;
end

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

data = cell(nRuns, 3);

for run = 1:nRuns
	runData = load(files{run});
	data{run, 1} = runData(:,1);
	data{run, 2} = runData(:,2);
	if alsoSaveFullState
		data{run, 3} = runData(:,3:end);
	else
		data{run, 3} = [];
	end
end

