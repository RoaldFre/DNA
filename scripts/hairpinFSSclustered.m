% Ns: list of system sizes
% filenamePrefix:   Prefix of the data files. The full datafile must be of the form
%                   <filenamePrefix>'N'<N>, with <N> the system size as per Ns, e.g.: 'dataN10'
% resultFilePrefix: Prefix of the file name where the resutling data will be saved to.
% variableName:     The name of the variable that holds the data in the data file
%                   each row in that matrix is assumed to be an independent run
%                   there is also assumed to be a 'time' variable with corresponding time values!
% clusterSize:      Do FSS on clusters of <clusterSize> adjacent sytem sizes
% dropLogFactor;    Transform linear time samples to log space time samples with dropDataLog() with this factor
% singleExponent:   If 'true', then do single exponent fit based on nu and guessBeta
function [clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = hairpinFSSclustered(Ns, filenamePrefix, variableName, resultFilePrefix, clusterSize, dropLogFactor, opt, plotopt, dataStartIndex)

addpath generic

more off

if nargin < 7; opt = struct(); end
if nargin < 8; plotopt = struct(); end
if nargin < 9; dataStartIndex = 1; end

opt.Ns = Ns;
numNs = numel(Ns);

% Full data of all runs
times  = cell(numNs, 1);
datas  = cell(numNs, 1);
timesDecim  = cell(numNs, 1);
datasDecim  = cell(numNs, 1);


for i = 1:numNs
	filename = [filenamePrefix,'N',num2str(Ns(i))];
	load(filename);

	if exist('bound')
		[runsWithUnboundXY, _] = find(bound < 0);
		runsWithUnboundXY = unique(runsWithUnboundXY)
		goodRuns = setdiff(1:numel(bound(:,1)), runsWithUnboundXY);
		nRuns = numel(goodRuns);
		bound = bound(goodRuns, :);
		%TODO take goodRuns  from 'variableName'!
	end

	if not(exist(variableName))
		error(["Can't find the variable with name '",variableName,"' in the data file '",filename,"'!"]);
	end

	times{i}  = time';
	eval(['datas{i} = ',variableName,';'])

	times{i} = times{i}(dataStartIndex:end);
	datas{i} = datas{i}(:, dataStartIndex:end);

	timesDecim{i}  = dropDataLogspace(times{i}, dropLogFactor);
	datasDecim{i} = dropDataLogspace(datas{i}, dropLogFactor);
end


% Normal (non-resampled) fit on clusters of 'adjacent' N values
numClusters = numNs - clusterSize + 1;
clustQuals = zeros(numClusters, 1);
clustQualErrs = zeros(numClusters, 1);
opts = cell(numClusters, 1);
for i = 1 : numClusters
	selectedNs = Ns(i : i+clusterSize-1);
	clustXsDecim = timesDecim(i : i+clusterSize-1);
	clustYsDecim = datasDecim(i : i+clusterSize-1);

	[exponentsAndOffsets, quality, exponentsAndOffsetsErr, qualityErr, newOpt] = ...
		finiteSizeScaling(selectedNs, clustXsDecim, clustYsDecim, [], opt, plotopt);

	
	clustQuals(i) = quality;
	clustQualErrs(i) = qualityErr;
	clustPs(i,:) = exponentsAndOffsets';
	clustPErrs(i,:) = exponentsAndOffsetsErr';
	opts{i} = newOpt;

	%clustNs
	%[xOffsets yOffsets]
end

% TODO QUICK HACK
opt = opts{1}; 
opt.Ns = Ns;

clustNs = zeros(numClusters, 1);
for i = 1 : numClusters
	clustNs(i) = mean(Ns(i : i+clusterSize-1), 'g'); % geometric mean
end

resultFile = hairpinFSSclusteredFilename(resultFilePrefix, clusterSize, dropLogFactor, opt)
save('-z', '-binary', resultFile, 'clustNs', 'clustPs', 'clustPErrs', 'clustQuals', 'clustQualErrs', 'Ns', 'clusterSize', 'opts', 'opt');


clf
finiteSizeResultplots(resultFilePrefix, clusterSize, dropLogFactor, opt);
