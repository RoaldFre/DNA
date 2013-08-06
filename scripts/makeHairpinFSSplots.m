% e.g.:
%
% variableName = 'totalEnergy';
% filenamePrefixBase = 'hairpinTotalEnergy';
% graphFilePrefixBase = 'energy';
% 
% plotopt.destDir = '../../thesis/latex/images/plots';
% plotopt.relImgDir = 'images/plots';
% plotopt.measurement = 'the base pairing energy';
% plotopt.quantity = '\Vbp';
% plotopt.delta = 1; or '\nu'
%
% dropLogFactor = 2;
%
% sqDevWithTime & sqDevWithSize: can be an array (ig [true, false] for both, or [] for no plots of that scaling [time or size])
function makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, dropLogFactor, delta, sqDevWithTime, sqDevWithSize, minN)
addpath generic

if nargin < 9; minN = 0; end;

generate = true;

dataStartIndex = 2; % ignore first sample at 't = 0' which is actually at 't = dt'
meanStartN = 20;

opt.eps = 1e-5;
opt.fixOffsetToZero = true;
opt.twoTimescales = false;
opt.loglog = true;

opt.bootstrapSamples = 15;
opt.bootstrapSamples = 10;
opt.simulatedAnnealing = false;
%opt.bootstrapSamples = 5;
%opt.simulatedAnnealing = true;

%opt.guessAlpha = delta/1.6; %should be x2 if squared dev...
opt.guessAlpha = 1; % just keep it simple
opt.guessAlpha2 = opt.guessAlpha * 0.8;
opt.guessBeta = 1.6;
opt.guessEta = 1;

opt.delta = delta;

baseOpt = opt;
basePlotopt = plotopt;

for T = [10 30]
    #BIT HACKY, but we need to get the temperature in the filename to discriminate the data sets...
    plotopt = basePlotopt;
    plotopt.resultFilePrefix = [plotopt.resultFilePrefix,'T',num2str(T)];

    for singleExponent = [false, true]

	if T == 10
		Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48];
		%Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40];
	elseif T == 30
		%Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57, 67, 80];
		Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48];
	else
		error "Unknown T"
	end

	Ns = Ns(find(Ns >= minN));

	plotopt.T = T;
	baseOpt.T = T;

	clusterSize = 2;

	filenamePrefix = ['./data/', filenamePrefixBase, 'T' ,num2str(T)];
	graphFilePrefix = hairpinFSSclusteredFilenamePrefix(graphFilePrefixBase, clusterSize, dropLogFactor, T);
	resultFilePrefix = ['data/',graphFilePrefix];
	plotopt.graphFilePrefix = graphFilePrefix;


	for squaredDeviation = sqDevWithTime
		% RESCALE WITH TIME
		opt = baseOpt;
		opt.squaredDeviation = squaredDeviation
		opt.singleExponent = singleExponent;
		opt.rescaleWithSize = false;
		if generate
			[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = ...
			hairpinFSSclustered(Ns, filenamePrefix, variableName, resultFilePrefix, clusterSize, dropLogFactor, opt, plotopt, dataStartIndex)
		else
			resultFile = hairpinFSSclusteredFilename(resultFilePrefix, clusterSize, dropLogFactor, opt)
			load(resultFile)
		end
		makeFiniteSizeResultplot(clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt, plotopt, meanStartN)


		% collapse everything
		thisNs = Ns;
		thisClusterSize = numel(thisNs);
		graphFilePrefix = hairpinFSSclusteredFilenamePrefix(graphFilePrefixBase, thisClusterSize, dropLogFactor, T);
		resultFilePrefix = ['data/',graphFilePrefix];
		thisPlotopt = plotopt;
		thisPlotopt.graphFilePrefix = graphFilePrefix;
		[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = ...
		hairpinFSSclustered(thisNs, filenamePrefix, variableName, resultFilePrefix, thisClusterSize, dropLogFactor, opt, thisPlotopt, dataStartIndex)

		if minN < meanStartN
			% collapse starting from meanStartN
			thisNs = Ns(find(Ns >= meanStartN));
			thisClusterSize = numel(thisNs);
			graphFilePrefix = hairpinFSSclusteredFilenamePrefix(graphFilePrefixBase, thisClusterSize, dropLogFactor, T);
			resultFilePrefix = ['data/',graphFilePrefix];
			thisPlotopt = plotopt;
			thisPlotopt.graphFilePrefix = graphFilePrefix;
			[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = ...
			hairpinFSSclustered(thisNs, filenamePrefix, variableName, resultFilePrefix, thisClusterSize, dropLogFactor, opt, thisPlotopt, dataStartIndex)
		end
	end


	for squaredDeviation = sqDevWithSize
		% RESCALE WITH SIZE
		opt = baseOpt;
		opt.squaredDeviation = squaredDeviation
		opt.singleExponent = singleExponent;
		opt.rescaleWithSize = true;
		if generate
			[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = ...
			hairpinFSSclustered(Ns, filenamePrefix, variableName, resultFilePrefix, clusterSize, dropLogFactor, opt, plotopt, dataStartIndex)
		else
			resultFile = hairpinFSSclusteredFilename(resultFilePrefix, clusterSize, dropLogFactor, opt)
			load(resultFile)
		end
		makeFiniteSizeResultplot(clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt, plotopt, meanStartN)


		% collapse everything
		thisNs = Ns;
		thisClusterSize = numel(thisNs);
		graphFilePrefix = hairpinFSSclusteredFilenamePrefix(graphFilePrefixBase, thisClusterSize, dropLogFactor, T);
		resultFilePrefix = ['data/',graphFilePrefix];
		thisPlotopt = plotopt;
		thisPlotopt.graphFilePrefix = graphFilePrefix;
		[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = ...
		hairpinFSSclustered(thisNs, filenamePrefix, variableName, resultFilePrefix, thisClusterSize, dropLogFactor, opt, thisPlotopt, dataStartIndex)

		if minN < meanStartN
			% collapse starting from meanStartN
			thisNs = Ns(find(Ns >= meanStartN));
			thisClusterSize = numel(thisNs);
			graphFilePrefix = hairpinFSSclusteredFilenamePrefix(graphFilePrefixBase, thisClusterSize, dropLogFactor, T);
			resultFilePrefix = ['data/',graphFilePrefix];
			thisPlotopt = plotopt;
			thisPlotopt.graphFilePrefix = graphFilePrefix;
			[clustNs, clustPs, clustPErrs, clustQuals, clustQualErrs, opt] = ...
			hairpinFSSclustered(thisNs, filenamePrefix, variableName, resultFilePrefix, thisClusterSize, dropLogFactor, opt, thisPlotopt, dataStartIndex)
		end
	end
    end
end

