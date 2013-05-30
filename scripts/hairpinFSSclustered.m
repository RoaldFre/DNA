% Ns: list of system sizes
% filenamePrefix:   Prefix of the data files. The full datafile must be of the form
%                   <filenamePrefix>'N'<N>, with <N> the system size as per Ns, e.g.: 'dataN10'
% variableName:     The name of the variable that holds the data in the data file
%                   each row in that matrix is assumed to be an independent run
%                   there is also assumed to be a 'time' variable with corresponding time values!
% clusterSize:      Do FSS on clusters of <clusterSize> adjacent sytem sizes
% dropLogFactor;    Transform linear time samples to log space time samples with dropDataLog() with this factor
% singleExponent:   If 'true', then do single exponent fit based on nu and guessBeta
% fixOffsetToZero:  Bool: fit without or with extra x,y offsets
% nu:               Used if singleExponent == true, then alpha = nu/beta
% scalingFunction:  Function handle to rescale data. e.g. 'finiteSizeRescaleWithTime' or 'finiteSizeRescaleWithSize'
% squaredDeviation: Take squared deviation of the data samples?
function [clustNs, alphas, alphaErrs, betas, betaErrs] = hairpinFSSclustered(Ns, filenamePrefix, variableName, clusterSize, dropLogFactor, singleExponent, fixOffsetToZero, nu, scalingFunction, squaredDeviation, guessBeta, guessAlpha, eps, bootstrapSamples, dataStartIndex)

addpath generic

more off

if nargin < 11; guessBeta = 1; end
if nargin < 12; guessAlpha = 1; end
if nargin < 13; eps = 1e-5; end
if nargin < 14; bootstrapSamples = 30; end
if nargin < 15; dataStartIndex = 1; end


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
	end

	if not(exist(variableName))
		error ["Can't find the variable with name '",variableName,"' in the data file '",filename,"'!"]
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
clustFits = zeros(numClusters, 2);
clustFitErrs = zeros(numClusters, 2);
for i = 1 : numClusters
	clustNs = Ns(i : i+clusterSize-1);
	clustXsDecim = timesDecim(i : i+clusterSize-1);
	clustYsDecim = datasDecim(i : i+clusterSize-1);

	if singleExponent
		exponentGuessOrNu = nu;
	else
		exponentGuessOrNu = guessAlpha;
	end
	[alpha, alphaErr, beta, betaErr, quality, xOffsets, xOffsetsErr, yOffsets, yOffsetsErr] = ...
		finiteSizeScaling(clustNs, clustXsDecim, clustYsDecim, [], exponentGuessOrNu, guessBeta, scalingFunction, eps, bootstrapSamples, fixOffsetToZero, singleExponent, squaredDeviation);

	clustQuals(i) = quality;
	clustFits(i,:) = [alpha; beta];
	clustFitErrs(i,:) = [alphaErr, betaErr];
	
	if not(fixOffsetToZero)
		clustFitxOfsets(:,i) = xOffsets;
		clustFityOfsets(:,i) = yOffsets;
		clustFitxOfsetErrs(:,i) = xOffsetsErr;
		clustFityOfsetErrs(:,i) = yOffsetsErr;
	end

	%clustNs
	%[xOffsets yOffsets]
end


clustNs = zeros(numClusters, 1);
for i = 1 : numClusters
	clustNs(i) = mean(Ns(i : i+clusterSize-1), 'g'); % geometric mean
end

clf; hold on
%betas:
betas = clustFits(:,2);
betaErrs = clustFitErrs(:,2);
ploterror(clustNs, betas, betaErrs)
%alphas:
alphas = clustFits(:,1);
alphaErrs = clustFitErrs(:,1);
ploterror(clustNs, alphas, alphaErrs, 'r')
hold off

