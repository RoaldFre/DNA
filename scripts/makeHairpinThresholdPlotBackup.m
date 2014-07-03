% Get exponents directly from the times when quantities (for different 
% system sizes) cross a certain threshold.
%
%  - thresholdParams: Each column corresponds to a set of parameters (one 
%                     row for each system size in Ns) that gets passed to 
%                     the findTauIndex function.
%                     The first column is used as the parameters for the 
%                     plot on the raw data. The other columns are used for 
%                     error estimate on the plot without the raw data.
%  - findTauIndex: Function: ind = findTauIndex(time, ys, thresholdParam, N) 
%                  that determines the index of the characteristic time 
%                  'tau' in the given (time,ys) data.
function makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs, xaxis, yaxis)

addpath generic
more off

if nargin < 13; withErrorBars = false; end
if nargin < 14; dataStartIndex = 1; end
if nargin < 15; plotTauYs = true; end
if nargin < 16; xaxis = []; end
if nargin < 17; yaxis = []; end

if isempty(plotopt)
	error "plotopt is empty!"
end

set(0, 'defaultlinelinewidth', 1.5);

bootstrapSamples = 100;
fitStartIndex = find(Ns >= fitStartN, 1, 'first');


destDir = plotopt.destDir;
relImgDir = plotopt.relImgDir;
graphFilePrefix = plotopt.graphFilePrefix;
measurement = plotopt.measurement;
quantity = plotopt.quantity;
units = plotopt.units;
unitsSq = plotopt.unitsSq;

numNs = numel(Ns);


numParamSets = size(thresholdParams, 2);
taus = zeros(numNs, bootstrapSamples + 1, numParamSets);
tauYs = zeros(numNs, bootstrapSamples + 1, numParamSets);
tauErrs = zeros(numNs, numParamSets);
tauYErrs = zeros(numNs, numParamSets);

clf
subplot(2,1,1);
hold on;
for i = 1:numNs
	N = Ns(i);
	filename = [filenamePrefix,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
	load(filename);

	if not(exist(variableName))
		error(["Can't find the variable with name '",variableName,"' in the data file '",filename,"'!"]);
	end

	eval(['data = ',variableName,';'])
	data = data(:,dataStartIndex:end);
	time = time(dataStartIndex:end);

	% Get the non-bootstrapped curve and tau data first
	if squaredDeviation
		[ys, dys, xs] = squaredMeanDeviation(data, time, false);
		sqDevPostfix = "_sqDev";
	else
		xs = time;
		[ys, dys] = meanOfSamples(data);
		sqDevPostfix = "";
	end
	% get tau data (for all parameters)
	for p = 1:numParamSets
		ind = findTauIndex(xs, ys, thresholdParams(i, p), N);
		if isempty(ind)
			error(['Could not find a tau index for threshold parameter ',num2str(thresholdParams(i,p)),' -- Threshold parameter too extreme?']);
		end
		taus(i, 1, p) = xs(ind);
		tauYs(i, 1, p) = ys(ind);
	end

	color = [1 1 1]*0.7;
	h = loglog(xs, ys, '-', 'linewidth', 1.1);
	set(h ,'Color', color);
	if withErrorBars
		dysFrac = 1 + dys ./ ys;
		h = loglog(xs, ys .* dysFrac, '-');
		set(h ,'Color', color);
		h = loglog(xs, ys ./ dysFrac, '-');
		set(h ,'Color', color);
	end



	% Do the boostrapping
	for b = 1:bootstrapSamples
		if squaredDeviation
			[ys, dys, xs] = bootstrapSampleSquaredDeviation(data, time);
		else
			xs = time;
			[ys, dys] = bootstrapSampleMean(data);
		end

		for p = 1:numParamSets
			ind = findTauIndex(xs, ys, thresholdParams(i, p), N);
			if isempty(ind)
				error(['Could not find a tau index for threshold parameter ',num2str(thresholdParams(i,p)),' for a bootstrap resample -- Threshold parameter too extreme or insufficient statistics?']);
			end
			taus(i, 1+b, p) = xs(ind);
			tauYs(i, 1+b, p) = ys(ind);
		end
	end
	for p = 1:numParamSets
		tauErrs(:, p) = std(squeeze(taus(:,:,p)), 0, 2);
		tauYErrs(:, p) = std(squeeze(tauYs(:,:,p)), 0, 2);
	end
end



% Determine the exponents from the tau data (tau in function of N)
exponents = zeros(1, numParamSets);
ctes = zeros(1, numParamSets);
exponentErrs = zeros(numParamSets, 1);
cteErrs = zeros(numParamSets, 1);
for p = 1:numParamSets
	% non bootstrapped result
	[cte, exponent] = loglogRegression(Ns(fitStartIndex:end), taus(fitStartIndex:end,1,p), 1, 1);
	ctes(p) = cte;
	exponents(p) = exponent;
	
	% use bootstrapped results for error
	cs = zeros(bootstrapSamples, 1);
	es = zeros(bootstrapSamples, 1);
	for b = 1:bootstrapSamples
		[cte, exponent] = loglogRegression(Ns(fitStartIndex:end), taus(fitStartIndex:end,b+1,p), 1, 1);
		cs(b) = cte;
		es(b) = exponent;
	end
	%exponentErrs
	%es
	%std(es)
	exponentErrs(p) = std(es);
	cteErrs(p) = std(cs);

end


% Determine the exponent of the threshold points (tauY in function of tau)
tauYctes = zeros(numParamSets,1);
tauYcteErrs = zeros(numParamSets,1);
tauYexponents = zeros(numParamSets,1);
tauYexponentErrs = zeros(numParamSets,1);
for p = 1:numParamSets
	[tauYcte, tauYexponent, tauYcteErr, tauYexponentErr] = loglogRegression(taus(fitStartIndex:end,:,p), tauYs(fitStartIndex:end, :, p), 1, 1);
	tauYctes(p) = tauYcte;
	tauYexponents(p) = tauYexponent;
	tauYcteErrs(p) = tauYcteErr;
	tauYexponentErrs(p) = tauYexponentErr;
end

% Adjust axes if needed
if not(isempty(xaxis))
	if not(isempty(yaxis))
		axis([xaxis(:); yaxis(:)]);
	else
		axis(xaxis);
	end
else
	if not(isempty(yaxis))
		axis([0; 1; yaxis(:)], 'autox');
	end
end







% Plot the (tau, tauY) points and fitted line for each parameter
tauYfit = [tauYexponents(:), tauYctes(:)]
plotHandles = zeros(1, numParamSets);
if plotTauYs
	theAxis = axis()
	plotXs = theAxis(1:2);
	colors = plotColors(numParamSets);
	for p = 1:numParamSets
		color = colors(p,:);
		h = loglog(taus(:,p,1), tauErrs(:,p), '.', 'linewidth', 3);
		set(h ,'Color', color);
		h = loglog(plotXs, tauYctes(p)*plotXs.^tauYexponents(p), '-');
		set(h ,'Color', color);
		set(h ,'linewidth', 1.8);
		plotHandles(p) = h;
	end
	axis(theAxis); % in case plotting the lines changed the axes.

	if not(isempty(thresholdLegendStrings))
		legend(plotHandles, thresholdLegendStrings, 'location', 'eastoutside')
	end
end

% End of plotting first graph
%hold off;
set(gca, 'xminortick', 'on', 'yminortick', 'on');





if squaredDeviation
	sqDevSuffix = 'sqDev';
	sqDevStr = 'the square of the mean deviation of ';
	observable = sprintf('\\sqdif{%s(t)}', quantity);
else
	sqDevSuffix = '';
	sqDevStr = '';
	observable = sprintf('\\expect{%s(t)}', quantity);
end

NsStr = '';
NsStrFile = '';
for N = Ns
	NsStr = [NsStr,num2str(N),', '];
	NsStrFile = [NsStrFile,num2str(N),'.'];
end
NsStr = NsStr(1:end-2);
NsStrFile = NsStrFile(1:end-1);

if withErrorBars
	errorBarsStr = 'Error bars are plotted as well.';
else
	errorBarsStr = '';
end

xlabel('$t$ (s)');
if isempty(units)
	ylabel(['$',observable,'$']);
else
	if squaredDeviation
		ylabel(['$',observable,'$ (',unitsSq,')']);
	else
		ylabel(['$',observable,'$ (',units,')']);
	end
end











% Make a second graph for the different threshold parameter values with error bars this time

subplot(2,1,2);
hold on
plotHandles = zeros(1, numParamSets);


colors = plotColors(numParamSets);
% First plot the points with error bars, to determine ideal axes
for p = 1:numParamSets
	color = colors(p,:);
	h = loglogerror(Ns, taus(:,1,p), tauErrs(:,p));
	set(h ,'Color', color);
end
theAxis = axis()
plotXs = theAxis(1:2);
% Then plot the fitted lines
for p = 1:numParamSets
	h = loglog(plotXs, ctes(p) * plotXs.^exponents(p), '--');
	set(h ,'Color', color);
	set(h ,'linewidth', 2);
	plotHandles(p) = h;
end
% Reset the axes, which may have changed by plotting the lines
axis(theAxis);

xlabel('$S$')
ylabel('$\tau$')
% End of plotting second graph
set(gca, 'xminortick', 'on', 'yminortick', 'on');
if not(isempty(thresholdLegendStrings))
	%legend(plotHandles, thresholdLegendStrings, 'location', 'northwest');
	legend(plotHandles, thresholdLegendStrings, 'location', 'eastoutside');
end
%hold off;

allExponents = exponents
[exponents(1) exponentErrs(1)]





% Make a table with results





caption = 'The fits in the top plot are of the form $';


ylabrule  = '-1.5cm';
xlab      = '';
ylab      = '';
width     = '1000';
height    = '1100';

graphFile = [graphFilePrefix,'_threshold_T',num2str(T),'N',NsStrFile,sqDevSuffix];

printf('\\input{images/plots/%s}\n', graphFile)

makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height,'TODO');

