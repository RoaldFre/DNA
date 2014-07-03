% Get exponents directly from the times when quantities (for different 
% system sizes) cross a certain threshold.
%
%  - thresholdParams: Vector of parameters. Each entry gets gets passed to 
%                     the findTauIndex function once.
%  - findTauIndex: Function: ind = findTauIndex(time, ys, thresholdParam, N) 
%                  that determines the index of the characteristic time 
%                  'tau' in the given (time,ys) data.
function makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYfits, xaxis, yaxis, insetPosition, legendLocation)

addpath generic
more off

if nargin < 13; withErrorBars = false; end
if nargin < 14; dataStartIndex = 1; end
if nargin < 15; plotTauYfits = true; end
if nargin < 16; xaxis = []; end
if nargin < 17; yaxis = []; end
%if nargin < 18; insetPosition = [0.55 0.10 0.38 0.3]; end
%if nargin < 19; legendLocation = 'northwest'; end
if nargin < 18; insetPosition = [0.18 0.62 0.35 0.3]; end
if nargin < 19; legendLocation = 'southeast'; end

if isempty(plotopt)
	error "plotopt is empty!"
end

set(0, 'defaultlinelinewidth', 1.5);

bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100; % TODO TODO TODO TODO            TODO     TODO
bootstrapSamples = 100;
bootstrapSamples = 10;
fitStartIndex = find(Ns >= fitStartN, 1, 'first');


destDir = plotopt.destDir;
relImgDir = plotopt.relImgDir;
graphFilePrefix = plotopt.graphFilePrefix;
measurement = plotopt.measurement;
quantity = plotopt.quantity;
units = plotopt.units;
unitsSq = plotopt.unitsSq;

numNs = numel(Ns);


numParamSets = numel(thresholdParams);
taus = zeros(numNs, bootstrapSamples + 1, numParamSets);
tauYs = zeros(numNs, bootstrapSamples + 1, numParamSets);
tauErrs = zeros(numNs, numParamSets);
tauYErrs = zeros(numNs, numParamSets);

% Store to plot at the end on top all other stuff
allXs = cell(numNs, 1);
allYs = cell(numNs, 1);

clf
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
		ind = findTauIndex(xs, ys, thresholdParams(p), N);
		if isempty(ind)
			error(['Could not find a tau index for threshold parameter ',num2str(thresholdParams(p)),' -- Threshold parameter too extreme?']);
		end
		taus(i, 1, p) = xs(ind);
		tauYs(i, 1, p) = ys(ind);
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
			ind = findTauIndex(xs, ys, thresholdParams(p), N);
			if isempty(ind)
				error(['Could not find a tau index for threshold parameter ',num2str(thresholdParams(p)),' for a bootstrap resample -- Threshold parameter too extreme or insufficient statistics?']);
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











% Make an inset graph for the different threshold parameter values with error bars this time
hax = axes();
set(hax, "position", insetPosition);

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

% Force x axis to something more sensible
%axis([10,60,0,1],'autoy');
axis([10,60,1e-10,2e-7]);
theAxis = axis();

% LABELS (AND TICKLABELS) DON'T WORK WITH INSET >_<
%xlabel('$S$')
%ylabel('$\tau$')
% End of plotting second graph
set(gca, 'xminortick', 'on', 'yminortick', 'on');
%if not(isempty(thresholdLegendStrings))
%	%legend(plotHandles, thresholdLegendStrings, 'location', 'northwest');
%	legend(plotHandles, thresholdLegendStrings, 'location', 'eastoutside');
%end
%hold off;


table = [exponents(:) exponentErrs(:)]

% Make a table with results
if numParamSets == 1
	resultsString = ['The fitted values are $\kappa = ',numErr2tex(tauYexponents(1), tauYexponentErrs(1)),'$ and $\beta = ',numErr2tex(exponents(1), exponentErrs(1)),'$'];
else
	% Make a table
	formatString = 'c|';
	header = '$\lambda$';
	kappas = '$\kappa$';
	betas = '$\beta$';
	for p=1:numParamSets
		formatString = [formatString, 'r@{$\pm$}l'];
		header = [header,'& \multicolumn{2}{c}{$',num2tex(thresholdParams(p)),'$}'];
		kappas = [kappas,'& $',numErr2tex(tauYexponents(p), tauYexponentErrs(p), '$&$'),'$']; 
		betas = [betas,'& $',numErr2tex(exponents(p), exponentErrs(p), '$&$'),'$']; 
	end
	header = [header, '\\\hline'];
	kappas = [kappas, '\\'];
	if not(plotTauYfits)
		kappas = '';
	end
	%betas = [betas, '\\'];
	resultsString = ['The results for the various values of $\lambda$ are\\[0.5em]',...
			'{\centering\begin{tabular}{',formatString,'}',...
			header,kappas,betas,'\end{tabular}}.'];
end

if numParamSets == 1
	mainInfoStr = 'The main plot shows the position of the threshold values that define the zipping times $\tau$. ';
else
	mainInfoStr = 'The main plot shows the position of the threshold values that define the zipping times $\tau$ for several threshold fractions $\lambda$. ';
end
insetInfoStr = ['The horizontal axis of the log-log inset plot shows the stem length $S$ and runs from ',num2str(theAxis(1)),' to ',num2str(theAxis(2)),', the vertical axis corresponds to the extracted zipping times $\tau$ (the horizontal position of the points on the main plot) and runs from $',num2tex(theAxis(3)),'$ to $',num2tex(theAxis(4)),'$. '];

if plotTauYfits
	if numParamSets == 1
		mainFitString = ['The fit in the main plot has exponent $\kappa$ with $\expect{',quantity,'(\tau)} \sim \tau^\kappa$. '];
	else
		mainFitString = ['The fits in the main plot have exponent $\kappa$ with $\expect{',quantity,'(\tau)} \sim \tau^\kappa$. '];
	end
else
	mainFitString = '';
end
if numParamSets == 1
	caption = ['Analysis of the zipping time $\tau$ as determined from simulations of ',sqDevStr,measurement,' $',quantity,'$ at $T=',num2str(T),'\degree$C. ',mainInfoStr,insetInfoStr,mainFitString,'The fit in the inset is of the form $\tau \sim S^{\beta}$. Curves for stem lengths $S = ',NsStr,'$ are plotted, but only data points with $S \ge ',num2str(fitStartN),'$ are used for the fits. ',resultsString];
else
	caption = ['Analysis of the zipping time $\tau$ as determined from simulations of ',sqDevStr,measurement,' $',quantity,'$ at $T=',num2str(T),'\degree$C. ',mainInfoStr,insetInfoStr,mainFitString,'The fits in the inset are of the form $\tau \sim S^{\beta}$. Curves for stem lengths $S = ',NsStr,'$ are plotted, but only data points with $S \ge ',num2str(fitStartN),'$ are used for the fits.  ',resultsString];
end




















hax = axes();
set(hax, "position", [0 0 1 1]);




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
		ind = findTauIndex(xs, ys, thresholdParams(p), N);
		if isempty(ind)
			error(['Could not find a tau index for threshold parameter ',num2str(thresholdParams(p)),' -- Threshold parameter too extreme?']);
		end
		taus(i, 1, p) = xs(ind);
		tauYs(i, 1, p) = ys(ind);
	end


	%allXs{i} = xs;
	%allYs{i} = ys;
	color = [1 1 1]*0.6;
	h = loglog(xs, ys, '-', 'linewidth', 1.3);
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
			ind = findTauIndex(xs, ys, thresholdParams(p), N);
			if isempty(ind)
				error(['Could not find a tau index for threshold parameter ',num2str(thresholdParams(p)),' for a bootstrap resample -- Threshold parameter too extreme or insufficient statistics?']);
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
theAxis = axis();
tauYfit = [tauYexponents(:), tauYctes(:)]
plotHandles = zeros(1, numParamSets);
plotXs = theAxis(1:2);
colors = plotColors(numParamSets);
for p = 1:numParamSets
	color = colors(p,:);
	h = loglog(taus(:,1,p), tauYs(:,1,p), 'o', 'linewidth', 10, 'markersize', 1);
	set(h ,'Color', color);
	if plotTauYfits
		h = loglog(plotXs, tauYctes(p)*plotXs.^tauYexponents(p), '--');
		set(h ,'Color', color);
		set(h ,'linewidth', 6);
	end
	plotHandles(p) = h;
end
axis(theAxis); % in case plotting the lines changed the axes.

if not(isempty(thresholdLegendStrings))
	legend(plotHandles, thresholdLegendStrings, 'location', legendLocation)
end

%% Now plot the raw data on top of everything else
%for i = 1:numNs
%	color = [1 1 1]*0.6;
%	h = loglog(allXs{i}, allYs{i}, '-', 'linewidth', 1.3);
%	set(h ,'Color', color);
%	if withErrorBars
%		error "blah"
%	end
%end
%axis(theAxis);


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






























ylabrule  = '-1.5cm';
xlab      = '';
ylab      = '';
width     = '1000';
height    = '1000';

graphFile = [graphFilePrefix,'_threshold_T',num2str(T),'N',NsStrFile,sqDevSuffix];

printf('\\input{images/plots/%s}\n', graphFile)

makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height,'TODO');

