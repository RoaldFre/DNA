% Makes plots of the raw hairpin data (optionally the squaredDevaiton of 
% the raw data).
function makeHairpinPlots(Ns, T, filenamePrefix, variableName, squaredDeviation, plotopt, withErrorBars, dataStartIndex)
addpath generic
more off

if nargin < 7; withErrorBars = false; end
if nargin < 8; dataStartIndex = 1; end

if isempty(plotopt)
	error "plotopt is empty!"
end

dropLogFactor = 200;
set(0, 'defaultlinelinewidth', 1.5);

destDir = plotopt.destDir;
relImgDir = plotopt.relImgDir;
graphFilePrefix = plotopt.graphFilePrefix;
measurement = plotopt.measurement;
quantity = plotopt.quantity;
units = plotopt.units;
unitsSq = plotopt.unitsSq;

numNs = numel(Ns);

clf;hold on;
colors = plotColors(numNs) ;
plotLegend = cell(1, numNs);
plotHandles = zeros(1, numNs);
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

	color = colors(i, :);

	if squaredDeviation
		[ys, dys, xs] = squaredMeanDeviation(data, time, false);
		sqDevPostfix = "_sqDev";
	else
		xs = time;
		[ys, dys] = meanOfSamples(data);
		sqDevPostfix = "";
	end

	h = loglog(xs, ys);
	set(h ,'Color', color);
	plotHandles(i) = h;
	if withErrorBars
		dysFrac = 1 + dys ./ ys;
		h = loglog(xs, ys .* dysFrac, 'r');
		set(h ,'Color', color);
		h = loglog(xs, ys ./ dysFrac, 'r');
		set(h ,'Color', color);
	end

	plotLegend{i} = ["$S$ = ",num2str(N)];

	% Save ascii file
	saveData = [xs(:), ys(:), dys(:)];
	saveFile = [filename,sqDevPostfix,'_ascii'];
	save('-ascii', saveFile, 'saveData');
end
hold off;
%legend(plotHandles, plotLegend, 'location', 'northwest');
legend(plotHandles, plotLegend, 'location', 'eastoutside');

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

graphFile = [graphFilePrefix,'T',num2str(T),'N',NsStrFile,sqDevSuffix];

if withErrorBars
	errorBarsStr = 'Error bars are plotted as well.';
else
	errorBarsStr = '';
end

%caption = sprintf('Simulation result of %s %s $%s$ at temperature $T = %d^\\circ$C for stem lengths $S$ = %s. %s', sqDevStr, measurement, quantity, T, NsStr, errorBarsStr);
%Stem lenghts are already in legend! Don't repeat ourselves or cause clutter!
caption = sprintf('Simulation result of %s %s $%s$ at temperature $T = %d^\\circ$C. %s', sqDevStr, measurement, quantity, T, errorBarsStr);
xlab = '$t$ (s)';
if isempty(units)
	ylab = ['$',observable,'$'];
else
	if squaredDeviation
		ylab = ['$',observable,'$ (',unitsSq,')'];
	else
		ylab = ['$',observable,'$ (',units,')'];
	end
end

ylabrule  = '-1.5cm';
%width     = '900';
width     = '1100';
height    = '900';

printf('\\input{images/plots/%s}\n', graphFile)

set(gca, 'xminortick', 'on', 'yminortick', 'on');
makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);
