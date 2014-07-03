% Generic plotter for force-vel-fric stuff.

function plotHairpinTension(fileNameBase, variableName, legendName, T, N, dropLogFactor, dataStartIndex, blockAvg, particleType, verticalLabel, direction, observable, yfactor)

set(0, 'defaultlinelinewidth', 1.5);

dataFile = [fileNameBase,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
load(dataFile);

if not(exist(variableName))
	error(["Can't find the variable with name '",variableName,"' in the data file '",dataFile,"'!"]);
end

if not(isempty(legendName))
	legendName = [legendName, ' '];
end



theAxis = [1e-10, 1.3e-7, 0, 1]; % TODO HARDCODED (use autoy)


eval(['data = ',variableName,';'])
nRuns = size(data, 1)
%nSamples = size(data, 2)
%nCols = size(data, 3)

plotMean = squeeze(mean(data,1));
plotErr  = squeeze(std(data,0,1)) / sqrt(nRuns - 1);
plotErrFrac = 1 + abs(plotErr./ plotMean);

time = time(dataStartIndex:end);
plotMean = plotMean(dataStartIndex:end, :);
plotErr = plotErr(dataStartIndex:end, :);
plotErrFrac = plotErrFrac(dataStartIndex:end, :);

clf;hold on;

subplot(2,1,1);
hold on;

%%enkel 1 kant -- kan symmetriseren met gespiegelde monomeren!
%plotCols = 1:(N+6);
plotCols = 1:(N+5); % TODO HARDCODED FOR SUGARS (they start at index 2)
%plotCols = 1:(2*N+12);
numPlots = numel(plotCols);
colors = plotColors(numPlots);
plotLegend = cell(1, numPlots);
%styles = {'s-', 'd-', '^-', 'v-', 'o-', 'p-'};
styles = {'s-', 's--', 'd-', 'd--', '^-', '^--', 'v-', 'v--', 'o-', 'o--', 'p-', 'p--'};
% - for solid, -. for dash-dot, : for dots, -- for dotted.
for i = 1:numPlots
	col = plotCols(i);
	style = styles{mod(i-1,numel(styles))+1};
	h = semilogx(blockAverage(time, blockAvg), blockAverage(plotMean(:,col)*yfactor, blockAvg), style);
	set(h ,'Color', colors(i,:));
	plotLegend{i} = [legendName,num2str(col + 1)]; % TODO COL + 1 HARDCODED!
end
hold off;
axis(theAxis, 'autoy');
set(gca, 'xminortick', 'on', 'yminortick', 'on');
%legend(plotLegend, 'location', 'eastoutside');
legend(plotLegend, 'location', 'east');
%title(['Variable "',variableName,'" from data set with T=',num2str(T),'C and dropLog=',num2str(dropLogFactor),', averaging ',num2str(blockAvg),' samples'])
%plotFile = [dataFile,'_plot_',variableName,'_blockAvg',num2str(blockAvg),'.eps'];
%print(plotFile, '-depsc2');

ylabel(verticalLabel);
xlabel("time $t$ (s)");

sleep(1e-9);





subplot(2,1,2);

%%enkel 1 kant -- kan symmetriseren met gespiegelde monomeren!
%plotCols = 1:(N+6);
%plotCols = (N+7):(2*N+12);
%plotCols = (N+7):(2*N+11); % ONE LESS BECAUSE THE FvP and FvS DATASETS LACK ONE MONOMER!
% XXX 
%plotCols = (N+7) : size(plotMean,2); % Always right :-)
plotCols = (N+6) : size(plotMean,2); % TODO HARDCODED FOR SUGARS

numPlots = numel(plotCols);
colors = plotColors(numPlots);
plotLegend = cell(1, numPlots);
%styles = {'s-', 'd-', '^-', 'v-', 'o-', 'p-'};
styles = {'s-', 's--', 'd-', 'd--', '^-', '^--', 'v-', 'v--', 'o-', 'o--', 'p-', 'p--'};
% - for solid, -. for dash-dot, : for dots, -- for dotted.
%clf;hold on;
hold on;
for i = 1:numPlots
	col = plotCols(i);
	style = styles{mod(i-1,numel(styles))+1};
	h = semilogx(blockAverage(time, blockAvg), blockAverage(plotMean(:,col)*yfactor, blockAvg), style);
	set(h ,'Color', colors(i,:));
	plotLegend{i} = [legendName,num2str(col + 1)]; % TODO COL + 1 HARDCODED!
end
hold off;
axis(theAxis, 'autoy');
set(gca, 'xminortick', 'on', 'yminortick', 'on');
%legend(plotLegend, 'location', 'eastoutside');
legend(plotLegend, 'location', 'east');
%title(['Variable "',variableName,'" from data set with T=',num2str(T),'C and dropLog=',num2str(dropLogFactor),', averaging ',num2str(blockAvg),' samples'])
%plotFile = [dataFile,'_plot_',variableName,'_blockAvg',num2str(blockAvg),'_mirror.eps'];
%print(plotFile, '-depsc2');

ylabel(verticalLabel);
xlabel("time $t$ (s)");

sleep(1e-9);








xlab = '';
ylab = '';
ylabrule  = '-1.5cm';
width     = '1250';
height    = '1200';
caption = ['Results for the ',observable,' the ',particleType,' in the ',direction,' direction for hairpins of stem length $S = ',num2str(N),'$ at a temperature $T = ',num2str(T),'\degree$C. Adjacent samples have been averaged to decrease the visual noise. The top plot shows the results for one side of the hairpin, the bottom plot shows the results for the the other halve. Particle indices are given as labels.'];


tmp = ['/',dataFile];
dataBaseName = tmp(find(tmp == '/', 1, 'last')+1 : end);
graphFile = [dataBaseName,'_plot_',variableName,'_blockAvg',num2str(blockAvg)];
destDir = '../../thesis/latex/images/plots/';
relImgDir = 'images/plots';
printf('\\input{images/plots/%s}\n', graphFile)

makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);
