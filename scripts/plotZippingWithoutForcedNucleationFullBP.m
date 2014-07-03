addpath "generic"

set(0, 'defaultlinelinewidth', 1.5);

destDir = '../../thesis/latex/images/plots/';
relImgDir = 'images/plots';


% HISTOGRAM OF:
% Log of unzipping times is pretty normally distributed
% Log of zipFromNucl times is pretty normally distributed
% however:
% log of zipping times is skewed to the right
% because:
% log of nucleation times is skewed to the right
%
% so it is natural to average the unzipping and zipFromNucl times in log 
% space(??)
% (likewise: when doing smallest squares in loglog, we are basically 'going 
% through the mean' in loglog space if the data is symmetrically 
% distributed -- in our case 'normally' distributed)
%
% Note: for the number of bound base pairs in time, we averaged in real 
% space before taking the loglog and doing linear regressions!



%if 0


dir="/home/other/roald/clusterdata/hairpinFormationLowTempFullBP_native_012/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/g5e11/dt10/zipT10C_unzipT90C_allowUnb0.2_nucl0.2/MCsweeps50000_relaxFactor0.3/"


NsTimes = [10 10;
           20 20;
%	   40 50;
	  ];
Ns = NsTimes(:,1);
times = NsTimes(:,2);

zipBoundBPsFraction = 0.7; % TODO TODO TODO DIFFERENT VALUE?
nucleationBoundBPs = 1;
NsFromNucl = Ns - nucleationBoundBPs;


numNs = numel(Ns);
numSamplesPerN = zeros(numNs, 1);
zippingTimes = cell(numNs, 1);
unzippingTimes = cell(numNs, 1);
zipFromNuclTimes = cell(numNs, 1);
nucleationTimes = cell(numNs, 1);
fullBound = cell(numNs, 1);
for i = 1 : numNs
	glob = [dir,'N',num2str(Ns(i)),'_minZipTime',num2str(times(i)),'/*itf11'];
	[zippingt, unzippingt, zipFromNuclt, nucleationt, preZipping, zipping, unzipping] = loadZippingTimeWithNucleationDataCustomThreshold(glob, nucleationBoundBPs, round(Ns(i) * zipBoundBPsFraction));

	zippingTimes{i} = zippingt;
	unzippingTimes{i} = unzippingt;
	zipFromNuclTimes{i} = zipFromNuclt;
	nucleationTimes{i} = nucleationt;

	numSamplesPerN(i) = numel(zippingt);

	fullBound{i} = cellfun(@(a,b) [a(:); b(:)], preZipping(:,2), zipping(:,2), 'UniformOutput', false);

	% Initial nucleation position histogram
	numRuns = numel(zippingt);
	averageNucleationState = zeros(1, Ns(i));
	for r = 1 : numRuns
		averageNucleationState += (zipping{r,3}(1,1:Ns(i)) == 1); % TODO Better way for 012 full BP?
	end
	averageNucleationState /= numRuns;
	subplot(3,1,i); % TODO HARDCODED FOR numNs == 3
	bar(averageNucleationState, 1)
	axis([0.5,Ns(i)+0.5,0,1],'autoy');
	title(["$S = ",num2str(Ns(i)),"$"]);
end
dt = preZipping{1,1}(3) - preZipping{1,1}(2);

xlab      = '';
ylab      = '';
ylabrule  = '-1.5cm';
width     = '900';
height    = '1200';
graphFileNuclPosition = 'zippingNucleationPositionHistogramsFullBP';

caption = "Histograms showing the fraction of successful nucleations per base. Low indices correspond to base pairs at the ends of the stem, large indices correspond to base pairs near the loop. For sufficiently long stem lengths, nearly all nucleation happens in the vicinity of the loop. All system sizes have around 200 independent runs, except $S = 57$, which has 64 runs.";

printf('\\input{images/plots/%s}\n', graphFileNuclPosition)

makeGraph(graphFileNuclPosition,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);










% Histogram of nucleation times
clf;
nBins = 20;
for i = 1 : numNs
	times = nucleationTimes{i};
	subplot(3,1,i); % TODO HARDCODED FOR numNs == 3
	%[NN, XX] = hist(log10(times), nBins);
	%bar(10.^XX, NN, 1);
	%hist(times, nBins);
	%set(gca,"xscale","log")
	normalizedHist(log10(times), nBins, -10, -5); % DANGEROUS! Hardcoded axes!
	axis([-10,-5,0,1]);% DANGEROUS! Hardcoded axes!
	title(["$S = ",num2str(Ns(i)),"$"]);
end


xlab      = '';
ylab      = '';
ylabrule  = '-1.5cm';
width     = '900';
height    = '1200';
graphFileNuclTime = 'zippingNucleationTimeHistogramsFullBP';

caption = "Histograms of nucleation times for different stem lengths $S$. Horizontal axes give the logarithm in base ten of the time till nucleation. Vertical axes are rescaled to give the histogram a surface area of unity. All system sizes have around 200 independent runs, except $S = 57$, which has 64 runs. We remark that these simulations were performed with $\\gamma = 5 \\cdot 10^{12}\\,$s$^{-1}$. Approximate comparison with other results in this thesis is possible after dividing the times reported here by a factor of ten.";

printf('\\input{images/plots/%s}\n', graphFileNuclTime);

makeGraph(graphFileNuclTime,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);





% Histogram of zipping-from-nucleation times
clf;
nBins = 40;
for i = 1 : numNs
	times = zipFromNuclTimes{i};
	subplot(3,1,i); % TODO HARDCODED FOR numNs == 3
	%[NN, XX] = hist(log10(times), nBins);
	%bar(10.^XX, NN, 1);
	%hist(times, nBins);
	%set(gca,"xscale","log")
	normalizedHist(log10(times), nBins, -10, -5); % DANGEROUS! Hardcoded axes!
	%axis([-10,-5,0,1]);% DANGEROUS! Hardcoded axes!
	title(["$S = ",num2str(Ns(i)),"$"]);
end


xlab      = '';
ylab      = '';
ylabrule  = '-1.5cm';
width     = '900';
height    = '1200';
graphFileFromNuclTime = 'zippingFromNucleationTimeHistogramsFullBP';

caption = ["Histograms of zippering times for different stem lengths $S$. The zipper time is defined as the time from successful nucleation till the time where 80\\% of the base pairs are bound (changing this threshold did not affect the shape of the histograms significantly). Horizontal axes give the logarithm in base ten of the zippering time and are chosen identical to Figure \\ref{",graphFileNuclTime,"} for easy comparison. Vertical axes are rescaled to give the histogram a surface area of unity. All system sizes have around 200 independent runs, except $S = 57$, which has 64 runs. We remark that these simulations were performed with $\\gamma = 5 \\cdot 10^{12}\\,$s$^{-1}$. Approximate comparison with other results in this thesis is possible after dividing the times reported here by a factor of ten."];

printf('\\input{images/plots/%s}\n', graphFileFromNuclTime);

makeGraph(graphFileFromNuclTime,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);













% Fully averaged bound ('convolved with wait-for-nucleation')
clf; hold on;

toPlot = (Ns != 57); % Too noisy!
plotNs = Ns(toPlot);
numPlots = numel(plotNs);
plotIndices = find(toPlot);

numToAverageForEquilibriumZipped = 10;
colors = plotColors(numPlots);
plotLegend = cell(numPlots, 1);
plotHandles = zeros(numPlots, 1);
for i = 1:numPlots
	N = plotNs(i);
	numSamples = 0;
	equilibriumZipped = 0;
	fullBounds = fullBound{plotIndices(i)};
	numRuns = numel(fullBounds);
	for r = 1:numRuns
		numSamples = max(numSamples, numel(fullBounds{r}));
		equilibriumZipped += mean(fullBounds{r}(end-numToAverageForEquilibriumZipped : end));
	end
	equilibriumZipped /= numRuns;
	meanTotalBound = zeros(numSamples, 1);
	for r = 1:numRuns
		b = fullBounds{r};
		meanTotalBound += [b; equilibriumZipped * ones(numSamples - numel(b), 1)];
	end
	meanTotalBound /= numRuns;
	t = (1:numSamples)*dt;
	h = plot(t/1e-6, meanTotalBound);
	set(h ,'Color', colors(i,:));
	plotLegend{i} = ['$S$ = ',num2str(N)];
	plotHandles(i) = h;
end
hold off;
legend(plotHandles, plotLegend, 'location', 'northwest');

axis([0,8,0,1],'autoy')

xlab      = 'time $t$ ($\mu$s)';
ylab      = 'mean number of bound base pairs $\expect{\nb(t)}$';
ylabrule  = '-1.5cm';
width     = '900';
height    = '800';
graphFileZippingFromQuench = 'zippingFromQuenchFullBP';

caption = "Mean number of bound base pairs during zipping after temperature quench from $T = 90\\degree$C to $T = 10\\degree$C. Each curve is the result of around 200 independent runs. Convergence to the mean is very slow due to the long hunt for nucleation, followed by rapid zippering. This is especially clear in the curve for $S = 40$, where we can observe the contributions of the zippering phases of individual runs. We remark that these simulations were performed with $\\gamma = 5 \\cdot 10^{12}\\,$s$^{-1}$. Approximate comparison with other results in this thesis is possible after dividing the times reported here by a factor of ten.";

printf('\\input{images/plots/%s}\n', graphFileZippingFromQuench);

makeGraph(graphFileZippingFromQuench,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);








%end






% number of bound base pairs of several individual runs for N = 40
clf
%dataSetIndex = find(Ns == 40); % TODO TODO
dataSetIndex = find(Ns == 20); % TODO TODO
dataSet = fullBound{dataSetIndex};
runIndices = [1 2 3 4 5 6 7 8 9 10];
runIndices = [6 8 10 80 100 110 120 130 140];
%runIndices = [1 8 10 5];
numPlots = numel(runIndices)
for i = 1:numPlots
	run = runIndices(i);
	subplot(numPlots,1,i);
	bound = dataSet{run};
	t = (1:numel(bound))*dt;
	h = plot(t/1e-6, bound);
	%axis([0,8,0,40]);

	% Quick hack to write data files that were requested by promotor
	%saveData = [t(:), bound(:)];
	%save('-ascii', ['data/individualRun',num2str(run)], 'saveData');
end



xlab      = 'time $t$ ($\mu$s)';
ylab      = '';
ylabrule  = '-1.5cm';
width     = '900';
height    = '1200';
graphFileZippingFromQuenchIndividual = 'zippingFromQuenchIndividualRunsFullBP';

caption = [""];
%caption = ["Number of bound base pairs during zipping after temperature quench from $T = 90\\degree$C to $T = 10\\degree$C of several individual simulations with $S = 40$. The vertical axis represents the number of bound base pairs $\\nb(t)$. The horizontal axes are identical to Figure \\ref{",graphFileZippingFromQuench,"} for easy comparison. The hunt for successful nucleation, followed by rapid zippering is clearly visible. We remark that these simulations were performed with $\\gamma = 5 \\cdot 10^{12}\\,$s$^{-1}$. Approximate comparison with other results in this thesis is possible after dividing the times reported here by a factor of ten."];

printf('\\input{images/plots/%s}\n', graphFileZippingFromQuenchIndividual);

makeGraph(graphFileZippingFromQuenchIndividual,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);



