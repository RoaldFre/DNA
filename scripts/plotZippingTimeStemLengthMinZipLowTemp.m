filename  = 'hairpinStemScaling' %destination file

addpath "generic"



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

parseData = true


if parseData

dir="/home/other/roald/clusterdata/hairpinFormationLowTemp_native/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/dt15/zipT10C_unzipT90C_allowUnb0.2_nucl0.2/MCsweeps100000_relaxFactor2/"


NsTimes = [10 100;
           14 100;
           20 200;
           28 300;
	   40 500;
%	   57 800;
	  ];
Ns = NsTimes(:,1);
times = NsTimes(:,2);

zipBoundBPsFraction = 0.5; % TODO TODO TODO DIFFERENT VALUE?
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
		averageNucleationState += zipping{r,3}(1,1:Ns(i));
	end
	averageNucleationState /= numRuns;
	bar(averageNucleationState, 1)
	title(["Nucleation histogram for S=",num2str(Ns(i))]);
	xlabel("base pair of stem (end of stem is left, loop is right)")
	ylabel("fraction of nucleations")
	print(["nucleationHistogram_N",num2str(Ns(i)),".eps"]);
	%sleep(3);
end

dt = preZipping{1,1}(3) - preZipping{1,1}(2);


%combine all samples
totalNumSamples = sum(numSamplesPerN);
allNs = zeros(totalNumSamples, 1);
allNsFromNucl = zeros(totalNumSamples, 1);
allZippingTimes = zeros(totalNumSamples, 1);
allUnzippingTimes = zeros(totalNumSamples, 1);
allZipFromNuclTimes = zeros(totalNumSamples, 1);
allNucleationTimes = zeros(totalNumSamples, 1);
zippingTimeMeans = zeros(numNs, 1);
zippingTimeErrs  = zeros(numNs, 1);
zippingTimeMeans = zeros(numNs, 1);
unzippingTimeErrs  = zeros(numNs, 1);
unzippingTimeMeans = zeros(numNs, 1);
zipFromNuclTimeErrs  = zeros(numNs, 1);
zipFromNuclTimeMeans = zeros(numNs, 1);
nucleationTimeErrs  = zeros(numNs, 1);
nucleationTimeMeans = zeros(numNs, 1);
i = 1;
for j = 1 : numNs
	allZippingTimes    (i : i+numSamplesPerN(j)-1) = zippingTimes{j};
	allUnzippingTimes  (i : i+numSamplesPerN(j)-1) = unzippingTimes{j};
	allZipFromNuclTimes(i : i+numSamplesPerN(j)-1) = zipFromNuclTimes{j};
	allNucleationTimes (i : i+numSamplesPerN(j)-1) = nucleationTimes{j};
	allNs              (i : i+numSamplesPerN(j)-1) = Ns(j);
	allNsFromNucl      (i : i+numSamplesPerN(j)-1) = NsFromNucl(j);
	i += numSamplesPerN(j);

	%zippingTimeMeans(j) = mean(zippingTimes{j});
	zippingTimeMeans(j) = exp(mean(log(zippingTimes{j})));
	zippingTimeErrs(j)  = std(zippingTimes{j}) / sqrt(numSamplesPerN(j) - 1);

	%unzippingTimeMeans(j) = mean(unzippingTimes{j});
	unzippingTimeMeans(j) = exp(mean(log(unzippingTimes{j})));
	unzippingTimeErrs(j)  = std(unzippingTimes{j}) / sqrt(numSamplesPerN(j) - 1);

	%zipFromNuclTimeMeans(j) = mean(zipFromNuclTimes{j});
	zipFromNuclTimeMeans(j) = exp(mean(log(zipFromNuclTimes{j})));
	zipFromNuclTimeErrs(j)  = std(zipFromNuclTimes{j}) / sqrt(numSamplesPerN(j) - 1);

	%nucleationTimeMeans(j) = mean(nucleationTimes{j});
	nucleationTimeMeans(j) = exp(mean(log(nucleationTimes{j})));
	nucleationTimeErrs(j)  = std(nucleationTimes{j}) / sqrt(numSamplesPerN(j) - 1);
end

%[Ns; zippingTimes; unzippingTimes]

end







[zipCte, zipExponent, zipCteErr, zipExponentErr] = loglogRegressionBootstrap(allNs, allZippingTimes, 1, 1);
zipC = [zipCte, zipCteErr]
zipE = [zipExponent, zipExponentErr]

[zipFromNuclCte, zipFromNuclExponent, zipFromNuclCteErr, zipFromNuclExponentErr] = loglogRegressionBootstrap(allNsFromNucl, allZipFromNuclTimes, 1, 1);
zipFromNuclC = [zipFromNuclCte, zipFromNuclCteErr]
zipFromNuclE = [zipFromNuclExponent, zipFromNuclExponentErr]

[unzipCte, unzipExponent, unzipCteErr, unzipExponentErr] = loglogRegressionBootstrap(allNs, allUnzippingTimes, 1, 1);
unzipC = [unzipCte, unzipCteErr]
unzipE = [unzipExponent, unzipExponentErr]

[nuclCte, nuclExponent, nuclCteErr, nuclExponentErr] = loglogRegressionBootstrap(allNs, allNucleationTimes, 1, 1);
nuclC = [nuclCte, nuclCteErr]
nuclE = [nuclExponent, nuclExponentErr]


factor = (max(Ns) / min(Ns))^(1/10);
fitNs = linspace(min(Ns) / factor, max(Ns) * factor, 100);
fitNsFromNucl = linspace(min(NsFromNucl) / factor, max(NsFromNucl) * factor, 100);
fitZipTimes = zipCte * fitNs.^zipExponent;
fitZipFromNuclTimes = zipFromNuclCte * fitNsFromNucl.^zipFromNuclExponent;
fitUnzipTimes = unzipCte * fitNs.^unzipExponent;
fitNuclTimes = nuclCte * fitNs.^nuclExponent;




clf;
hold on;

col='b'
%col='c'
h = loglogerr(Ns, zippingTimeMeans, zippingTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);
loglog(fitNs, fitZipTimes, col, "linewidth", 4);


col='g'
h = loglogerr(NsFromNucl, zipFromNuclTimeMeans, zipFromNuclTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);
loglog(fitNsFromNucl, fitZipFromNuclTimes, col, "linewidth", 4);


col='r'
%col='m'
h = loglogerr(Ns, unzippingTimeMeans, unzippingTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);
loglog(fitNs, fitUnzipTimes, col, "linewidth", 4);


col='m'
h = loglogerr(Ns, nucleationTimeMeans, nucleationTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);
loglog(fitNs, fitNuclTimes, col, "linewidth", 4);



%axis([4,100,0,1], 'autoy');

hold off;


%caption   = 'Scaling of the hairpin zipping (blue) and unzipping (red) times for an (A)$_N$(T)$_N$ hairpin with the length $N$ varying between 5 and 100. The scaling exponent is fitted, with exclusion of the first two data points (at the short lenghts of $N = 5$ and $N = 10$), yielding a scaling behaviour $\tau \sim N^{1.33 \pm 0.14}$ for the zipping time and a scaling of $\tau \sim N^{2.57 \pm 0.10}$ for the unzipping time.';

destdir   = '../report/images';
relImgDir = 'images'; %relative to where your latex project root directory is
ylabrule  = '2.0cm';
xlab      = 'Number of base pairs, $N$';
ylab      = 'Hairpin (un)zipping time (seconds)'
width     = '700';
height    = '500';

%makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);

%presentationDir = '../presentation/images';
%makeGraphPresentation(filename,presentationDir,xlab,ylab,ylabrule,width,height);

end

