%filename  = 'hairpinStemScaling_combined' %destination file
filename  = 'hairpinStemScaling' %destination file

dir="/home/other/roald/clusterdata/hairpinFormation/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/relaxT90C_zipT10C_unzipT90C_allowUnb0.1_zippedRel100/dt15/";

%dir="/home/other/roald/clusterdata/oldHairpinFormation/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time10000*/";

%NO EXCLUSION, 180C
%dir="/home/other/roald/clusterdata/oldHairpinFormation/hairpinFormationNoExcl/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/";
%NO EXCLUSION, 100C
%dir="/home/other/roald/clusterdata/oldHairpinFormation/hairpinFormationNoExcl/formation_zipT20_unzipT100_allowUnb2_allowB2_zippedRel10/dt15_time100000/";
%Ns = [5 10 15 20 30 40 50 60 80 100];
Ns = [10 15 20 30 40 50 60 80];
fitstart = 1; %start at this index for fitting
%Ns = [15 20 30 40 50 60 80];
%Ns = [30 40 50 60 80];
zippingTimes = zeros(size(Ns));
zippingTimeErrs = zeros(size(Ns));
unzippingTimes = zeros(size(Ns));
unzippingTimeErrs = zeros(size(Ns));
for i = 1 : numel(Ns)
	[zipTime, zipErr, unzipTime, unzipErr, zipFromNuclTime, zipFromNuclErr] = averageZippingTime([dir, "N", num2str(Ns(i)), "*/formation*itf11"]);
	zippingTimes(i) = zipTime;
	zippingTimeErrs(i) = zipErr;
	unzippingTimes(i) = unzipTime;
	unzippingTimeErrs(i) = unzipErr;
	zippingFromNuclTimes(i) = zipFromNuclTime;
	zippingFromNuclTimeErrs(i) = zipFromNuclErr;
end

%[Ns; zippingTimes; unzippingTimes]


[zipCte, zipExponent, zipCteErr, zipExponentErr] = loglogRegression(Ns'(fitstart:end), zippingTimes'(fitstart:end), 1, 1, zippingTimeErrs'(fitstart:end));
zipC = [zipCte, zipCteErr]
zipE = [zipExponent, zipExponentErr]

[zipFromNuclCte, zipFromNuclExponent, zipFromNuclCteErr, zipFromNuclExponentErr] = loglogRegression(Ns'(fitstart:end), zippingFromNuclTimes'(fitstart:end), 1, 1, zippingFromNuclTimeErrs'(fitstart:end));
zipFromNuclC = [zipFromNuclCte, zipFromNuclCteErr]
zipFromNuclE = [zipFromNuclExponent, zipFromNuclExponentErr]

[unzipCte, unzipExponent, unzipCteErr, unzipExponentErr] = loglogRegression(Ns'(fitstart:end), unzippingTimes'(fitstart:end), 1, 1, unzippingTimeErrs'(fitstart:end));
unzipC = [unzipCte, unzipCteErr]
unzipE = [unzipExponent, unzipExponentErr]


fitNs = linspace(min(Ns), max(Ns), 100);
fitZipTimes = zipCte * fitNs.^zipExponent;
fitZipFromNuclTimes = zipFromNuclCte * fitNs.^zipFromNuclExponent;
fitUnzipTimes = unzipCte * fitNs.^unzipExponent;


clf;
hold on;

col='b'
%col='c'
h = loglogerr(Ns, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);
loglog(fitNs, fitZipTimes, col, "linewidth", 4);


col='g'
h = loglogerr(Ns, zippingFromNuclTimes, zippingFromNuclTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);
loglog(fitNs, fitZipFromNuclTimes, col, "linewidth", 4);


col='r'
%col='m'
h = loglogerr(Ns, unzippingTimes, unzippingTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);
loglog(fitNs, fitUnzipTimes, col, "linewidth", 4);

axis([4,100,0,1], 'autoy');

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

