%filename  = 'hairpinLoopScaling_combined' %destination file
filename  = 'hairpinLoopScaling' %destination file

dir="/home/other/roald/clusterdata/loopFormation/formation_S=CG_L=A_zipT20_unzipT180_allowUnb0_allowB0_zippedRel10/dt15_time100000/";

#Ls = [2 4 6 8 10 14 18 22 26 30 40 50 60 70 80 90 100];
Ls = [2 4 6 8 10 14    22       34 38 40      60];
fitstart = 1; %start at this index for fitting
zippingTimes = zeros(size(Ls));
zippingTimeErrs = zeros(size(Ls));
unzippingTimes = zeros(size(Ls));
unzippingTimeErrs = zeros(size(Ls));
for i = 1 : numel(Ls)
	[zipTime, zipErr, unzipTime, unzipErr] = averageZippingTime([dir, "L", num2str(Ls(i)), "*/formation*"]);
	zippingTimes(i) = zipTime;
	zippingTimeErrs(i) = zipErr;
	unzippingTimes(i) = unzipTime;
	unzippingTimeErrs(i) = unzipErr;
end

%[Ls; zippingTimes; unzippingTimes]


[zipCte, zipExponent, zipCteErr, zipExponentErr] = loglogRegression(Ls'(fitstart:end), zippingTimes'(fitstart:end), 1, 1, zippingTimeErrs'(fitstart:end));
zipC = [zipCte, zipCteErr]
zipE = [zipExponent, zipExponentErr]

[unzipCte, unzipExponent, unzipCteErr, unzipExponentErr] = loglogRegression(Ls'(fitstart:end), unzippingTimes'(fitstart:end), 1, 1, unzippingTimeErrs'(fitstart:end));
unzipC = [unzipCte, unzipCteErr]
unzipE = [unzipExponent, unzipExponentErr]


fitLs = linspace(min(Ls), max(Ls), 100);
fitZipTimes = zipCte * fitLs.^zipExponent;
fitUnzipTimes = unzipCte * fitLs.^unzipExponent;


clf;
hold on;

col='b';
%col='c';
h = loglogerr(Ls, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitLs, fitZipTimes, col, "linewidth", 4);

col='r'
%col='m'
h = loglogerr(Ls, unzippingTimes, unzippingTimeErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitLs, fitUnzipTimes, col, "linewidth", 4);

axis([4,100,0,1], 'autoy');

hold off;




return



caption   = 'Scaling of the hairpin zipping (blue) and unzipping (red) times for an (A)$_N$(T)$_N$ hairpin with the length $N$ varying between 5 and 100. The scaling exponent is fitted, with exclusion of the first two data points (at the short lenghts of $N = 5$ and $N = 10$), yielding a scaling behaviour $\tau \sim N^{1.33 \pm 0.14}$ for the zipping time and a scaling of $\tau \sim N^{2.57 \pm 0.10}$ for the unzipping time.';

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

