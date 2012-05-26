addpath("..");

destdir   = '../../report/images/';
relImgDir = 'images'; %relative to where your latex project root directory is
width     = '800';
height    = '600';
ylabrule  = '-1.5cm';

Ns = [100 200 300 400 500 600 700 800 900 1000];

bestB = zeros(size(Ns));

clf; hold on;
for i = 1 : numel(Ns)
	data = load(["ideal",num2str(Ns(i))]); % load the files ideal<NumOfMonomers>
	avgs = mean(data(:, 3:end)')';
	normalized = avgs ./ data(:, 2);
	bestIndex = find(avgs == min(avgs));
	bestB(i) = data(bestIndex(1), 1);
	bestTime(i) = min(normalized);
	singleBoxTime(i) = normalized(1);
	plot(data(:,1), normalized, "linewidth", 3);
end
hold off;
filename  = 'performancePerNumBox';
caption   = 'Performance in function of the number of boxes for various strand lengths $N$ of ssDNA. The bottom curve corresponds to $N = 100$, each successively higher curve represents an additional 100 monomers, up to the top curve, where $N = 1000$.';
xlab      = 'Number of boxes per dimension';
ylab      = 'CPU time in seconds to simulate one nanosecond.';
width     = '800';
height    = '800';
makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);

[Ns; bestB]

[cte, exponent, cteStddev, exponentStddev] = loglogRegression(Ns', bestB', 1, 1/3)

%fitNs = logspace(log10(min(Ns)), log10(max(Ns)), 50);
fitNs = logspace(2, 3, 50);
fit = cte * fitNs.^exponent;

figure;
loglog(Ns, bestB,"*", "linewidth", 4, fitNs, fit, "linewidth", 4);

filename  = 'numBoxExponentFit';
caption   = 'Ideal number of boxes per dimension in function of ssDNA strand length (logarithmic axes). The fitted curve has exponent $0.43 \pm 0.07$.';
xlab      = 'Number of monomers';
ylab      = 'Ideal number of boxes per dimension';
width     = '700';
height    = '500';
makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);



figure;
plot(Ns, bestTime, "*-", "linewidth", 4, Ns, singleBoxTime, "r*-", "linewidth", 4);

filename  = 'spacePartVsNoSpacePart';
caption   = 'Performance of spatial partitioning with ideal number of boxes (blue curve) compared to naive quadratic pair search (corresponding to only a single box for the entire world, red curve).';
xlab      = 'Number of monomers';
ylab      = 'CPU time in seconds per simulated nanosecond.';
width     = '700';
height    = '500';
makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
