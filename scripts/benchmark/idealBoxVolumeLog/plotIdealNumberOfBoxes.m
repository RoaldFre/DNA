addpath ../../generic

set(0, 'defaultlinelinewidth', 1.5);

Ns = [10 14 20 28 40 57 80 113 160 226 320];

fitStartN = 160;


destdir   = '../../../../thesis/latex/images/plots/'
relImgDir = 'images/plots'; %relative to where your latex project root directory is

filenameAllBoxes  = 'idealNumberOfBoxes'
filenameComplexity  = 'SPcomplexityFit'

spScaled = zeros(size(Ns));
spScaledErr = zeros(size(Ns));
noSpScaled = zeros(size(Ns));
noSpScaledErr = zeros(size(Ns));

numNs = numel(Ns);
colors = plotColors(numNs);
plotLegend = cell(1, numNs);
clf; hold on;
for i = 1 : numNs
	N = Ns(i);
	data = load(["ideal_",num2str(N)]);
	boxes = data(:,1);
	periods = data(:,2);
	worldsizes = data(:,3);
	boxVolumes = (worldsizes ./ boxes).^3;
	times = data(:,4:end);
	nSamples = size(times, 2);
	avgTimes = mean(times, 2);
	errTimes = std(times, 0, 2) / sqrt(nSamples - 1);

	% scaled: CPU seconds per simualted nanosecond per monomer
	scaled = avgTimes ./ periods / N;
	scaledErr = errTimes ./ periods / N;

	indices = find(mod(boxes, 2) == 1); % only uneven number of boxes...
	h = semilogyerror(boxes(indices), scaled(indices), scaledErr(indices))
	plotColorsI = round(size(colors, 1) * i / numNs);
	set(h ,'Color', colors(plotColorsI,:));
	plotLegend{i} = ["$N$ = ",num2str(N)];
	
	noSpScaled(i) = scaled(1);
	spScaled(i) = scaled(end);
	noSpScaledErr(i) = scaledErr(1);
	spScaledErr(i) = scaledErr(end);
end
hold off
legend(plotLegend);

%axis([1,100,0.1 3]) % TODO XXX HARDCODED AXES!
axis([0,100,0.1 3]) % TODO XXX HARDCODED AXES!

% TODO ASSUMED IN CAPTION: worldsize = 2000
caption   = ['CPU time per simulated nanosecond per monomer for various strand lengths $N$ in function of number of boxes per dimension. In all cases, the simulated world was a cube of length 200\,$\mu$m. For a single box and sufficiently large $N$, the CPU time per monomer increases linearly with the size $N$ (note that the $N$ sizes scale exponentially, and that the vertical axis is logarithmic). In the limit of $N \to \infty$, the CPU time per monomer reaches a constant for the maximum number of boxes (boxes with length equal to the potential truncation length). These two extreme cases are shown in Figure \ref{',filenameComplexity,'}. All simulations were run on an Intel Core i5-3570 CPU with the Langevin molecular dynamics integrator at a time step of 10\,fs.'];

ylabrule  = '-1.5cm';
xlab      = 'Number of boxes per dimension';
ylab      = 'CPU time (seconds) per simulated nanosecond per monomer';
width     = '900';
height    = '700';

set (gca, 'xminortick', 'on', 'yminortick', 'on');
makeGraph(filenameAllBoxes,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);








clf

%MAIN FIGURE
ylabrule  = '-1.5cm';
xlab      = 'Number of monomers $N$';
ylab      = 'CPU time (seconds) per simulated nanosecond';

fitIs = find(Ns >= fitStartN);
%fitNs = linspace(10,1000,1e3);
fitNs = linspace(0,400,1e3);

hold on

hNoSp = ploterror(Ns, noSpScaled.*Ns, noSpScaledErr.*Ns, 'r')
[cteNoSp, exponentNoSp, cteStddevNoSp, exponentStddevNoSp] = loglogRegression(Ns(fitIs), noSpScaled(fitIs).*Ns(fitIs), 1, 1, noSpScaledErr(fitIs).*Ns(fitIs))
plot(fitNs, cteNoSp*fitNs.^exponentNoSp, 'r')

hSp = ploterror(Ns, spScaled.*Ns, spScaledErr.*Ns)
[cteSp, exponentSp, cteStddevSp, exponentStddevSp] = loglogRegression(Ns(fitIs), spScaled(fitIs).*Ns(fitIs), 1, 1, spScaledErr(fitIs).*Ns(fitIs))
plot(fitNs, cteSp*fitNs.^exponentSp)
legend([hNoSp, hSp], {'without space partitioning', 'with space partitioning'});
%legend('location', 'northwest');
legend('location', 'north');

%axis([min(fitNs),max(fitNs),0.5 1e4]) % TODO XXX HARDCODED AXES!
axis([min(fitNs),max(fitNs),0, 800]) % TODO XXX HARDCODED AXES!

xlabel(xlab);
ylabel(['\rule{0pt}{',ylabrule,'}',ylab]);



% INSET
hax = axes();
set(hax, "position", [0.16 0.46 0.38 0.4]);
hold on;
fitNs = [9 1000];
hNoSp = loglogerror(Ns, noSpScaled.*Ns, noSpScaledErr.*Ns, 'r')
[cteNoSp, exponentNoSp, cteStddevNoSp, exponentStddevNoSp] = loglogRegression(Ns(fitIs), noSpScaled(fitIs).*Ns(fitIs), 1, 1, noSpScaledErr(fitIs).*Ns(fitIs))
loglog(fitNs, cteNoSp*fitNs.^exponentNoSp, 'r')

hSp = loglogerror(Ns, spScaled.*Ns, spScaledErr.*Ns)
[cteSp, exponentSp, cteStddevSp, exponentStddevSp] = loglogRegression(Ns(fitIs), spScaled(fitIs).*Ns(fitIs), 1, 1, spScaledErr(fitIs).*Ns(fitIs))
loglog(fitNs, cteSp*fitNs.^exponentSp)
%legend([hNoSp, hSp], {'without space partitioning', 'with space partitioning'});
%legend('location', 'northwest');
axis([min(fitNs),max(fitNs),1 1e3]) % TODO XXX HARDCODED AXES!

set(gca, 'xminortick', 'on', 'yminortick', 'on');

% TODO: looks like the inset gets no 'number labels' on the axes!

hold off




caption   = sprintf('CPU time per simulated nanosecond with and without space partitioning, in function of strand length $N$. In all cases, the simulated world was a cube of length 200\\,$\\mu$m, with 1 and 100 boxes per dimension for the curve without and with space partitioning, respectively (see Figure \\ref{%s} for the behaviour with intermediate number of boxes). The plotted curves are fitted power laws (fitted on the data points with $N \\ge %d$ to capture the asymptotic regime) with exponent $%s$ for the times without space partitioning, and $%s$ for the times with space partitioning. In the limit of $N \\to \\infty$, these exponents are expected to converge to two and one, respectively. The inset shows a log-log plot. All simulations were run on an Intel Core i5-3570 CPU with the Langevin molecular dynamics integrator at a time step of 10\\,fs.', filenameAllBoxes, fitStartN, numErr2tex(exponentNoSp, exponentStddevNoSp), numErr2tex(exponentSp, exponentStddevSp));

width     = '900';
height    = '700';

sleep(1e-9);
makeGraph(filenameComplexity,caption,destdir,relImgDir,'','',ylabrule,width,height);
