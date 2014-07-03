addpath "generic"

loadFromCache = true;

Ns = [7 10 14 20 28 40 57];
T = 10;

set(0, 'defaultlinelinewidth', 1.5);
clf
numNs = numel(Ns);
colors = plotColors(numNs);
plotLegend = cell(numNs, 1);
plotHandles = zeros(numNs, 1);
for i = 1:numNs
	N = Ns(i);
	h = plotHairpinState(N, T, loadFromCache)
	set(h ,'Color', colors(i,:));
	plotHandles(i) = h;
	plotLegend{i} = ['$S = ',num2str(N),'$'];
end
legend(plotHandles, plotLegend, 'location', 'northwest');



xlab = 'time $t$ (s)'
ylab = 'mean number of bound base pairs $\expect{\nb(t)}$'
ylabrule  = '-1.5cm';
width     = '900';
height    = '800';
caption = 'Mean number of bound base pairs when quenching strands with forced nucleation centers from $T = 90\degree$C to $T = 10\degree$C. The four XY base pairs of the nucleation center are not taken into account to determine $\nb(t)$. We remark that these simulations were performed with $\gamma = 5 \times 10^{12}\,$s$^{-1}$. Approximate comparison with other results in this thesis is possible after dividing the times reported here by a factor of ten.';

graphFile = 'zippingWithNucleationAndQuench';
destDir = '../../thesis/latex/images/plots/';
relImgDir = 'images/plots';
printf('\\input{images/plots/%s}\n', graphFile)

set(gca, 'xminortick', 'on', 'yminortick', 'on');
makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);
