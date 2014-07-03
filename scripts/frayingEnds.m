function frayingEnds(N, T, dropLogFactor, plotopt)
addpath generic

if nargin < 1; N = 20; end;
if nargin < 2; T = 10; end;
if nargin < 3; dropLogFactor = 200; end;

factor = 1.7;

set(0, 'defaultlinelinewidth', 1.5);

filenamePrefix = 'data/hairpin';
filename = [filenamePrefix,'StateT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
load(filename);
nbTime = time;
[nb, nbErr] = meanOfSamples(bound);


filenamePrefix = 'data/hairpin';
filename = [filenamePrefix,'StateNoDih',num2str(factor),'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
load(filename);
nbNoDihTime = time;
[nbNoDih, nbNoDihErr] = meanOfSamples(bound);


clf; hold on
loglog(nbTime, nb, 'b')
loglog(nbNoDihTime, nbNoDih, 'r')
legend({'with dihedral interaction','without dihedral interaction'}, 'location', 'northwest');
hold off



xlab = 'time $t$ ($\mu$s)'
ylab = 'mean number of bound base pairs $\expect{\nb(t)}$'
ylabrule  = '-1.5cm';
width     = '800';
height    = '700';
caption = ['Effect of the dihedral interaction on the mean number of bound base pairs during zipping of a hairpin with $S = ',num2str(N),'$ at $T = ',num2str(T),'\degree$C. The base pairing interaction was increased with a factor ',num2str(factor),' for the curve without dihedral interaction to attain a similar equilibrium zipped value $\expect{\nb(\infty)}$ as for the original model.'];

if N == 20
	axis([2e-11, 4e-8, 0.07, 20]);
elseif N == 40
	axis([2e-11, 1e-7, 0.1, 40]);
end

graphFile = ['frayingEndsS',num2str(N)];
destDir = '../../thesis/latex/images/plots/';
relImgDir = 'images/plots';
printf('\\input{images/plots/%s}\n', graphFile)

set(gca, 'xminortick', 'on', 'yminortick', 'on');
makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);
