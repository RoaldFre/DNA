% Plot the end-to-end distance based on the number of bound base pairs 
% under the assumption that the arms would have time to equilibrate.
function equilibriumEndToEnd(N, T, dropLogFactor, plotopt)
addpath generic

if nargin < 1; N = 48; end;
if nargin < 2; T = 10; end;
if nargin < 3; dropLogFactor = 200; end;

set(0, 'defaultlinelinewidth', 1.5);

filenamePrefix = 'data/hairpin';
filename = [filenamePrefix,'StateT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
load(filename);
nbTime = time;
[nb, nbErr] = meanOfSamples(bound);

filename = [filenamePrefix,'EndToEndT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
load(filename);
eteTime = time;
[ete, eteErr] = meanOfSamples(dists);



zippedFraction = nb / mean(nb(end-3:end));
armsFraction = 1-zippedFraction;

equilibriumEte = ete(end) + (ete(1) - ete(end)) .* armsFraction.^(0.588);

clf; hold on
plot(eteTime/1e-6, ete/max(ete), 'r');
plot(nbTime/1e-6, equilibriumEte/max(equilibriumEte), 'g')
plot(nbTime/1e-6, nb/max(nb), 'b');
legend({'\smaller $\expect{\Rete}$','\smaller $\Rete^{\text{eq}}(\expect{\nb})$','\smaller $\expect{\nb}$'}, 'location', 'east');
hold off



xlab = 'time $t$ ($\mu$s)'
ylab = 'arbitrary units'
ylabrule  = '-1.5cm';
width     = '800';
height    = '700';
caption = ['Mean end-to-end distance (red) and number of bound base pairs (blue) during zipping of a hairpin with $S = ',num2str(N),'$ at $T = ',num2str(T),'\degree$C. Notice that the end-to-end distance only starts to respond to the zipping when a significant fraction of the monomers have already zipped. For comparison, the green $\Rete^\text{eq}(\expect{\nb})$ curve corresponds to how the end-to-end distance would behave if it could equilibrate while it zipped (with the zipping rate still determined by the measured $\expect{\nb}$). This curve is generated as $(\expect{\nb(\infty)} - \expect{\nb(t)})^\nu$ with the amplitude and added offset from the actual $\Rete$ measurement.'];

axis([0, 0.08, 0, 1]);

graphFile = 'endToEndAndNumbound';
destDir = '../../thesis/latex/images/plots/';
relImgDir = 'images/plots';
printf('\\input{images/plots/%s}\n', graphFile)

set(gca, 'xminortick', 'on', 'yminortick', 'on');
makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);

