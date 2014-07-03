function bpAndRg(N, T, dropLogFactor, plotopt)
addpath generic

if nargin < 1; N = 20; end;
if nargin < 2; T = 10; end;
if nargin < 3; dropLogFactor = 200; end;

set(0, 'defaultlinelinewidth', 1.5);

filenamePrefix = 'data/hairpin';
filename = [filenamePrefix,'StateT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
load(filename);
nbTime = time;
[nb, nbErr] = meanOfSamples(bound);

filename = [filenamePrefix,'GyradT',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
load(filename);
gyrTime = time;
[gyr, gyrErr] = meanOfSamples(gyrad);



zippedFraction = nb / mean(nb(end-3:end));
armsFraction = 1-zippedFraction;

%factor = 5; % TODO hardcoded in caption!
gyrExagerated = gyr/max(gyr);
factor = 1/(1 - min(gyrExagerated))
gyrExagerated = 1 - factor*(1 - gyrExagerated);

clf; hold on
plot(gyrTime/1e-6, gyr/max(gyr), 'r');
plot(gyrTime/1e-6, gyrExagerated, 'r--');
plot(nbTime/1e-6, nb/max(nb), 'b');
legend({'\smaller $\expect{\Rg}$','\smaller $\expect{\Rg}$ (stretched)','\smaller $\expect{\nb}$'}, 'location', 'east');
hold off



xlab = 'time $t$ ($\mu$s)'
ylab = 'arbitrary units'
ylabrule  = '-1.5cm';
width     = '800';
height    = '700';
caption = ['Mean gyration radius (red) and number of bound base pairs (blue) during zipping of a hairpin with $S = ',num2str(N),'$ at $T = ',num2str(T),'\degree$C. Notice that the gyration radius reaches its minimum at approximately the point where half of the base pairs are zipped. An additional $\Rg$ curve that is stretched to fill the full plot window is shown, in order to help in visually locating the minimum.'];

axis([0, 0.03, 0, 1]); %TODO HARDCODED FOR S=20

graphFile = 'nbAndRg';
destDir = '../../thesis/latex/images/plots/';
relImgDir = 'images/plots';
printf('\\input{images/plots/%s}\n', graphFile)

set(gca, 'xminortick', 'on', 'yminortick', 'on');
makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);
