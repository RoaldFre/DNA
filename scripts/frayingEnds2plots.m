function frayingEnds2plots(N1, N2, T, dropLogFactor, plotopt)
addpath generic

if nargin < 1; N1 = 20; end;
if nargin < 2; N2 = 40; end;
if nargin < 3; T = 10; end;
if nargin < 4; dropLogFactor = 200; end;

Ns = [N1 N2]

factor = 1.7;

set(0, 'defaultlinelinewidth', 1.5);

clf
for i = 1:2
	N = Ns(i);
	subplot(1,2,i)
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


	hold on
	loglog(nbTime, nb, 'b')
	loglog(nbNoDihTime, nbNoDih, 'r')
	%legend({'\smaller dihedral','\smaller no dihedral'}, 'location', 'northwest');
	legend({'\smaller with $\Vdih$','\smaller without $\Vdih$'}, 'location', 'northwest');
	hold off

	xlabel('time $t$ ($\mu$s)');
	if i == 1
		ylabel('mean number of bound base pairs $\expect{\nb(t)}$');
	end

	if N == 20
		axis([2e-11, 1e-7, 0.1, 50]);
	elseif N == 40
		axis([2e-11, 1e-7, 0.1, 50]);
	end
	set(gca, 'xminortick', 'on', 'yminortick', 'on');
end



xlab = '';
ylab = '';
ylabrule  = '-1.5cm';
width     = '1200';
height    = '700';
caption = ['Effect of the dihedral interaction on the mean number of bound base pairs during zipping of a hairpin with $S = ',num2str(N1),'$ (left) and $S = ',num2str(N2),'$ (right) at $T = ',num2str(T),'\degree$C. The base pairing interaction was increased with a factor ',num2str(factor),' for the curve without dihedral interaction to attain a similar equilibrium zipped value $\expect{\nb(\infty)}$ as for the original model.'];


graphFile = ['frayingEndsS',num2str(N1),'andS',num2str(N2)];
destDir = '../../thesis/latex/images/plots/';
relImgDir = 'images/plots';
printf('\\input{images/plots/%s}\n', graphFile)

makeGraph(graphFile,caption,destDir,relImgDir,xlab,ylab,ylabrule,width,height);
