Ts = [0 10 20 30 40 50 60 70 80 90];
%Ts = [10 30 50 70];
N = 6;

clf; hold on;

middle="TTAAAAAA"
plotZippingTemperature(middle, N, Ts, "b");

middle="TGAAAACA"
plotZippingTemperature(middle, N, Ts, "r");

hold off;


filename  = 'meltingTemperature';
caption   = 'Melting transition of two DNA hairpins: 3''CCTA\textbf{TT}AAAA\textbf{AA}TAGG5'' (blue) and 3''CCTA\textbf{TG}AAAA\textbf{CA}TAGG5'' (red). The simulation is performed at an ionic strenght of 115\,mM, compatible with \cite{vallone1999melting}. The data is fitted by a tanh function, giving a critical temperature of $(55.4 \pm 0.5)\degree$C and $(50.3 \pm 0.8)\degree$C for the red (TG/CA) and blue (TT/AA) curve, respectively.';

%filename  = 'meltingTemperatureKnotts';
%caption   = 'Melting transition of the same DNA hairpin sequences as in figure \ref{meltingTemperature}, but with base pairing interaction constants as per Knotts \etal \cite{knotts2007coarse}, while keeping the interaction btween matching hairpin base pairs only. This leads to very unphysical melting temperatures of ($8.7 \pm 0.3)\degree$C for the red (TG/CA) curve and ($-13 \pm 2)\degree$C for the blue (TT/AA) curve.';

destdir   = '../report/images';
relImgDir = 'images'; %relative to where your latex project root directory is
ylabrule  = '-1.5cm';
xlab      = 'Temperature ($^\circ$C)';
ylab      = 'Fraction of bound base pairs in the stem';
width     = '700';
height    = '500';

%makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);

presentationDir = '../presentation/images';
makeGraphPresentation(filename,presentationDir,xlab,ylab,ylabrule,width,height);
