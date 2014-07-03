function makeHairpinGyradFSSplots(suff, descr)
if nargin < 1; suff = ''; end
if nargin < 2; descr = ''; end

addpath generic

generate = true;

variableName = 'dists';
filenamePrefixBase = ['hairpinEndToEnd',suff];
graphFilePrefixBase = ['endToEnd',suff];

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.resultFilePrefix = ['results/endToEnd',suff];
plotopt.measurement = ['the end-to-end distance ',descr];
plotopt.quantity = '\Rete';
plotopt.delta = '\nu';
delta = 0.587597;


sqDevWithTime = [true];
sqDevWithSize = [true];

%minN = 20;
minN = 0;

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, delta, sqDevWithTime, sqDevWithSize, minN)

