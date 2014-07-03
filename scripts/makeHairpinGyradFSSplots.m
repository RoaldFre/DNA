function makeHairpinGyradFSSplots(suff, descr)
if nargin < 1; suff = ''; end
if nargin < 2; descr = ''; end

addpath generic

generate = true;

variableName = 'gyrad';
filenamePrefixBase = ['hairpinGyrad',suff];
graphFilePrefixBase = ['gyrad',suff];

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.resultFilePrefix = ['results/gyrad',suff];
plotopt.measurement = ['the gyration radius ',descr];
plotopt.quantity = '\Rg';
plotopt.delta = '\nu';
delta = 0.587597;

sqDevWithTime = [true];
sqDevWithSize = [true, false];

%minN = 20;
minN = 0;

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, delta, sqDevWithTime, sqDevWithSize, minN)


