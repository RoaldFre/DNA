function makeHairpinStateFSSplots(suff, descr)
if nargin < 1; suff = ''; end
if nargin < 2; descr = ''; end

addpath generic

generate = true;

variableName = 'bound';
filenamePrefixBase = ['hairpinState',suff];
graphFilePrefixBase = ['state',suff];

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.resultFilePrefix = ['results/state',suff];
plotopt.measurement = ['the number of bound base pairs ',descr];
plotopt.quantity = '\nb';
plotopt.delta = 1;
delta = 1;

sqDevWithTime = [false];
sqDevWithSize = [false];

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, delta, sqDevWithTime, sqDevWithSize)

