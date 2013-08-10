addpath generic

generate = true;

variableName = 'dists';
filenamePrefixBase = 'hairpinEndToEnd';
graphFilePrefixBase = 'endToEnd';

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.resultFilePrefix = 'results/endToEnd';
plotopt.measurement = 'the end-to-end distance';
plotopt.quantity = '\Rete';
plotopt.delta = '\nu';
delta = 0.587597;

dropLogFactor = 3;

%sqDevWithTime = [true];
sqDevWithTime = [];
sqDevWithSize = [true];

%minN = 20;
minN = 0;

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, dropLogFactor, delta, sqDevWithTime, sqDevWithSize, minN)

