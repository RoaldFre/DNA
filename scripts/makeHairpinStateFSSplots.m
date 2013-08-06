addpath generic

generate = true;

variableName = 'bound';
filenamePrefixBase = 'hairpinState';
graphFilePrefixBase = 'state';

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.resultFilePrefix = 'results/state';
plotopt.measurement = 'the number of bound base pairs';
plotopt.quantity = '\nb';
plotopt.delta = 1;
delta = 1;

dropLogFactor = 3;

sqDevWithTime = [false];
sqDevWithSize = [false];

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, dropLogFactor, delta, sqDevWithTime, sqDevWithSize)

