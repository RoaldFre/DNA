addpath generic

generate = true;

variableName = 'totalEnergy';
filenamePrefixBase = 'hairpinTotalEnergy';
graphFilePrefixBase = 'energy';

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.resultFilePrefix = 'results/energy';
plotopt.measurement = 'the base pairing energy';
plotopt.quantity = '\Vbp';
plotopt.delta = 1;
delta = 1;

dropLogFactor = 3;

sqDevWithTime = [true];
sqDevWithSize = [true];

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, dropLogFactor, delta, sqDevWithTime, sqDevWithSize)
