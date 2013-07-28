addpath generic

generate = true;

variableName = 'gyrad';
filenamePrefixBase = 'hairpinGyrad';
graphFilePrefixBase = 'gyrad';

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.measurement = 'the gyration radius';
plotopt.quantity = '\Rg';
plotopt.delta = '\nu';
delta = 0.587597;

dropLogFactor = 2;

%sqDevWithTime = [true];
sqDevWithTime = [];
sqDevWithSize = [true, false];

minN = 20;

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, dropLogFactor, delta, sqDevWithTime, sqDevWithSize, minN)


