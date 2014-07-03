function makeHairpinGyradFSSplots(suff, descr)
if nargin < 1; suff = ''; end
if nargin < 2; descr = ''; end

addpath generic

generate = true;

variableName = 'totalEnergy';
filenamePrefixBase = ['hairpinTotalEnergy',suff];
graphFilePrefixBase = ['energy',suff];

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.resultFilePrefix = ['results/energy',suff];
plotopt.measurement = ['the base pairing energy ',descr];
plotopt.quantity = '\Vbp';
plotopt.delta = 1;
delta = 1;

sqDevWithTime = [false,true];
sqDevWithSize = [false,true];

makeHairpinFSSplots(variableName, filenamePrefixBase, graphFilePrefixBase, plotopt, delta, sqDevWithTime, sqDevWithSize)
