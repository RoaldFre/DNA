
factor = 1.7

gamma = '5e11';
dt = 10;
sweeps = 20e3;
%dropLogFactor = 3;
dropLogFactor = 5;
%dropLogFactor = 10;
dropLogFactor = 200;
DLdataStartIndex = 2;
color = 'b';
basedir = ['/home/other/roald/clusterdata_properVelocityInit_bpRepulsiveCORRECT/hairpinState_lowTemp_gamma_noDih/baseFact',num2str(factor),'/'];
suff = ['NoDih',num2str(factor)];

loadFromCache = false;

T = 10;
%plotHairpinStateLowTempGamma(5, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(6, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(7, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(8, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(10, T, 20, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(12, T, 20, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(14, T, 20, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(17, T, 25, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
plotHairpinStateLowTempGamma(20, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
plotHairpinStateLowTempGamma(24, T, 45, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
plotHairpinStateLowTempGamma(28, T, 60, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
plotHairpinStateLowTempGamma(34, T, 75, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
plotHairpinStateLowTempGamma(40, T, 85, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(48, T, 100, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(57, T, 120, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(67, T, 160, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(80, T, 200, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);


T = 30;
%plotHairpinStateLowTempGamma(5, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(6, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);   
%plotHairpinStateLowTempGamma(7, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(8, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);   
%plotHairpinStateLowTempGamma(10, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(12, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);  
%plotHairpinStateLowTempGamma(14, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(17, T, 40, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);  
%plotHairpinStateLowTempGamma(20, T, 50, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(24, T, 60, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);  
%plotHairpinStateLowTempGamma(28, T, 75, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(34, T, 100, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);  
%plotHairpinStateLowTempGamma(40, T, 120, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(48, T, 160, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff); 
%plotHairpinStateLowTempGamma(57, T, 200, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
%plotHairpinStateLowTempGamma(67, T, 280, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff); 
%plotHairpinStateLowTempGamma(80, T, 350, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
