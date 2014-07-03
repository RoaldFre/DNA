function plotRelaxation(T, gamma, dt, sweeps, dropLogFactor, DLdataStartIndex, N, totTime, wait1, wait2)
% WARNING, THIS OVERWRITES THE CACHED FILE!
% WARNING, THIS OVERWRITES THE CACHED FILE!
% WARNING, THIS OVERWRITES THE CACHED FILE!
% WARNING, THIS OVERWRITES THE CACHED FILE!
% WARNING, THIS OVERWRITES THE CACHED FILE!


loadFromCache = false;

clf
plotHairpinStateLowTempGamma(N, T, totTime, sweeps, wait1, gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, 'r');
sleep(1e-9);
plotHairpinStateLowTempGamma(N, T, totTime, sweeps, wait2, gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, 'b');
sleep(1e-9);

%% RESET THE CACHE TO THE FULL DATA SET
%plotHairpinStateLowTempGamma(N, T, totTime, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, 'k');
