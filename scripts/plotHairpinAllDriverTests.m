gamma = '5e11';
dt = 10;
sweeps = 20e3;

dropLogFactor = 200
%dropLogFactor = 4
DLdataStartIndex = 2;
#loadFromCache = false;
loadFromCache = true;

#funcs = {@plotHairpinStateLowTempGamma, @plotHairpinEnergyLowTempGamma, @plotHairpinGyradLowTempGamma, @plotEndToEndZipping};
#funcs = {@plotHairpinGyradLowTempGamma};
funcs = {@plotEndToEndZipping};

if true
for f = funcs
	f = f{1};
	T = 10;
	
	%f(40, T, 70, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(34, T, 60, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(28, T, 45, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(24, T, 40, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(20, T, 30, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);

	f(22, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);

	continue

	%%f(80, T, 200, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%%f(67, T, 160, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%%f(57, T, 120, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(48, T, 100, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(44, T, 85, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(40, T, 70, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(37, T, 65, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(34, T, 60, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(31, T, 55, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(28, T, 45, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(26, T, 43, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(24, T, 40, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(22, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(20, T, 30, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(17, T, 25, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(14, T, 20, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(12, T, 20, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(10, T, 20, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(8, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(7, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(6, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(5, T, 10, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);

	T = 30;
	%%f(80, T, 350, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%%f(67, T, 280, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache); 
	%%f(57, T, 200, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(48, T, 160, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache); 
	f(44, T, 140, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(40, T, 120, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(37, T, 110, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(34, T, 100, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);  
	f(31, T, 85, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);  
	%f(28, T, 75, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	f(26, T, 70, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);  
	%f(24, T, 60, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);  
	f(22, T, 55, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(20, T, 50, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(17, T, 40, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);  
	%f(14, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(12, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);  
	%f(10, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(8, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);   
	%f(7, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
	%f(6, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);   
	%f(5, T, 15, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache);
end
end



if false
color = 'b';
factors = [1.7];

for factor = factors
    for f = funcs
	basedir = ['/home/other/roald/clusterdata_properVelocityInit_bpRepulsiveCORRECT/hairpinState_lowTemp_gamma_noDih/baseFact',num2str(factor),'/'];
	suff = ['NoDih',num2str(factor)];

	f = f{1};
	T = 10;

	%f(48, T, 100, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
	%
	f(40, T, 85, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
	%f(34, T, 75, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
	%f(28, T, 60, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
	%f(24, T, 45, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
	%f(20, T, 35, sweeps, '*', gamma, dt, dropLogFactor, DLdataStartIndex, loadFromCache, color, basedir, suff);
    end
end
end
