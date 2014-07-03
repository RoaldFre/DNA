function plotHairpinStatesWithError(Ns, loadFromCache, timeShiftFactor)

if nargin < 2
	loadFromCache = false;
end
if nargin < 3
	timeShiftFactor = 0;
end

for N = Ns
	plotHairpinState(N, loadFromCache, timeShiftFactor)
	sleep(1e-9);
end
