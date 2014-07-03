function plotHairpinFullState(N, loadFromCache)

if nargin < 2
	loadFromCache = true;
end

timeShift = 0;
if N == 7
	totTime = 300;
elseif N == 10
	totTime = 300;
elseif N == 14
	totTime = 300;
elseif N == 20
	totTime = 300;
elseif N == 28
	totTime = 500;
elseif N == 40
	totTime = 800;
elseif N == 57
	totTime = 800;
else
	error blah
end



clear bound time state
if not(loadFromCache)
	[bound, time, state] = parseHairpinStateToArray(["/home/other/roald/clusterdata/hairpinState_native/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/relaxT90C_sampleT10C_time",num2str(totTime),"_NXY4/dt15/N",num2str(N),"/*"], true);
	save('-binary', '-z', ["./data/hairpinStateN",num2str(N)], 'bound', 'time', 'state')
	%save('-z', ["./data/hairpinStateN",num2str(N)], 'bound', 'time', 'state')
else
	load(["./data/hairpinStateN",num2str(N)]);
	whos
end



timeCorr = time - time(1) + timeShift;

numRuns = size(state)(1)
numMonomers = size(state)(2);
numSamples = size(state)(3);
meanState = reshape(mean(state), [numMonomers, numSamples]);

%colmap = linspace(1,0,256)' * [1 1 1];
colmap = jet(100);
colormap(colmap);
imagesc(meanState)

hold on;
plot(numMonomers - mean(bound))
hold off;


%save('-binary', '-z', ["./data/hairpinStateN",num2str(N)], 'bound', 'time', 'meanState')
