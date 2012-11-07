function [timesTillZippingFromNucleation, timesTillUnzipping] = parseHairpinFormationAverageBound(filesglob);

files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end

nRuns = numel(files);

zippingSamples = 100000; %TODO must be bigger than largest time 'from nucleation to zipping confirmation'
plotRelaxSamples = 10;

accumUnzippingBound = [0]
accumZippingNuclBound = zeros(zippingSamples, 1);

zippingPlot = figure();
title("zipping")
hold on;
unzippingPlot = figure();
title("unzipping")
hold on;

for run = 1:nRuns
	timesTillZipping(run) = load(files{run}, "timeTillZipping").timeTillZipping;
	timesTillUnzipping(run) = load(files{run}, "timeTillUnzipping").timeTillUnzipping;

	[zippingTime,   zippingState,   zippingBound ...
         relaxTime,     relaxState,     relaxBound, ...
         unzippingTime, unzippingState, unzippingBound] ...
			 = parseHairpinFormationState(files{run});
	N = numel(zippingState(1,:));
	nucleationThreshold = round(N * 0.10);
	nucleationIndex = find(zippingBound < nucleationThreshold, 1, 'last');
	timesTillZippingFromNucleation(run) = zippingTime(end) - zippingTime(nucleationIndex);
	
	% Accumulate the number of zipped base pairs during zipping 
	% (starting from nucleation). Use the relaxation phase data to 
	% continue the data set for a total of 'zippingSamples' samples
	% If the relaxation dataset is too short, pad the data with the 
	% average of the last half of relaxationBound
	zippingBoundFromNucl = zippingBound(nucleationIndex:end);
	paddedWithRelax = [zippingBoundFromNucl; relaxBound];
	if (numel(paddedWithRelax) < numel(accumZippingNuclBound))
		paddedWithRelax = [paddedWithRelax; mean(relaxBound(floor(end/2) : end)) * ones(numel(accumZippingNuclBound) - numel(paddedWithRelax), 1)];
	end
	accumZippingNuclBound += paddedWithRelax(1:numel(accumZippingNuclBound));

	% Accumulate the number of zipped base pairs during unzipping. 
	% Assume that we remain at 0 bound base pairs after the data ends
	if (numel(unzippingBound) > numel(accumUnzippingBound))
		% new data is bigger then accumulated data
		tmp = unzippingBound;
		tmp(1:numel(accumUnzippingBound)) += accumUnzippingBound;
		accumUnzippingBound = tmp;
	else
		accumUnzippingBound(1:numel(unzippingBound)) += unzippingBound;
	end

	figure(unzippingPlot);
	loglog(N -  unzippingBound, 'r')
	figure(zippingPlot);
	%loglog([zippingBound; relaxBound(1:plotRelaxSamples)], 'r')
	plot(zippingBoundFromNucl, 'r')

end

figure(unzippingPlot);
loglog(N - accumUnzippingBound/nRuns)
figure(zippingPlot);
loglog(accumZippingNuclBound/nRuns)

hold off;
