% N = number of matching base pairs in the stem
function [timesTillZippingFromNucleation] = parseHairpinFormationAverageBoundFromCell(file, N);

more off;

load(file); % gives cell array 'data'

nRuns = size(data)(1);

% determine longest run (in number of samples) for zipping and unzipping
zippingSamples = 0;
unzippingSamples = 0;
for run = 1:nRuns
	zippingBound = data{run, 3};
	relaxBound = data{run, 6};
	zippingSamples = max(zippingSamples, numel(zippingBound) + numel(relaxBound));
	unzippingBound = data{run, 9};
	unzippingSamples = max(unzippingSamples, numel(unzippingBound));
end

accumUnzippingBound = zeros(unzippingSamples, 1);
accumZippingBound = zeros(zippingSamples, 1);

zippingPlot = figure();
title("zipping")
hold on;
unzippingPlot = figure();
title("unzipping")
hold on;

printf("\n");
for run = 1:nRuns
	%zippingTime    = data{run, 1};
	%zippingState   = data{run, 2};
	zippingBound   = data{run, 3};
	%relaxTime      = data{run, 4};
	%relaxState     = data{run, 5};
	relaxBound     = data{run, 6};
	%unzippingTime  = data{run, 7};
	%unzippingState = data{run, 8};
	unzippingBound = data{run, 9};

	% Accumulate the number of zipped base pairs during zipping 
	% (starting from nucleation). Use the relaxation phase data to 
	% continue the data set for a total of 'zippingSamples' samples
	% If the relaxation dataset is too short, pad the data with the 
	% average of the last half of relaxationBound
	paddedWithRelax = [zippingBound; relaxBound];
	paddedWithRelax = [paddedWithRelax; mean(relaxBound(floor(end/2) : end)) * ones(numel(accumZippingBound) - numel(paddedWithRelax), 1)];
	accumZippingBound += paddedWithRelax;

	% Accumulate the number of zipped base pairs during unzipping. 
	% Assume that we remain at 0 bound base pairs after the data ends
	accumUnzippingBound(1:numel(unzippingBound)) += unzippingBound;

	figure(unzippingPlot);
	plot(N - unzippingBound - 4, 'r')
	figure(zippingPlot);
	plot(zippingBound - 4, 'r')

	printf("\rLoaded run %d of %d", run, nRuns);
end
printf("\n");

figure(unzippingPlot);
plot(N - accumUnzippingBound/nRuns - 4);
figure(zippingPlot);
plot(accumZippingBound/nRuns - 4);

hold off;
more on;
