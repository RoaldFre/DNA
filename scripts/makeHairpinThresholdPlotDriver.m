plotopt.destDir = '../../thesis/latex/images/plots/';
plotopt.relImgDir = 'images/plots';
withErrorBars = false;


withErrorBars = false;
dataStartIndex = 2;
dropLogFactor = 200;

%fitStartN = 20;
%fitStartN = 24;

for T = [-10 10 30];
#for T = [10 30];
%for T = 10
    if T == -10
	    fitStartN = 24;
    else
	    fitStartN = 24;
    end
    %suffs = {'','NoDih1.4','NoDih1.7'}
    suffs = {'NoDih1.4','NoDih1.7'}
    for a = 1:numel(suffs);
	suff = suffs{a};
	noDih = not(isempty(suff))
	if noDih
		if T == 10
			Ns = [20 24 28 34 40 48];
		else
			continue
		end
		plotopt.additionalInfoStr = ['Dihedral interactions were disabled and base pairing was increased by a factor of ',num2str(factor),'. '];
	else
		if T == -10
			Ns = [20, 24, 28, 34, 40];
		elseif T == 10
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48];
			Ns = [12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48];
			% TODO QUICK TEST TO SPEED THINGS UP!!
			%Ns = [20, 24, 28];
		elseif T == 30
			%Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48, 57, 67, 80];
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48];
			Ns = [12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48];
		else
			error "unknown T"
		end
		suff = '';
		plotopt.additionalInfoStr = ''
	end
	numNs = numel(Ns);




	filenamePrefix = ['./data/hairpinTotalEnergy',suff];
	variableName = 'totalEnergy';
	plotopt.graphFilePrefix = ['energy',suff];
	plotopt.measurement = 'the base pairing energy';
	plotopt.quantity = '-\Vbp';
	plotopt.units = 'J';
	plotopt.unitsSq = 'J$^2$';
	squaredDeviation = false;
	thresholdFractions = [0.8 0.7 0.6 0.5 0.4];
	thresholdParams = thresholdFractions;
	findTauIndex = @(time, energy, fraction, N) find(energy >= mean(energy(end-3:end))*fraction, 1, 'first');

	thresholdInformationString = 'thresholdInformationString';
	thresholdLegendStrings = arrayfun(@(r) ['$\lambda = ',num2str(r),'$'], thresholdFractions, 'UniformOutput', false);
	plotTauYs = true;

	if noDih
		xaxis = [1e-9, 1e-7];
		yaxis = [1e-19, 3e-18];
		insetAxes = [10 60 1e-9 1e-7]
	else
		if T == 10 || T == -10
			xaxis = [1e-9, 1e-7];
			yaxis = [6e-20, 1.2e-18];
			insetAxes = [10 60 1e-9 1e-7]
		elseif T == 30
			xaxis = [6e-10, 1.5e-7];
			yaxis = [4e-20, 1.0e-18];
			insetAxes = [10 60 1e-9 1e-7]
		end
	end

	makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs, xaxis, yaxis, insetAxes);








	
	filenamePrefix = ['./data/hairpinState',suff];
	variableName = 'bound';
	plotopt.graphFilePrefix = ['state',suff];
	plotopt.measurement = 'the number of bound base pairs';
	plotopt.quantity = '\nb';
	plotopt.units = '';
	plotopt.unitsSq = '';
	squaredDeviation = false;
	thresholdFractions = [0.8 0.7 0.6 0.5 0.4];
	thresholdParams = thresholdFractions;
	findTauIndex = @(time, bound, fraction, N) find(bound >= N*fraction, 1, 'first');

	thresholdInformationString = 'thresholdInformationString';
	thresholdLegendStrings = arrayfun(@(r) ['$\lambda = ',num2str(r),'$'], thresholdFractions, 'UniformOutput', false);
	plotTauYs = true;

	if noDih
		xaxis = [3e-9, 1e-7];
		yaxis = [7, 50]
		insetAxes = [10 60 1e-9 1e-7]
	else
		if T == 10 || T == -10
			xaxis = [1e-9, 1e-7];
			yaxis = [3, 60]
			insetAxes = [10 60 1e-9 1e-7]
		elseif T == 30
			xaxis = [6e-10, 1.5e-7];
			yaxis = [2, 60]
			insetAxes = [10 60 1e-9 1e-7]
		end
	end

	makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs, xaxis, yaxis, insetAxes);










	variableName = 'dists';
	filenamePrefix = ['./data/hairpinEndToEnd',suff];
	plotopt.graphFilePrefix = ['endToEnd',suff];
	plotopt.measurement = 'the end-to-end distance';
	plotopt.quantity = '\Rete';
	plotopt.units = 'm';
	plotopt.unitsSq = 'm$^2$';

if 0 % NON SQUARE DEVIATION
	squaredDeviation = false;

	%thresholdRete = [1.5e-9, 1e-9, 2e-9, 3e-9];
	%thresholdParams = thresholdRete;
	%thresholdLegendStrings = arrayfun(@(x)['$',num2tex(x),'$'], thresholdRete, 'UniformOutput', false);
	%findTauIndex = @(time, Rete, threshold, N) find(Rete <= threshold, 1, 'first');

	thresholdFractions = [0.5 0.4 0.3 0.2 0.1 0.05 0.02];
	thresholdParams = thresholdFractions;
	thresholdLegendStrings = arrayfun(@(x)['$\lambda = ',num2tex(x),'$'], thresholdFractions, 'UniformOutput', false);
	findTauIndex = @(time, Rete, fraction, N) find((Rete - mean(Rete(end-3:end))) <= (Rete(1) - mean(Rete(end-3:end)))*fraction, 1, 'first');

	thresholdInformationString = 'thresholdInformationString';
	plotTauYs = true;

	if T == 10 || T == -10
		xaxis = [1e-11, 1e-7];
		yaxis = [5e-10, 2e-8];
	elseif T == 30
		xaxis = [1e-11, 1.5e-7];
		yaxis = [5e-10, 2e-8];
	end
	
	%insetPosition = [0.17 0.62 0.35 0.3];
	insetPosition = [0.17 0.1 0.35 0.3];
	legendLocation = 'east';


	makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs, xaxis, yaxis, insetAxes, insetPosition, legendLocation);

	% SINGLE CONSTANT THRESHOLD
	plotopt.graphFilePrefix = ['endToEndConstantThreshold',suff];
	thresholdParams = [1.3 1.5 1.7 1.9]*1e-9;
	thresholdLegendStrings = arrayfun(@(x)['$\lambda = ',num2tex(x),'$'], thresholdFractions, 'UniformOutput', false);
	findTauIndex = @(time, Rete, threshold, N) find(Rete <= threshold, 1, 'first');

	thresholdInformationString = 'thresholdInformationString';
	plotTauYs = true;

	if T == 10 || T == -10
		xaxis = [1e-11, 1e-7];
		yaxis = [5e-10, 2e-8];
	elseif T == 30
		xaxis = [1e-11, 1.5e-7];
		yaxis = [5e-10, 2e-8];
	end
	
	%insetPosition = [0.17 0.62 0.35 0.3];
	insetPosition = [0.17 0.18 0.35 0.3];
	legendLocation = 'southeast';


	makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs, xaxis, yaxis, insetAxes, insetPosition, legendLocation);
end % END OF NON SQUARED DEVIATION



	variableName = 'dists';
	filenamePrefix = ['./data/hairpinEndToEnd',suff];
	plotopt.graphFilePrefix = ['endToEnd',suff];
	plotopt.measurement = 'the end-to-end distance';
	plotopt.quantity = '\Rete';
	plotopt.units = 'm';
	plotopt.unitsSq = 'm$^2$';

	% SQUARED DEVIATION

	squaredDeviation = true;

	thresholdFractions = [0.9 0.8 0.7 0.6 0.5];
	thresholdParams = thresholdFractions;
	thresholdLegendStrings = arrayfun(@(x)['$\lambda = ',num2tex(x),'$'], thresholdFractions, 'UniformOutput', false);
	findTauIndex = @(time, ReteSqDev, fraction, N) find(ReteSqDev >= mean(ReteSqDev(end-3:end))*fraction, 1, 'first');

	if noDih
		xaxis = [3e-9, 1e-7];
		yaxis = [1e-17, 3e-16];
	else
		if T == 10 || T == -10
			xaxis = [2e-9, 1e-7];
			yaxis = [1e-17, 3e-16];
			insetAxes = [10 60 2e-9 1e-7]
		elseif T == 30
			xaxis = [2e-9, 1.5e-7];
			yaxis = [1e-17, 3e-16];
			insetAxes = [10 60 3e-9 2e-7]
		end
	end

	makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs, xaxis, yaxis, insetAxes);










	variableName = 'gyrad';
	filenamePrefix = ['./data/hairpinGyrad',suff];
	plotopt.graphFilePrefix = ['gyrad',suff];
	plotopt.measurement = 'the gyration radius';
	plotopt.quantity = '\Rg';
	plotopt.units = 'm';
	plotopt.unitsSq = 'm$^2$';

	squaredDeviation = false;

	thresholdParams = 1;
	thresholdLegendStrings = {};
	findTauIndex = @(time, gyrad, _, N) find(gyrad == min(gyrad), 1, 'first');

	thresholdInformationString = 'thresholdInformationString';
	%plotTauYs = false;
	plotTauYs = true;

	if noDih
		xaxis = [1e-11, 1e-7];
		yaxis = [];
	else
		if T == 10 || T == -10
			xaxis = [1e-11, 1e-7];
			yaxis = [];
			insetAxes = [10 60 3e-10 7e-8]
		elseif T == 30
			xaxis = [1e-11, 1.5e-7];
			yaxis = [];
			insetAxes = [10 60 3e-10 7e-8]
		end
	end

	makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs, xaxis, yaxis, insetAxes);


	% GYRAD WITH SQDEV

%	squaredDeviation = true;
%
%	thresholdFractions = [0.9 0.8 0.7 0.6 0.5];
%	thresholdParams = thresholdFractions;
%	thresholdLegendStrings = arrayfun(@(x)['$\lambda = ',num2tex(x),'$'], thresholdFractions, 'UniformOutput', false);
%	findTauIndex = @(time, gyradSqDev, fraction, N) find(gyradSqDev >= mean(gyradSqDev(end-3:end))*fraction, 1, 'first');
%	% TODO: probleem: S=48 is onvoldoende ver gesimuleerd -- Rg nog niet ge-equilibreerd op einde!!
%
%	thresholdInformationString = 'thresholdInformationString';
%	%plotTauYs = false;
%	plotTauYs = true;
%
%	%yaxis = [0.5, 100]
%	%if T == 10
%	%	xaxis = [1e-10, 1e-7];
%	%elseif T == 30
%	%	xaxis = [1e-10, 1e-7];
%	%end
%
%	makeHairpinThresholdPlot(Ns, thresholdParams, findTauIndex, fitStartN, T, filenamePrefix, variableName, dropLogFactor, squaredDeviation, plotopt, thresholdInformationString, thresholdLegendStrings, withErrorBars, dataStartIndex, plotTauYs); %, xaxis, yaxis, insetAxes);
	
	
    end
end
