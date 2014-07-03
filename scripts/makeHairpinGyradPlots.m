variableName = 'gyrad';
filenamePrefix = './data/hairpinGyrad';
plotopt.graphFilePrefix = 'gyrad';

plotopt.destDir = '../../thesis/latex/images/plots';
plotopt.relImgDir = 'images/plots';
plotopt.measurement = 'the gyration radius';
plotopt.quantity = '\Rg';
plotopt.units = 'm';
plotopt.unitsSq = 'm$^2$';

dataStartIndex = 2;
withErrorBars = false;

for T = [10 30];
	for squaredDeviation = [true false]

		if T == 10
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48];
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48];
		elseif T == 30
			%Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57, 67, 80];
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48];
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48];
		else
			error "unknown T"
		end

		makeHairpinPlots(Ns, T, filenamePrefix, variableName, squaredDeviation, plotopt, withErrorBars, dataStartIndex)
	end
end
