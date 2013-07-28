filenamePrefix = './data/hairpinState';
variableName = 'bound';
plotopt.graphFilePrefix = 'state';

plotopt.destDir = '../../thesis/latex/images/plots/';
plotopt.relImgDir = 'images/plots';
plotopt.measurement = 'the number of bound base pairs';
plotopt.quantity = '\nb';
plotopt.units = '';
plotopt.unitsSq = '';

dataStartIndex = 2;
withErrorBars = false;

for T = [10 30];
	for squaredDeviation = [true false]

		if T == 10
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48];
		elseif T == 30
			Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57, 67, 80];
		else
			error "unknown T"
		end

		makeHairpinPlots(Ns, T, filenamePrefix, variableName, squaredDeviation, plotopt, withErrorBars, dataStartIndex)
	end
end
