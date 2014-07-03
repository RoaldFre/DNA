dropLogFactor = 200;
filenamePrefix = './data/hairpinEndToEnd';
variableName = 'dists';
dataIndex = 1;

Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 22, 24, 26, 28, 31, 34, 37, 40, 44, 48];

for T = [10 30]
	equilibriumDists = zeros(1,numel(Ns));
	equilibriumDistsErr = zeros(1,numel(Ns));

	for i = 1:numel(Ns)
		N = Ns(i);
		filename = [filenamePrefix,'T',num2str(T),'N',num2str(N),'droplog',num2str(dropLogFactor)];
		load(filename);

		if not(exist(variableName))
			error(["Can't find the variable with name '",variableName,"' in the data file '",filename,"'!"]);
		end

		eval(['data = ',variableName,';'])
		data = data(:,dataIndex);
		equilibriumDists(i) = mean(data);
		equilibriumDistsErr(i) = std(data) / sqrt(numel(data) - 1);
	end

	filename = ['./data/endToEndEquilibriumScaling','T',num2str(T)];
	asciiData = [Ns(:), equilibriumDists(:), equilibriumDistsErr(:)];
	save('-ascii', filename, 'asciiData');
end





