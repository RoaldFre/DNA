function [avgFraction, errFraction, temperatures] = zippingTemperature(filesglob, N)

[timesTillZipping, temperature, temperatures, data, measureTime, timestep, allowedUnboundBPs] = parseHairpins(filesglob, true);

numTemps = numel(temperatures)
numRuns = numel(data(1,1,1,:))

for i = 1:numTemps
	for r = 1:numRuns
		fractionZipped(i,r) = mean(data(2,:,i,r));
	end
end

avgFraction = zeros(numTemps, 1);
errFraction = zeros(numTemps, 1);
for i = 1:numTemps
	avgFraction(i) = mean(fractionZipped(i,:));
	errFraction(i) = std(fractionZipped(i,:)) / sqrt(numRuns);
end

avgFraction /= N;
errFraction /= N;
