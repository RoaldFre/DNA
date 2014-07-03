zippedThreshold = 0.5; %vanaf N=20: 1.60 +- 0.05
zippedThreshold = 0.8; %vanaf N=20: 1.60 +- 0.06, veel mooiere grafiek, lijkt zelfs recht voor lengtes tot 10~12!

dataStartIndex = 2;

T = 30;
Ns = [10 14 20 28 40 57 80];
Ns = [10 14 20 28 40];
Ns = [5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40, 48, 57, 67, 80];
fitStartN = 20;
fitStartIndex = find(Ns >= fitStartN, 1);
numNs = numel(Ns);
times = cell(numNs, 1);
bounds = cell(numNs, 1);
for i = 1:numNs
	load(["~/DNA/scripts/data/hairpinStateT",num2str(T),"N", num2str(Ns(i))]);
	times{i} = time(dataStartIndex:end); % - time(1);
	bounds{i} = bound(:,dataStartIndex:end);
end

%nuclBoundBPs = round(Ns * 0.5);
%zippedBoundBPs = round(Ns * zippedThreshold);
nuclBoundBPs = (Ns * 0.5);
zippedBoundBPs = (Ns * zippedThreshold);
[exponent, exponentError] = plotZippingLoglogThresholdBootstrap(Ns, times, bounds, nuclBoundBPs, zippedBoundBPs, 'g', fitStartIndex);
