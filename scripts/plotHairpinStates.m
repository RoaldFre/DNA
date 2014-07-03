Ns = [10 20 40];
colors = ['r'; 'g'; 'b'];
timeShifts = [2e-9 1e-8 1e-8];
timeShifts = [0 0 0];
%Ns = [7 10 14 20 28 40 57];
%colors = ['c'; 'r'; 'm'; 'g'; 'k'; 'b'; 'k'];
%timeShifts = [0 0 0 0 0 0 0];

clf;hold on;
for i = 1:numel(Ns)
	N = Ns(i)
	load(["./data/hairpinStateN",num2str(N),"_onlyBound"]);

	time = time - time(1) + timeShifts(i);
	bound = bound - 4; % correct for XY pairs

	% remove runs where XY hasn't zipped yet
	goodRuns = find(bound(:,1) >= 0);
	nRunsAll = numel(bound(:,1))
	nRuns = numel(goodRuns)
	bound = bound(goodRuns, :);
	meanBound = mean(bound);

	loglog(time, meanBound, colors(i));
end


