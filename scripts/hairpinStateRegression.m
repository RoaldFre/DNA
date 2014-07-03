function hairpinStateRegression(T, NsAndTimes, threshold, sweeps)

for i = 1:numel(NsAndTimes(:,1))
	N = NsAndTimes(i,1);
	totTime = NsAndTimes(i,2);

	dir = ["/home/other/roald/clusterdata/hairpinState_lowTemp/T",num2str(T),"C_time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_NXY4/dt15/N",num2str(N)]
	%dir = ["/home/other/roald/clusterdata/hairpinState_lowTemp/T",num2str(T),"C_time",num2str(totTime),"_sweeps",num2str(sweeps,20),"_relaxFactor100_NXY4/dt15/N",num2str(N)]
	if not(loadFromCache)
		[bound, time] = parseHairpinStateToArray([dir, '/*.itf11']);
		save('-binary', '-z', [dir,"/hairpinStateN",num2str(N),"_onlyBound"], 'bound', 'time')
	else
		load([dir,"/hairpinStateN",num2str(N),"_onlyBound"]);
	end


	time = time - time(1);

	bound = bound - 4; % correct for XY pairs

	% remove runs where XY hasn't zipped yet
	goodRuns = find(bound(:,1) >= 0);
	nRunsAll = numel(bound(:,1))
	nRuns = numel(goodRuns)
	bound = bound(goodRuns, :);

	meanBound = mean(bound);


	nucleationThreshold = N * threshold;
	nucleationIndex = find(meanBound <= nucleationThreshold, 1, 'last');
	nucleationTime = time(nucleationIndex);



	eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
end

	
