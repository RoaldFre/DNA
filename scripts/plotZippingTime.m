clf;
hold on;

dir="/home/other/roald/clusterdata/hairpin/formation_T20_dt15_time2000_allowUnb1/";
Ns = [4 6 8 10 12 14 16 18 20];
numNs = numel(Ns);
zippingTimes = zeros(numNs, 1);
zippingTimeErrs = zeros(numNs, 1);
for i = 1:numNs
	[time, err] = averageZippingTime([dir, "N", num2str(Ns(i)), "*/formation*"]);
	zippingTimes(i) = time;
	zippingTimeErrs(i) = err;
end

h = loglogerr(Ns, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "linestyle", "none");



dir="/home/other/roald/clusterdata/hairpin/formation_T20_dt15_time2000_allowUnb2/";
Ns = [25 30 35 40];
numNs = numel(Ns);
zippingTimes = zeros(numNs, 1);
zippingTimeErrs = zeros(numNs, 1);
for i = 1:numNs
	[time, err] = averageZippingTime([dir, "N", num2str(Ns(i)), "*/formation*"]);
	zippingTimes(i) = time;
	zippingTimeErrs(i) = err;
end

h = loglogerr(Ns, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "linestyle", "none");



hold off;
