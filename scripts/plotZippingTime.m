Ns = [];

dir="/home/other/roald/clusterdata/hairpin/formation_T20_dt15_time2000_allowUnb1/";
theseNs = [4 6 8 10 12 14 16 18 20];
prevNum = numel(Ns);
Ns = [Ns theseNs];
currNum = numel(Ns);
for i = (prevNum + 1) : currNum
	[time, err] = averageZippingTime([dir, "N", num2str(Ns(i)), "*/formation*"]);
	zippingTimes(i) = time;
	zippingTimeErrs(i) = err;
end



dir="/home/other/roald/clusterdata/hairpin/formation_T20_dt15_time2000_allowUnb2/";
theseNs = [25 30 35 40 45 50];
prevNum = numel(Ns);
Ns = [Ns theseNs];
currNum = numel(Ns);
for i = (prevNum + 1) : currNum
	[time, err] = averageZippingTime([dir, "N", num2str(Ns(i)), "*/formation*"]);
	zippingTimes(i) = time;
	zippingTimeErrs(i) = err;
end




[cte, exponent, cteErr, exponentErr] = loglogregression(Ns, zippingTimes, 1, 1, zippingTimeErrs)
[cte, cteErr]
[exponent, exponentErr]


%[cte, exponent, shift, cteErr, exponentErr, shiftErr] = loglogregressionWithShift(Ns, zippingTimes, 1, 1, zippingTimeErrs)
%[cte, cteErr]
%[exponent, exponentErr]
%[shift, shiftErr]


fitNs = linspace(min(Ns), max(Ns), 100);
fitTimes = cte * fitNs.^exponent;


close all; clf;
hold on;

h = loglogerr(Ns, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "linestyle", "none");

loglog(fitNs, fitTimes);

hold off;


figure;
hold on;
h = errorbar(Ns, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "linestyle", "none");

plot(fitNs, fitTimes);
hold off;


