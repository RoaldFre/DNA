Ns = [];

dir="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time10000/";
Ns = [5 10 15 20 30 40];
for i = 1 : numel(Ns)
	[zipTime, zipErr, unzipTime, unzipErr] = averageZippingTime([dir, "N", num2str(Ns(i)), "*/formation*"]);
	zippingTimes(i) = zipTime;
	zippingTimeErrs(i) = zipErr;
	unzippingTimes(i) = unzipTime;
	unzippingTimeErrs(i) = unzipErr;
end

[Ns; zippingTimes; unzippingTimes]


[zipCte, zipExponent, zipCteErr, zipExponentErr] = loglogRegression(Ns', zippingTimes', 1, 1, zippingTimeErrs')
[zipCte, zipCteErr]
[zipExponent, zipExponentErr]

[unzipCte, unzipExponent, unzipCteErr, unzipExponentErr] = loglogRegression(Ns', unzippingTimes', 1, 1, unzippingTimeErrs')
[unzipCte, unzipCteErr]
[unzipExponent, unzipExponentErr]


fitNs = linspace(min(Ns), max(Ns), 100);
fitZipTimes = zipCte * fitNs.^zipExponent;
fitUnzipTimes = unzipCte * fitNs.^unzipExponent;


hold on;

h = loglogerr(Ns, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "linestyle", "none");

loglog(fitNs, fitZipTimes);

h = loglogerr(Ns, unzippingTimes, unzippingTimeErrs);
set(h, "marker", ".");
set(h, "color", "r");
set(h, "linestyle", "none");

loglog(fitNs, fitUnzipTimes, "r");

hold off;





return





figure;
hold on;
h = errorbar(Ns, zippingTimes, zippingTimeErrs);
set(h, "marker", ".");
set(h, "linestyle", "none");

plot(fitNs, fitZipTimes);

h = errorbar(Ns, unzippingTimes, unzippingTimeErrs);
set(h, "marker", ".");
set(h, "linestyle", "none");

plot(fitNs, fitunzipTimes, "r");

hold off;


