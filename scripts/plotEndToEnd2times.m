more off;

tau1s = [];
tau2s = [];
tau1Errs = [];
tau2Errs = [];
offsets = [];
offsetErrs = [];
amplFracs = [];
amplFracErrs = [];

dir="~/clusterdata/endToEnd/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/T90C"

Ss = 30:10:50; % Stem lengths (ignore S=20)
Ss = []
decimateFactor = 11;

for i = 1 : numel(Ss)
	Ss(i)
	[tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, amplFrac, amplFracStddev] = averageEndToEndDist2times([dir, "/dt15_time500/N", num2str(Ss(i)), "/endToEnd*itf11"], decimateFactor);
	tau1s(i) = tau1;
	tau1Errs(i) = tau1Stddev;
	tau2s(i) = tau2;
	tau2Errs(i) = tau2Stddev;
	offsets(i) = offset;
	offsetErrs(i) = offsetStddev;
	amplFracs(i) = amplFrac;
	amplFracErrs(i) = amplFracStddev;
end

% Quick addition: from N=60 onwards, we have simulated a longer period, so 
% append these results separately here.
extraSs = 60:10:100;
for i = 1 : numel(extraSs)
	extraSs(i)
	[tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, amplFrac, amplFracStddev] = averageEndToEndDist2times([dir, "/dt15_time2000/N", num2str(extraSs(i)), "/endToEnd*itf11"], decimateFactor);
	tau1s = [tau1s tau1];
	tau1Errs = [tau1Errs tau1Stddev];
	tau2s = [tau2s tau2];
	tau2Errs = [tau2Errs tau2Stddev];
	offsets = [offsets offset];
	offsetErrs = [offsetErrs offsetStddev];
	amplFracs = [amplFracs amplFrac];
	amplFracErrs = [amplFracErrs amplFracStddev];
end

Ss = [Ss extraSs];

Ns = 2*Ss + 4; %Total number of monomers!!





[tau1Cte, tau1Exponent, tau1CteErr, tau1ExponentErr] = loglogRegression(Ns', tau1s', 1, 1, tau1Errs');
tau1C = [tau1Cte, tau1CteErr]
tau1E = [tau1Exponent, tau1ExponentErr]

[tau2Cte, tau2Exponent, tau2CteErr, tau2ExponentErr] = loglogRegression(Ns', tau2s', 2, 2, tau2Errs');
tau2C = [tau2Cte, tau2CteErr]
tau2E = [tau2Exponent, tau2ExponentErr]

[offsetCte, offsetExponent, offsetCteErr, offsetExponentErr] = loglogRegression(Ns', offsets', 1, 1, offsetErrs');
offsetC = [offsetCte, offsetCteErr]
offsetE = [offsetExponent, offsetExponentErr]


fitNs = linspace(min(Ns), max(Ns), 100);
fitTau1s = tau1Cte * fitNs.^tau1Exponent;
fitTau2s = tau2Cte * fitNs.^tau2Exponent;
fitOffsets = offsetCte * fitNs.^offsetExponent;


clf;
hold on;

col='b';
h = loglogerr(Ns, tau1s, tau1Errs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitNs, fitTau1s, col, "linewidth", 4);


col='g';
h = loglogerr(Ns, tau2s, tau2Errs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitNs, fitTau2s, col, "linewidth", 4);


col='r';
h = loglogerr(Ns, offsets, offsetErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitNs, fitOffsets, col, "linewidth", 4);

%axis([4,100,0,1], 'autoy');

hold off;
