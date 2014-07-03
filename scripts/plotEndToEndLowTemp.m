more off;

taus = [];
tauErrs = [];
offsets = [];
offsetErrs = [];
amplFracs = [];
amplFracErrs = [];
ampls = [];
amplErrs = [];

rootdir="~/clusterdata/endToEnd/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/T10C/"

%table of simulation Time and Stem lengths
timeSs = [...
          %300,  10;...
          500,  20;...
          500,  24;...
          800,  28;...
          800,  33;...
          2000, 40;...
          2000, 48;...
          4000, 57;...
          ];
decimateFactor = 11;

for i = 1 : numel(timeSs(:,1))
	N = timeSs(i,2)
	[tau, tauStddev, offset, offsetStddev, ampl, amplStddev] = averageEndToEndDist([rootdir, "/dt15_time",num2str(timeSs(i,1)),"/N",num2str(timeSs(i,2)),"/endToEnd*itf11"], decimateFactor);
	taus(i) = tau;
	tauErrs(i) = tauStddev;
	offsets(i) = offset;
	offsetErrs(i) = offsetStddev;
	ampls(i) = ampl;
	amplStddevs(i) = amplStddev;
end

Ns = 2*timeSs(:,2)' + 4; %Total number of monomers!!



[tauCte, tauExponent, tauCteErr, tauExponentErr] = loglogRegression(Ns', taus', 1, 1, tauErrs');
tauC = [tauCte, tauCteErr]
tauE = [tauExponent, tauExponentErr]

[offsetCte, offsetExponent, offsetCteErr, offsetExponentErr] = loglogRegression(Ns', offsets', 1, 1, offsetErrs');
offsetC = [offsetCte, offsetCteErr]
offsetE = [offsetExponent, offsetExponentErr]


fitNs = linspace(min(Ns), max(Ns), 100);
fitTaus = tauCte * fitNs.^tauExponent;
fitOffsets = offsetCte * fitNs.^offsetExponent;


clf;
hold on;

col='b';
h = loglogerr(Ns, taus, tauErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitNs, fitTaus, col, "linewidth", 4);


col='r';
h = loglogerr(Ns, offsets, offsetErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitNs, fitOffsets, col, "linewidth", 4);

%axis([4,100,0,1], 'autoy');

hold off;
