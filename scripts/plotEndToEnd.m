%dir="/home/other/roald/clusterdata/endToEnd/CTTTTAAACGCAACTGGAGCCGAACCTACAGGAGATAGACGTCTGTTCCGGGAAGGTGCAGCGCTAGAGTATATATGTGACAGTCTGAAGTACTGCAGGG/T80/dt15_time200/";
dir="~/clusterdata/endToEnd/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/T90C"

Ss = 30:10:50; % Stem lengths (ignore S=20)
decimateFactor = 10;

taus = zeros(size(Ss));
tauErrs = zeros(size(Ss));
offsets = zeros(size(Ss));
offsetErrs = zeros(size(Ss));
for i = 1 : numel(Ss)
	Ss(i)
	[tau, tauStddev, offset, offsetStddev] = averageEndToEndDist([dir, "/dt15_time500/N", num2str(Ss(i)), "/endToEnd*itf11"], decimateFactor);
	taus(i) = tau;
	tauErrs(i) = tauStddev;
	offsets(i) = offset;
	offsetErrs(i) = offsetStddev;
end

% Quick addition: from N=60 onwards, we have simulated a longer period, so 
% append these results separately here.
extraSs = 60:10:100;
for i = 1 : numel(extraSs)
	extraSs(i)
	[tau, tauStddev, offset, offsetStddev] = averageEndToEndDist([dir, "/dt15_time2000/N", num2str(extraSs(i)), "/endToEnd*itf11"], decimateFactor);
	taus = [taus tau];
	tauErrs = [tauErrs tauStddev];
	offsets = [offsets offset];
	offsetErrs = [offsetErrs offsetStddev];
end

Ss = [Ss extraSs];

Ns = 2*Ss + 4; %Total number of monomers!!





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

col='b'
%col='c'
h = loglogerr(Ns, taus, tauErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitNs, fitTaus, col, "linewidth", 4);

col='r'
%col='m'
h = loglogerr(Ns, offsets, offsetErrs);
set(h, "marker", ".");
set(h, "color", col);
set(h, "linestyle", "none");
set(h, "linewidth", 4);

loglog(fitNs, fitOffsets, col, "linewidth", 4);

%axis([4,100,0,1], 'autoy');

hold off;
