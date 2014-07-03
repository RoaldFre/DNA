more off;

tau1s = [];
tau2s = [];
tau1Errs = [];
tau2Errs = [];
offsets = [];
offsetErrs = [];
amplFracs = [];
amplFracErrs = [];

rootdir="~/clusterdata/endToEndMC/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/T10C/"

sweeps = 1000;
Ss = [10 14 20 24 40];
decimateFactor = 7;

for i = 1 : numel(Ss)
	N = Ss(i)
	glob = [rootdir, "/MCsweeps",num2str(sweeps),"/N",num2str(Ss(i)),"/endToEnd*itf11"];
	%FIXED INITIAL LENGTH
	%[tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, ampl1, ampl1Stddev, ampl2, ampl2Stddev] = averageEndToEndDist2timesMC(glob, decimateFactor)
	[tau1, tau1Stddev, tau2, tau2Stddev, offset, offsetStddev, amplFrac, amplFracStddev] = averageEndToEndDist2timesMC(glob, decimateFactor);
	%[tau, tauStddev, offset, offsetStddev, ampl, amplFracStddev] = averageEndToEndDist2timesMC(glob, decimateFactor);
	%taus(i) = tau;
	%tauErrs(i) = tauStddev;
	tau1s(i) = tau1;
	tau1Errs(i) = tau1Stddev;
	tau2s(i) = tau2;
	tau2Errs(i) = tau2Stddev;
	offsets(i) = offset;
	offsetErrs(i) = offsetStddev;
	%amplFracs(i) = amplFrac;
	%amplFracErrs(i) = amplFracStddev;
end

Ns = 2*Ss + 4; %Total number of monomers!!


tau1s
tau2s
mean(tau1s)
mean(tau2s)
