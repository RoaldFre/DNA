addpath "generic"

temperature = 30;
K = 0.02;
%K = 0.005;

filesglob = ["~/clusterdata/eteFreeEnergy/GAGTCAACGTACTGATCACGCTGGATCCTATTTTTAGGATCCAGCGTGATCAGTACGTTGACTC/dt15/T",num2str(temperature),"C/K",num2str(K),"/time1000/*.itf11"]

startIndex = 10000; %TODO check relaxation
nBins = 1000;

T = temperature + 273.15;
beta = 1/(1.3806503e-23 * T);


files = glob(filesglob);
if (isempty(files))
	error "No files match the given glob!"
end
nRuns = numel(files);

hists = zeros(nRuns, nBins);
combinedHist = zeros(nBins, 1);
forces = zeros(nBins, nRuns);
Ks = zeros(nRuns, 1);
Rrefs = zeros(nRuns, 1);
Ubias = zeros(nRuns, nBins);
numSamples = zeros(nRuns, 1);

moreWasOn = page_screen_output;
more off;
maxRref = -Inf;
minRref = Inf;
for run = 1:nRuns
	printf("\rloading Rref from file %d of %d", run, nRuns);
	load(files{run});
	maxRref = max(maxRref, Rref);
	minRref = min(minRref, Rref);
end
deltaR = (maxRref - minRref) / nRuns; % assuming uniform Rref distribution over runs
minR = minRref + 0.25 * deltaR;
maxR = maxRref - 0.25 * deltaR;
printf("\nUsing R from %d to %d Angstrom\n", minR * 1e10, maxR * 1e10);
bins = linspace(minR, maxR, nBins)';

for run = 1:nRuns
	printf("\rloading data from file %d of %d", run, nRuns);
	load(files{run});
	ete = ete(startIndex : end);

	Ubias(run, :) = K/2*(bins - Rref).^2;

	%plot(ete); sleep(1e-9);

	hists(run,:) = hist(ete, bins);
	combinedHist += hists(run,:)';
	m = mean(ete);
	v = var(ete);
	forces(:,run) = (bins - m) / (beta*v) - K * (bins - Rref);
	numSamples(run) = numel(ete);
end
printf("\n");
if moreWasOn
	more on;
end




% UMBRELLA INTEGRATION
% Combine all forces for all windows
force = zeros(nBins, 1);
for run = 1:nRuns
	force += forces(:,run) .* hists(run,:)';
end
force = force ./ combinedHist;

%% XXX THROW AWAY BINS WITH NO SAMPLES!
%goodSamples = find(combinedHist > 0);
%combinedHist = combinedHist(goodSamples);
%force = force(goodSamples);
%bins = bins(goodSamples);
%nBins = numel(bins)

freeEnergy = cumtrapz(bins, force);

plot(bins / 1e-9, freeEnergy * 6.022e20, 'r');
xlabel("end-to-end distance (nm)");
ylabel('Free energy (kJ/mol)');




return




% WEIGHTED HISTOGRAM ANALYSIS METHOD
% http://www.google.be/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDQQFjAA&url=http%3A%2F%2Fmembrane.urmc.rochester.edu%2Fsites%2Fdefault%2Ffiles%2Fwham%2Fwham_talk.pdf&ei=gapRUciANbOy0AXR94DgAg&usg=AFQjCNE1knjZJZybc3HOyAwDiOXqUSPaSw&sig2=y9-GGijC45BG4mYCLH4mAw&bvm=bv.44158598,d.d2k&cad=rja

combinedHistNorm = combinedHist / trapz(bins, combinedHist);

% Initial prob density and F:
P = ones(nBins, 1); P = P/trapz(bins, P);
F = zeros(nRuns, 1);

% Bias weights for each run, matrix (run, bin)
UbiasWeight = exp(-beta * Ubias);

diff = Inf;

warning("off", "Octave:broadcast");
while diff > (1e-2)^2
	prevP = P;
	for j = 1:100
		F = -1/beta * log(UbiasWeight * P); %column vector (run)

		% a(run, bin)
		a = numSamples .* exp(beta*F) .* UbiasWeight; % broadcasting in second .*
		% P(bin)
		P = combinedHistNorm ./ sum(a)';
		P = P/trapz(bins, P);
		%trapz(bins, P)
	end

	plot(bins, P); sleep(1e-9)
	diff = sum(((prevP - P) ./ (prevP + P)).^2);
end
warning("on", "Octave:broadcast");

freeEnergyWHAM = -1/beta * log(P);
