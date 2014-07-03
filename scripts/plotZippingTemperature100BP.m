function plotZippingTemperature100BP(N, Ts, color)
if nargin < 3
	color = 'b';
end
dir="/home/other/roald/clusterdata/melting/CTTTTAAACGCAACTGGAGCCGAACCTACAGGAGATAGACGTCTGTTCCGGGAAGGTGCAGCGCTAGAGTATATATGTGACAGTCTGAAGTACTGCAGGG/meltingTemp_dt15_wait1000_measTime500_relaxTime50_Tstart100_Tstep-10_nSteps8_salt115/";

start = 1;

clear avgFractions;
clear errFractions;
for i = 1:numel(Ts)
	T = Ts(i);
	filesglob = ([dir,"/N",num2str(N),"/meltingTemp*"]);
	[avgFraction, errFraction, temperatures] = zippingTemperature(filesglob, N);
	avgFraction = avgFraction(start:end);
	errFraction = errFraction(start:end);
	avgFractions(:,i) = avgFraction;
	errFractions(:,i) = errFraction;
	combinedFraction(i) = mean(avgFraction);
	combinedErr(i) = norm([norm(errFraction); std(avgFraction)])  / sqrt(numel(errFraction));
end

%hold on;
%for i = 1:numel(avgFractions(:,1))
%	errorbar(Ts, avgFractions(i,:), errFractions(i,:));
%end
%hold off

%[A, B, Astddev, Bstddev] = meltingTempRegression(Ts, combinedFraction, 0.05, 50, combinedErr)
%[A, B, Astddev, Bstddev] = meltingTempRegression(Ts, combinedFraction, 0.05, 50)

%for knotts (wrong):
[A, B, Astddev, Bstddev] = meltingTempRegression(Ts, combinedFraction, 0.05, 10)

TsFit = linspace(0, 100, 100);
fit = 0.5 - 0.5*tanh(A*(TsFit - B));

hold on;
plot(TsFit, fit, color, "linewidth", 4);
%hold off;

h1 = errorbar(Ts, combinedFraction, combinedErr);
set(h1, "color", color);
set(h1, "marker", ".");
set(h1, "linestyle", "none");
set(h1, "linewidth", 5);
