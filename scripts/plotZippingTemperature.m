function plotZippingTemperature(middle, N, Ts, color)
dir="/home/other/roald/clusterdata/hairpinMeltingSeq2FlorJoy/";
%dir="/home/other/roald/clusterdata/hairpinMeltingSeq2/";

start = 2; %ignore first reading (not yet relaxed)

clear avgFractions;
clear errFractions;
for i = 1:numel(Ts)
	T = Ts(i);
	filesglob = ([dir,"/T",num2str(T),"/*/",middle,"/meltingTemp*"]);
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
