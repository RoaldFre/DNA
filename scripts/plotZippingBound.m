% FOR MIN ZIP, LOW TEMP
function plotZippingBound(N, minTime)

addpath "generic"

dir="/home/other/roald/clusterdata/hairpinFormationLowTemp_native/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/dt15/zipT10C_unzipT90C_allowUnb0.2_nucl0.2/MCsweeps100000_relaxFactor2/"

nucleationBoundBPs = round(N * 0.2);


glob = [dir,'N',num2str(N),'_minZipTime',num2str(minTime),'/*itf11'];
[zippingt, unzippingt, zipFromNuclt, nucleationt, zipping, relax, unzipping] = loadZippingTimeWithNucleationData(glob);

numRuns = numel(zippingt);

numSamplesPerRun = zeros(numRuns, 1);
nucleationIndex = zeros(numRuns, 1);
for r = 1 : numRuns
	nucleationIndex(r) = find(zipping{r,1} - zipping{r,1}(1) >= nucleationt(r))(1);
	numSamplesPerRun(r) = numel(zipping{r,1}(nucleationIndex(r):end)) + numel(relax{r,1});
end
numSamples = min(numSamplesPerRun);

zippingBound = zeros(numRuns, numSamples);
meanZippingState = zeros(numSamples, N);
for r = 1 : numRuns
	zippingBound(r,:) = [zipping{r,2}(nucleationIndex(r):end); relax{r,2}](1 : numSamples)';
	% ignore the base pairs in the loop (so from 1:N instead of 1:(N+2)):
	meanZippingState += [zipping{r,3}(nucleationIndex(r):end,1:N); relax{r,3}(:,1:N)](1:numSamples, :);
end
meanZippingState /= numRuns;

% XXX ASSUMED TO BE THE SAME FOR ALL MEASUREMENTS! XXX
dt = zipping{1}(2) - zipping{1}(1);
time = (0 : (numSamples - 1))*dt;


%loglog(time, mean(zippingBound));


% XXX Correcting for the already "nucleation bounded" base pairs to get SOMEWHAT proper power law without TOO MUCH bias XXX TODO NOT FULLY WITHOUT BIAS!!!!
% N(t) = c t^z
% N(t) - N(t_0) = c (t^z - t_0^z)      !=     c (t - t_0)^z !!!!!
zippingBound = zippingBound - nucleationBoundBPs;



init_t1 = -1; % guard
end_t1 = -1; % guard
if N == 7
	t1=4.0e-9;
	t2=1.8e-8;
elseif N == 10
	%init_t1 = 2e-10;
	%init_t2 = 3e-9;
	%t1=4.0e-9;
	%t2=1.8e-8;
	t1 = 2e-10;
	t2 = 5.5e-9;
elseif N == 14
	init_t1 = 2.5e-10;
	init_t2 = 3e-9;
	t1=4.0e-9;
	t2=2.2e-8;
elseif N == 20
	init_t1 = 1.3e-10;
	init_t2 = 5.1e-9;
	t1 = 8.5e-9;
	t2 = 4.4e-8;
elseif N == 28
	init_t1 = 4.3e-10;
	init_t2 = 1.3e-8;
	t1 = 2.6e-8;
	t2 = 6.8e-8;
elseif N == 40
	%init_t1 = 4.0e-10;
	%init_t2 = 2.4e-8;
	t1 = 7.5e-8;
	t2 = 1.7e-7;
	%end_t1 = 1.6e-7;
	%end_t2 = 2.0e-7;
elseif N == 57
	init_t1 = 4.0e-10;
	init_t2 = 3.3e-8;
	t1 = 5.0e-8;
	t2 = 3.4e-7;
	end_t1 = 3.6e-7;
	end_t2 = 5.0e-7;
else
	error blah
end



meanBound = mean(zippingBound);
errBound = std(zippingBound) / sqrt(numRuns - 1);

errFrac = 1 + errBound ./ meanBound;

%clf;
hold on;
loglog(time, meanBound, 'b');
loglog(time, meanBound .* errFrac, 'r');
loglog(time, meanBound ./ errFrac, 'r');

% last part
is=find(t1<time & time<t2);
xs=time(is);
ys=meanBound(is);
yerrs=errBound(is)';
[cte, exponent, cteStddev, exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
loglog(xs, cte*xs.^exponent, 'k', 'linewidth', 3)


% initial part
if init_t1 > 0
	is=find(init_t1<time & time<init_t2);
	xs=time(is);
	ys=meanBound(is);
	yerrs=errBound(is)';
	[init_cte, init_exponent, init_cteStddev, init_exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
	loglog(xs, init_cte*xs.^init_exponent, 'g', 'linewidth', 3)
end

% final part
if end_t1 > 0
	is=find(end_t1<time & time<end_t2);
	xs=time(is);
	ys=meanBound(is);
	yerrs=errBound(is)';
	[end_cte, end_exponent, end_cteStddev, end_exponentStddev] = loglogRegression(xs, ys', 1e5, 0.6, yerrs)
	loglog(xs, end_cte*xs.^end_exponent, 'm', 'linewidth', 3)
end


hold off;


%figure
%clf; hold on;
%plot(time, meanBound, 'b');
%plot(time, meanBound + errBound, 'r');
%plot(time, meanBound - errBound, 'r');
%hold off;

if init_t1 > 0
	[1/init_exponent, init_exponentStddev/init_exponent^2]
end
[1/exponent, exponentStddev/exponent^2]
if end_t1 > 0
	[1/end_exponent, end_exponentStddev/end_exponent^2]
end



return

figure
imagesc(meanZippingState');
figure
hold on
plot(meanZippingState(1,:), 'r')
plot(meanZippingState(round(end*0.25),:), 'g')
plot(meanZippingState(round(end*0.5),:), 'b')
hold off
