% function [D, Derr, meanSqDisplacement, sqDisplacement] = plotDiffusion(datafilename, fraction)
%
% datafilename: a file that holds a 3D matrix.
% data(:,:,i) = data of i'th run, which is:
%   [ time(1)  x(1)  y(1)  z(1);
%     time(2)  x(2)  y(2)  z(2);
%      ...     ...   ...   ...  ]
%
% It is assumed that the timedifference is the same everywhere (so the 
% first columns are pretty much redundant)
function [D, Derr, meanSqDisplacement, sqDisplacement] = plotDiffusion(datafilename, fraction)

% command used for plot in report: plotDiffusion("data/strand_dt10_wait50_time1000_N12_g5e12", 0.9999)

if (nargin < 1)
	error("Not enough required arguments!");
end
if (nargin == 1)
	fraction = 0.5; %make every subtrajectory this fraction of the total simulation time
end

data = load(datafilename);
data = data.data;

nRuns = size(data)(3)

time = data(:,1,1);
dt = time(3) - time(2); %XXX HARDCODED FOR SAMPLING AT CONSTANT INTERVAL!
nSamples = numel(time);

nSamplesInResult = floor(nSamples * fraction);
stride = ceil(nSamples * 0.002);
nTrajectories = 1 + floor((nSamples - nSamplesInResult - 1) / stride); %number of sub trajectories

sqDisplacement = zeros(nSamplesInResult, nRuns);
for r = 1:nRuns
	positions = data(:, 2:4, r);
	for i = 1 : stride : (nSamples - nSamplesInResult)
		sqDisplacement(:,r) += squaredDisplacements(positions(i:i+nSamplesInResult, :));
	end
	sqDisplacement(:,r) /= nTrajectories;
end

ts = linspace(dt, nSamplesInResult*dt, nSamplesInResult)';
if (nRuns == 1)
	meanSqDisplacement = sqDisplacement;
	D = mean(meanSqDisplacement ./ (6 * ts))
	Derr = inf;
else
	meanSqDisplacement = mean(sqDisplacement')';
	individualDs = sqDisplacement ./ (6 * (ts * ones(1,nRuns)));
	individualDstds = std(individualDs);
	individualMeanDs = mean(individualDs);
	D = mean(individualMeanDs);
	%TODO check error analysis
	statisticalErrBecauseOfIndividualStds = norm(individualDstds) / sqrt(nRuns - 1)
	statisticalErrBecauseOfStdsInD = std(individualMeanDs)
	Derr = norm([statisticalErrBecauseOfStdsInD, statisticalErrBecauseOfIndividualStds]);
end

[D, Derr]

tsFit = [min(ts), max(ts)];
fit = 6*D*tsFit;


stride = 200;
for i = 1 : numel(ts)/stride
	tsCompressed(i) = ts(i * stride);
	sqDisplacementCompressed(i, :) = sqDisplacement(i * stride, :);
end

clf; hold on;
plot(tsCompressed * 1e6, sqDisplacementCompressed * 1e15)
%plot(ts * 1e6, meanSqDisplacement * 1e15, "g", "linewidth", 4)
plot(tsFit * 1e6, fit * 1e15, "k", "linewidth", 4);
hold off;



filename  = 'diffusion'
caption   = 'Diffusion of a 12 monomer single strand of Adenine. Shown are the squared displacements of 45 individual simulations at 300\,K and with friction $\gamma = 5 \times 10^{12}$\,s$^{-1}$. The black line is the $6Dt$ mean squared displacement curve for the fitted value of $D = (1.3 \pm 0.1) \times 10^{-10}$\,m$^2$/s.';
% XXX caption hardcoded for data/strand_dt10_wait50_time1000_N12_g5e12 dataset! XXX

destdir   = '../report/images';
relImgDir = 'images'; %relative to where your latex project root directory is
ylabrule  = '-1.5cm';
xlab      = 'Time ($\mu$s)';
ylab      = 'Squared displacement ($10^{15}$\,m$^2$)';
width     = '1000';
height    = '800';

%makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
presentationDir = '../presentation/images';
%makeGraphPresentation(filename,presentationDir,xlab,ylab,ylabrule,width,height);
