function front = zippingFront(meanState, threshold)

if nargin < 2
	threshold = 0.5;
end

samples = size(meanState)(2);

front = zeros(1, samples);

for i = 1:samples
	aboveInd = find(meanState(:,i) > threshold, 1, 'first');
	above = meanState(aboveInd, i);
	if aboveInd <= 1
		front(i) = 1; % or zero to indicate 'error'?
		continue
	end
	belowInd = aboveInd - 1;
	below = meanState(belowInd, i);

	% linear interpolation to find fractional "threshold index"
	front(i) = belowInd + (aboveInd - belowInd) * (threshold - below) / (above - below);
end
