% Don't garble data when 'not' decimating (ie when decimateFactor == 1)
% Also works for 2D matricess, in which case each column will get decimated
function decimated = decimateWrapper(data, decimateFactor)

addpath('octave-forge');

if (decimateFactor == 1)
	decimated = data;
	return
end

dataSize = size(data)
if (dataSize(1) == 1 || dataSize(2) == 1)
	% It's a 1D array
	decimated = decimate(data, decimateFactor);
else
	% It's a matrix: decimate each column
	for i = 1:numel(data(1,:))
		decimated(:,i) = decimateWrapper(data(:,i), decimateFactor);
		% XXX alloc in advance to avoid quadratic(?) complexity?
	end
end
