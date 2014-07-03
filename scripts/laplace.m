% simple O(n^2) implementation of laplace transforms
function Fs = laplace(ts, fs, ss)

%ss=linspace(5/(time(end)-time(1)), 5/(time(2)-time(1)), 10000)

if numel(fs) != numel(ts)
	error "time and function value arrays must have same length"
end

n = numel(ss);
Fs = zeros(size(ss));
for i=1:n
	s = ss(i);
	Fs(i) = sum(fs .* exp(-s * ts));
end
