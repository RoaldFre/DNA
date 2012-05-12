function plotDiffusion(filename)

data=load(filename);

time = data(:,1);
rsq = data(:,2:end);

if (size(rsq)(2) == 1)
	rsqmean = rsq;
else
	rsqmean = mean(rsq')';
end


D = mean(rsqmean(2:end) ./ (6 * time(2:end)))

fit = 6*D*time;

plot(time, rsqmean, time, fit, 'g');
