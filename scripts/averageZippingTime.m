function [mean, err] = averageZippingTime(filesglob)

[timesTillZipping, temperature, temperatures, data, measureTime, timestep, allowedUnboundBPs] = parseHairpins(filesglob, false);

mean = mean(timesTillZipping);
err = std(timesTillZipping)/sqrt(numel(timesTillZipping));
