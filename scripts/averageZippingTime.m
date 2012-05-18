function [meanZip, errZip, meanUnzip, errUnzip] = averageZippingTime(filesglob)

[timesTillZipping, timesTillUnzipping, temperature, temperatures, data, measureTime, timestep, allowedUnboundBPs, allowedBoundBps] = parseHairpins(filesglob, false);

meanZip = mean(timesTillZipping);
errZip = std(timesTillZipping)/sqrt(numel(timesTillZipping));

meanUnzip = mean(timesTillUnzipping);
errUnzip = std(timesTillUnzipping)/sqrt(numel(timesTillUnzipping));
