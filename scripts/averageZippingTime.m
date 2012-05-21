function [meanZip, errZip, meanUnzip, errUnzip] = averageZippingTime(filesglob)

[timesTillZipping, timesTillUnzipping, temperature, temperatures, data, timestep, allowedUnboundBPs, allowedBoundBps] = parseHairpinFormation(filesglob, false);

meanZip = mean(timesTillZipping);
errZip = std(timesTillZipping)/sqrt(numel(timesTillZipping));

meanUnzip = mean(timesTillUnzipping);
errUnzip = std(timesTillUnzipping)/sqrt(numel(timesTillUnzipping));
