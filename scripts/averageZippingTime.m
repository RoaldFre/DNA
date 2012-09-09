function [meanZip, errZip, meanUnzip, errUnzip] = averageZippingTime(filesglob)

[timesTillZipping, timesTillUnzipping, temperature, temperatures, data, timestep, allowedUnboundBPs, allowedBoundBps] = parseHairpinFormation(filesglob, false);

meanZip = mean(timesTillZipping);
errZip = std(timesTillZipping)/sqrt(numel(timesTillZipping) - 1);

meanUnzip = mean(timesTillUnzipping);
errUnzip = std(timesTillUnzipping)/sqrt(numel(timesTillUnzipping) - 1);
