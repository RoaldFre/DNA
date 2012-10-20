function [meanZip, errZip, meanUnzip, errUnzip, meanZipFromNucleation, errZipFromNucleation] = averageZippingTime(filesglob)

[timesTillZipping, timesTillUnzipping, timesTillZippingFromNucleation] = parseHairpinFormation(filesglob);

meanZip = mean(timesTillZipping);
errZip = std(timesTillZipping)/sqrt(numel(timesTillZipping) - 1);

meanZipFromNucleation = mean(timesTillZippingFromNucleation);
errZipFromNucleation = std(timesTillZippingFromNucleation)/sqrt(numel(timesTillZippingFromNucleation) - 1);

meanUnzip = mean(timesTillUnzipping);
errUnzip = std(timesTillUnzipping)/sqrt(numel(timesTillUnzipping) - 1);
