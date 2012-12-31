% if numZippable < 0, then the *fromNucleation data is not calculated
function [meanZip, errZip, meanUnzip, errUnzip, meanZipFromNucleation, errZipFromNucleation] = averageZippingTime(filesglob, numZippable)

[timesTillZipping, timesTillUnzipping, timesTillZippingFromNucleation] = parseHairpinFormation(filesglob, numZippable);

meanZip = mean(timesTillZipping);
errZip = std(timesTillZipping)/sqrt(numel(timesTillZipping) - 1);

if numZippable >= 0
	meanZipFromNucleation = mean(timesTillZippingFromNucleation);
	errZipFromNucleation = std(timesTillZippingFromNucleation)/sqrt(numel(timesTillZippingFromNucleation) - 1);
end

meanUnzip = mean(timesTillUnzipping);
errUnzip = std(timesTillUnzipping)/sqrt(numel(timesTillUnzipping) - 1);
