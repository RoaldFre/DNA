function [zippingTime, nucleationTime, zipFromNuclTime, zipIndex, nucleationIndex] = getZippingTimesFromBound(time, numBound, nuclBoundBPs, zippedBoundBPs)

zipIndex = find(numBound > zippedBoundBPs, 1, "first");
nucleationIndex = find(numBound(1:zipIndex) < nuclBoundBPs, 1, "last") + 1;

if isempty(zipIndex) || isempty(nucleationIndex)
	zipIndex
	nucleationIndex
	size(numBound)
	nuclThreshold = [min(numBound), nuclBoundBPs]
	zipThreshold = [max(numBound), zippedBoundBPs]
	error "Couldn't find nucleation or zipping index. Thresholds too strict?"
end

zippingTime = time(zipIndex) - time(1);
nucleationTime = time(nucleationIndex) - time(1);
zipFromNuclTime = zippingTime - nucleationTime;

