function toNativeData2D(filesglob, destinationFile)

data = parse2DRuns(filesglob);
nRuns = size(data)(3)
save("-binary", [destinationFile], "data");

