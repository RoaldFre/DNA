function toNativeData2D(dataDirectory, destinationFile)

data = parse2DRuns(dataDirectory);
save("-binary", destinationFile, "data");

