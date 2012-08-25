function diffusionToNativeData2D(filesglob, destinationFile)

data = diffusionParse2DRuns(filesglob);
nRuns = size(data)(3)
save("-binary", [destinationFile], "data");

