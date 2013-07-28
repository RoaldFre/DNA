function resultFile = hairpinFSSclusteredFilename(resultFilePrefix, clusterSize, dropLogFactor, opt)

if nargin != 4
	error "hairpinFSSclusteredFilename got wrong number of arguments"
end

addpath generic

newPrefix = hairpinFSSclusteredFilenamePrefix(resultFilePrefix, clusterSize, dropLogFactor, opt.T);
resultFile = finiteSizeScalingFilename(newPrefix, opt);

