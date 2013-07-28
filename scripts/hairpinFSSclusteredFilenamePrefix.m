function newPrefix = hairpinFSSclusteredFilenamePrefix(resultFilePrefix, clusterSize, dropLogFactor, T)

newPrefix = [resultFilePrefix,'_T',num2str(T),'_clustSize',num2str(clusterSize),'_dropLog',num2str(dropLogFactor)];
