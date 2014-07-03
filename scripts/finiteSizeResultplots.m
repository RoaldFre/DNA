% DEPRECATED
function finiteSizeResultplots(resultFilePrefix, clusterSize, dropLogFactor, opt, meanStartN, color);

if nargin < 5; meanStartN = 0; end;
if nargin < 6; color = 'b'; end;

resultFile = hairpinFSSclusteredFilename(resultFilePrefix, clusterSize, dropLogFactor, opt)
load(resultFile);

printf("Qualities:\n");
finiteSizeResultplot(clustNs, clustQuals, clustQualErrs, meanStartN, color);
figure
printf("Parameters:\n");
finiteSizeResultplot(clustNs, clustPs, clustPErrs, meanStartN, color);

