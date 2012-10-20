function [zippingTime,   zippingState,   zippingBound, ...
          relaxTime,     relaxState,     relaxBound, ...
          unzippingTime, unzippingState, unzippingBound] ...
	  	= parseHairpinFormationState(file, decimateFactor)

addpath('octave-forge');

if (nargin < 2)
	decimateFactor = 1;
end

%[_, tmpfile] = system('mktemp'); %load breaks on this somehow below
tmpfile = "parseHairpinFormationState.temp";

system(['sed -n "s/## \[waiting to zip\] //p" ',file,' > ',tmpfile]);
zipping = load(tmpfile);
zippingTime  = decimate(zipping(:,1), decimateFactor);
zippingBound = decimate(zipping(:,2), decimateFactor);
fullZippingState = zipping(:,3:end);
for i = 1:numel(fullZippingState(1,:))
	zippingState(:,i) = decimate(fullZippingState(:,i), decimateFactor);
end



system(['sed -n "s/## \[zipped relaxation\] //p" ',file,' > ',tmpfile]);
relax = load(tmpfile);
relaxTime  = decimate(relax(:,1), decimateFactor);
relaxBound = decimate(relax(:,2), decimateFactor);
fullRelaxState = relax(:,3:end);
for i = 1:numel(fullZippingState(1,:))
	relaxState(:,i) = decimate(fullRelaxState(:,i), decimateFactor);
end



system(['sed -n "s/## \[waiting to unzip\] //p" ',file,' > ',tmpfile]);
unzipping = load(tmpfile);
unzippingTime  = decimate(unzipping(:,1), decimateFactor);
unzippingBound = decimate(unzipping(:,2), decimateFactor);
fullUnzippingState = unzipping(:,3:end);
for i = 1:numel(fullZippingState(1,:))
	unzippingState(:,i) = decimate(fullUnzippingState(:,i), decimateFactor);
end

system(['rm ',tmpfile]);
