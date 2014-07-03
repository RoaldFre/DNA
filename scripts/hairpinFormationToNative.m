function hairpinFormationToNative(dataFileName, makeLogical)
if nargin < 2
	makeLogical = true;
end

load(dataFileName); % Grab all the variables

% Grab the data in comments

tmpfile = [dataFileName,'.temp'];

system(['sed -n "s/## \[waiting to zip\] //p" ',dataFileName,' > ',tmpfile]);
zipping = load(tmpfile);
zippingTime  = zipping(:,1);
zippingBound = zipping(:,2);
if makeLogical
	zippingState = logical(zipping(:,3:end));
else
	zippingState = single(zipping(:,3:end));
end
clear zipping;

system(['sed -n "s/## \[zipped relaxation\] //p" ',dataFileName,' > ',tmpfile]);
relax = load(tmpfile);
relaxTime  = relax(:,1);
relaxBound = relax(:,2);
if makeLogical
	relaxState = logical(relax(:,3:end));
else
	relaxState = single(relax(:,3:end));
end
clear relax;

system(['sed -n "s/## \[waiting to unzip\] //p" ',dataFileName,' > ',tmpfile]);
unzipping = load(tmpfile);
unzippingTime  = unzipping(:,1);
unzippingBound = unzipping(:,2);
if makeLogical
	unzippingState = logical(unzipping(:,3:end));
else
	unzippingState = single(unzipping(:,3:end));
end
clear unzipping;

system(['rm ',tmpfile]);
clear tmpfile

clear argn ans
save(dataFileName, "-z", "-binary");
% somewhat forced variable name 'dataFileName' to avoid collisions when 
% loading the data file somewhere, because we can't clear dataFileName 
% before we save the variables!
