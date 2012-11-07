function [zippingTime,   zippingState,   zippingBound, ...
          relaxTime,     relaxState,     relaxBound, ...
          unzippingTime, unzippingState, unzippingBound] ...
	  	= parseHairpinFormationState(filePath, decimateFactor)

if (nargin < 2)
	decimateFactor = 1;
end

[_, baseDir] = system(['dirname ',filePath]);
[_, fileName] = system(['basename ',filePath]);
baseDir = baseDir(1:end-1); % hacky hacky hack to skip newline at end
fileName = fileName(1:end-1); % hacky hacky hack to skip newline at end

% See if we already parsed and cached this data file
cacheDir = [baseDir,filesep,'parseHairpinCache'];
cacheFile = [cacheDir, filesep, fileName];
[_, err] = stat(cacheFile);

if (err == 0)
	% Found cached file! :-)
	load(cacheFile)
else
	% Not cached yet!
	% Parse data and generate fast cache file
	system(['mkdir -p ',cacheDir]); % make dir if not existent yet
	tmpfile = 'parseHairpinFormationState.temp';

	system(['sed -n "s/## \[waiting to zip\] //p" ',filePath,' > ',tmpfile]);
	zipping = load(tmpfile);
	zippingTime  = zipping(:,1);
	zippingBound = zipping(:,2);
	zippingState = zipping(:,3:end);

	system(['sed -n "s/## \[zipped relaxation\] //p" ',filePath,' > ',tmpfile]);
	relax = load(tmpfile);
	relaxTime  = relax(:,1);
	relaxBound = relax(:,2);
	relaxState = relax(:,3:end);


	system(['sed -n "s/## \[waiting to unzip\] //p" ',filePath,' > ',tmpfile]);
	unzipping = load(tmpfile);
	unzippingTime  = unzipping(:,1);
	unzippingBound = unzipping(:,2);
	unzippingState = unzipping(:,3:end);

	system(['rm ',tmpfile]);

	save(cacheFile, "-z", ...
		"zippingTime",   "zippingState",   "zippingBound", ...
		"relaxTime",     "relaxState",     "relaxBound",   ...
		"unzippingTime", "unzippingState", "unzippingBound");
end

% Decimate
zippingTime  = decimateWrapper(zippingTime,  decimateFactor);
zippingBound = decimateWrapper(zippingBound, decimateFactor);
zippingState = decimateWrapper(zippingState, decimateFactor);

relaxTime  = decimateWrapper(relaxTime,  decimateFactor);
relaxBound = decimateWrapper(relaxBound, decimateFactor);
relaxState = decimateWrapper(relaxState, decimateFactor);

unzippingTime  = decimateWrapper(unzippingTime,  decimateFactor);
unzippingBound = decimateWrapper(unzippingBound, decimateFactor);
unzippingState = decimateWrapper(unzippingState, decimateFactor);
