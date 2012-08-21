addpath('octave-forge');


%stays one rather homogenous 'chunk'
%file="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12355.itf11";

%file="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12366.itf11";

%this one 'gets stuck' again:
%file="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12365.itf11";

%file="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12354.itf11";


%NO EXCLUSION
%file="/home/other/roald/clusterdata/hairpinFormationNoExcl/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12856.itf11";
%Really funky and short:
%file="/home/other/roald/clusterdata/hairpinFormationNoExcl/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12914.itf11";
file="/home/other/roald/clusterdata/hairpinFormationNoExcl/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12712.itf11";

%fileNameBase = "stateInTimeN100_1";
fileNameBase = "stateInTimeN100_4_noExl";



%[_, tmpfile] = system('mktemp'); %load breaks on this somehow below
tmpfile = "plotZippingInTimeTemp";
system(['sed -n "s/## \[waiting to zip\] //p" ',file,' > ',tmpfile]);
zipping = load(tmpfile);
%XXX Decimate because otherwise gives errors when plotting (data set too large??)
decimateFactor = 3;
zippingTimes = decimate(zipping(:,1) - zipping(1,1), decimateFactor);
zippingBound = decimate(zipping(:,2), decimateFactor);
fullZippingState = zipping(:,3:end);
clear('zippingState');
for i = 1:numel(fullZippingState(1,:))
	zippingState(:,i) = decimate(fullZippingState(:,i), decimateFactor);
end

system(['sed -n "s/## \[waiting to unzip\] //p" ',file,' > ',tmpfile]);
unzipping = load(tmpfile);
unzippingTimes = unzipping(:,1) - unzipping(1,1);
unzippingBound = unzipping(:,2);
unzippingState = unzipping(:,3:end);

system(['rm ',tmpfile]);

N = numel(zippingState(1,:));

presentationDir = "../presentation/images/";
xlab = "Time (ns)";
ylab = "Base pair";
ylabrule = "1.0cm";
width = '1200';
height = '400';


colmap = [1 1 1; 0 0 0];
%colmap = linspace(1,0,100)' * [1 1 1];

colormap(colmap);
imagesc(unzippingTimes * 1e9, [1, N], unzippingState');
axis([0,180,0,1], 'autoy');
makeGraphPresentation([fileNameBase,'_unzipping'],presentationDir,xlab,ylab,ylabrule,width,height);

%figure

colormap(colmap);
imagesc(zippingTimes * 1e9, [1, N], zippingState');
makeGraphPresentation([fileNameBase,'_zipping'],presentationDir,xlab,ylab,ylabrule,width,height);
