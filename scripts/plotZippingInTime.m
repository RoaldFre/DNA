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
%file="/home/other/roald/clusterdata/hairpinFormationNoExcl/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12712.itf11";

file = "~/clusterdata/hairpinFormation/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/relaxT90C_zipT10C_unzipT90C_allowUnb0.1_zippedRel100_NXY4/dt15/N80/formation_49977.itf11";

%fileNameBase = "stateInTimeN100_1";
%fileNameBase = "stateInTimeN100_4_noExl";


%XXX Decimate because otherwise gives errors when plotting (data set too large??)
decimateFactor = 3;

[zippingTime,   zippingState,   zippingBound, ...
 relaxTime,     relaxState,     relaxBound, ...
 unzippingTime, unzippingState, unzippingBound] = parseHairpinFormationState(file, decimateFactor);

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
imagesc(unzippingTime * 1e9, [1, N], unzippingState');
%axis([0,180,0,1], 'autoy');
%makeGraphPresentation([fileNameBase,'_unzipping'],presentationDir,xlab,ylab,ylabrule,width,height);

figure

colormap(colmap);
imagesc(zippingTime * 1e9, [1, N], zippingState');
%makeGraphPresentation([fileNameBase,'_zipping'],presentationDir,xlab,ylab,ylabrule,width,height);
