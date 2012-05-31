file="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12366.itf11";

%this one 'gets stuck' again:
%file="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12365.itf11";

%file="/home/other/roald/clusterdata/hairpinFormation/formation_zipT20_unzipT180_allowUnb2_allowB2_zippedRel10/dt15_time100000/N100/formation_12366.itf11";


%[_, tmpfile] = system('mktemp'); %load breaks on this somehow below
tmpfile = "plotZippingInTimeTemp";
system(['sed -n "s/## \[waiting to zip\] //p" ',file,' > ',tmpfile]);
zipping = load(tmpfile);
zippingTimes = zipping(:,1) - zipping(1,1);
zippingBound = zipping(:,2);
zippingState = zipping(:,3:end);

system(['sed -n "s/## \[waiting to unzip\] //p" ',file,' > ',tmpfile]);
unzipping = load(tmpfile);
unzippingTimes = unzipping(:,1) - unzipping(1,1);
unzippingBound = unzipping(:,2);
unzippingState = unzipping(:,3:end);

system(['rm ',tmpfile]);

N = numel(zippingState(1,:));

colormap([1 1 1; 0 0 0]);
imagesc(zippingTimes * 1e9, [1, N], zippingState');
xlabel("time (ns)");
ylabel("base pair");
print("zipping.png")

figure
colormap([1 1 1; 0 0 0]);
imagesc(unzippingTimes * 1e9, [1, N], unzippingState');
xlabel("time (ns)");
ylabel("base pair");
print("unzipping.png")



