dir="/home/other/roald/clusterdata/hairpinFormation/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/relaxT90C_zipT10C_unzipT90C_allowUnb0.1_zippedRel100_NXY4/dt15/";

savedir = [dir, "/nativeCells/"]
system(['mkdir -p ',savedir,' &>/dev/null']);

Ns = [30 40 50 60 80 100];

more off;
for i = 1 : numel(Ns)
	printf("N = %d\n", Ns(i));
	%data = parseHairpinFormationToCell([dir, "N", num2str(Ns(i)), "*/formation*itf11"]);
	%save("-binary", "-z", [savedir, "N", num2str(Ns(i))], "data")

	data = parseHairpinFormationToCellNoState([dir, "N", num2str(Ns(i)), "*/formation*itf11"]);
	save("-binary", "-z", [savedir, "N", num2str(Ns(i)), "_noState"], "data")
end
more on;
