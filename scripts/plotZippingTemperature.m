#dir="/home/other/roald/clusterdata/hairpin/";
dir="/home/other/roald/hairpinMelt/meltingTemp_dt15_measTime300_relaxTime20_Tstart10_Tstep40_nSteps3_salt50/";

N = 10;
filesglob = ([dir, "*N", num2str(N), "*/meltingTemp*"]);

[avgFraction, errFraction, temperatures] = zippingTemperature(filesglob, N) 

errorbar(temperatures, avgFraction, errFraction);
