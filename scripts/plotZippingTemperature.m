dir="/home/other/roald/clusterdata/hairpin/";

N = 10;
filesglob = ([dir, "*N", num2str(N), "*/hairpin*"]);
%filesglob = "/home/other/roald/clusterdata/hairpin/N4_dt15_time1450_measTime150_Tsample20_Tstart20_Tstep40_nSteps3/hairpin_wait0_interval10_allowUnb1_9929.itf11"

[avgFraction, errFraction, temperatures] = zippingTemperature(filesglob, N) 

errorbar(temperatures, avgFraction, errFraction);
