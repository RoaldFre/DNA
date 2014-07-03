function eteFreeEnergyToNative(dataFileName)

data = load(dataFileName);

K = data(1,1);
Rref = data(1,2);
time = data(2:end, 1);
ete = data(2:end, 2);

save(dataFileName, "-z", "-binary", "K", "Rref", "time", "ete");
