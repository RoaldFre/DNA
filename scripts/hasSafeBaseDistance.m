% Print 0 if OK, 1 if not OK.
function hasSafeBaseDistance(filename)

data = load(filename);

energies = data(1, 1:end-3);

if sum(energies > 0) == 0
	printf("0\n");
else
	printf("1\n");
end
