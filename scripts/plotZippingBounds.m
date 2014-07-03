NsTimes = [10 100;
           14 100;
           20 200;
           28 300;
	   40 500;
	  ];
Ns = NsTimes(:,1);
minTimes = NsTimes(:,2);

numNs = numel(Ns);
hold on;
for i = 1 : numNs
	plotZippingBoundNuclThreshold(Ns(i), minTimes(i), 1);
end
hold off;
