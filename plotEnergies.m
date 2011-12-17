data = load("/tmp/data.txt");

t  = data(:,1);
E  = data(:,2);
K  = data(:,3);
Vb = data(:,4);
Va = data(:,5);
Vd = data(:,6);
Vs = data(:,7);
V  = Vb + Va + Vd + Vs;

hold on;
plot(t, E, 'k');
plot(t, K, 'r');
plot(t, V, 'b');
hold off;
figure;
hold on;
plot(t, Vb, 'k');
plot(t, Va, 'r');
plot(t, Vd, 'g');
plot(t, Vs, 'b');
hold off;
