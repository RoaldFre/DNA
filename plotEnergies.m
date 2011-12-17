close all;
data = load("/tmp/data.txt");


t  = data(:,1);
E  = data(:,2);
K  = data(:,3);
Vb = data(:,4);
Va = data(:,5);
Vd = data(:,6);
Vs = data(:,7);

% smoothen data by moving average
wndw = ceil(length(t) / 100); %sliding window size
E  = filter(ones(wndw,1)/wndw, 1, E);
K  = filter(ones(wndw,1)/wndw, 1, K);
Vb = filter(ones(wndw,1)/wndw, 1, Vb);
Va = filter(ones(wndw,1)/wndw, 1, Va);
Vd = filter(ones(wndw,1)/wndw, 1, Vd);
Vs = filter(ones(wndw,1)/wndw, 1, Vs);

% cut off the initial bit as it's distorted by the window
t  = t(wndw:end);
E  = E(wndw:end);
K  = K(wndw:end);
Vb = Vb(wndw:end);
Va = Va(wndw:end);
Vd = Vd(wndw:end);
Vs = Vs(wndw:end);


V  = Vb + Va + Vd + Vs;


hold on;
plot(t, E, 'k');
plot(t, K, 'r');
plot(t, V, 'b');
axis([t(1), t(end), 0, 1.1*max(E)]);
hold off;
figure;
hold on;
plot(t, Vb, 'k');
plot(t, Va, 'r');
plot(t, Vd, 'g');
plot(t, Vs, 'b');
axis([t(1), t(end), 0, 1.1 * max([Vb;Va;Vd;Vs])]);
hold off;
