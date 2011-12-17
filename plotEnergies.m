function plotEnergies(filename_prefix, caption_prefix)
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

destdir   = 'latex/images';
relImgDir = 'images';
ylabrule  = '-1.5cm';
xlab      = 'time (s)';
ylab      = 'Energy (eV)';
width     = '800';
height    = '600';


hold on;
plot(t, E, 'k;Total;');
plot(t, K, 'r;Kinetic;');
plot(t, V, 'b;Potential;');
axis([t(1), t(end), 0, 1.2*max(E)]);

name      = [filename_prefix,'-EKV'];
caption   = [caption_prefix,' Total, kinetic and potential energy.'];
makeGraph(name,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
hold off;


figure;
hold on;
plot(t, Vb, 'k;Bond;');
plot(t, Va, 'r;Angle;');
plot(t, Vd, 'g;Dihedral;');
plot(t, Vs, 'b;Stacking;');
axis([t(1), t(end), 0, 1.2 * max([Vb;Va;Vd;Vs])]);

name      = [filename_prefix,'-V'];
caption   = [caption_prefix,' Different potential energies.'];
makeGraph(name,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
hold off;
