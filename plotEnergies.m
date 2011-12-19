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
tns = t * 1e9; %got to nanoseconds

destdir   = 'latex/images';
relImgDir = 'images';
ylabrule  = '-1.5cm';
xlab      = 'time (ns)';
ylab      = 'Energy (eV)';
width     = '1000';
height    = '800';


hold on;
plot(tns, E, 'k;Total;');
plot(tns, K, 'r;Kinetic;');
plot(tns, V, 'b;Potential;');
axis([tns(1), tns(end), 0, 1.2*max(E)]);

name      = [filename_prefix,'-EKV'];
caption   = [caption_prefix,' Total, kinetic and potential energy.'];
makeGraph(name,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
hold off;


figure;
hold on;
plot(tns, Vb, 'k;Bond;');
plot(tns, Va, 'r;Angle;');
plot(tns, Vd, 'g;Dihedral;');
plot(tns, Vs, 'b;Stacking;');
axis([tns(1), tns(end), 0, 1.2 * max([Vb;Va;Vd;Vs])]);

name      = [filename_prefix,'-V'];
caption   = [caption_prefix,' Different potential energies.'];
makeGraph(name,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
hold off;
