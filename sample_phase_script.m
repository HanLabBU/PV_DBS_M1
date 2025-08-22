clc;

t = 0:1/1024:1;
x = sin(2*pi*60*t);
y = hilbert(x);

phases = rad2deg(angle(y));

figure;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, x);
xlim([0 0.1]);
title('Original Signal');
nexttile;
plot(t, phases, '.');
xlim([0 0.1]);
ylim([-180 180]);
ylabel('Degrees');
title('Phases in Degrees');

% Determine what are negative radians in the polarhistogram
figure;
phases = [-pi];
polarhistogram(phases);
