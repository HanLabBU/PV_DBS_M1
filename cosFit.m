clear all;
close all;
clc;

x = [0:50:1000];
y([1:2:length(x)]) = 1;
x1 = 0:max(x);
y1 = interp1(x, y, x1);

figure;
plot(x1, y1, '.');
hold on;
plot(y, y, 'o');

amp = range(y1);
yz = y1-max(y1)+(amp/2);
zx = x1(yz .* circshift(yz,[0 1]) <= 0);
per = 2*mean(diff(zx));
ymean = mean(y1);

fit = @(b,x) b(1).*(cos(2*pi*x./b(2) + 2*pi/b(3))) + b(4);
fcn = @(b) sum((fit(b,x1) - y1).^2);

s = fminsearch(fcn, [amp; per; 0; ymean]);

figure;
plot(x1, y1, 'r', x1, fit(s, x1), 'b');
