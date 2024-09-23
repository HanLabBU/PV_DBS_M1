%% Create a noisey signal

% Create a noisey signal
Fs = 1000;
T = 1/Fs;
L = 1500;
t = (0:L - 1)*T;

% Sample signal with 50 Hz and 120 Hz sin component
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

% Add some random noise
X = S + 2*randn(size(t));
%X = S;

% Show raw signal
figure;
plot(1000*t, X);
hold on;
title('Signal corrupted with Zero-mean Random Noise');

% Filter each signal and show the bands for each plots
butter_sig = [];
for i = [1:1:200] 
    butter_sig(i, :) = Multi_func.raw_filt(X', i, Fs);
end

% Filter each signal and show the bands for each plots
fir_sig = [];
for i = [1:1:200] 
    fir_sig(i, :) = Multi_func.fir_filt(X', i, Fs);
end


% Butterworth filtered signal
figure;
plot(1000*t, X, 'DisplayName', 'Original');
hold on;
plot(1000*t, butter_sig(50, :), 'DisplayName', '50 Hz filtered');
hold on;
plot(1000*t, butter_sig(80, :), 'DisplayName', '80 Hz filtered');
legend();
xlim([900 1100]);
title('butter signal');

% FIR filtered signal
figure;
plot(1000*t, X);
hold on;
plot(1000*t, fir_sig(50, :), 'DisplayName', '50 Hz filtered');
hold on;
plot(1000*t, fir_sig(80, :), 'DisplayName', '80 Hz filtered');
legend();
xlim([900 1100]);
title('FIR signal');
